#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
})

build_option_list <- function() {
  list(
    make_option("--seurat-rds", type = "character", help = "Path to processed Seurat .rds file"),
    make_option("--marker-file", type = "character", help = "Path to marker table (.tsv or .csv)"),
    make_option("--cluster-column", type = "character", default = "seurat_clusters",
                help = "Metadata column containing cluster labels [default %default]"),
    make_option("--assay", type = "character", default = NULL,
                help = "Assay to use. Defaults to Seurat DefaultAssay(object)"),
    make_option("--slot", type = "character", default = "data",
                help = "Expression slot to score from [default %default]"),
    make_option("--output-dir", type = "character", default = "results/cluster_annotation",
                help = "Directory for output files [default %default]"),
    make_option("--min-pct-expressing", type = "double", default = 0.1,
                help = "Expression threshold for percent-expressing calculations [default %default]"),
    make_option("--score-weight-expression", type = "double", default = 0.7,
                help = "Weight for average expression component [default %default]"),
    make_option("--score-weight-detection", type = "double", default = 0.3,
                help = "Weight for percent-expressing component [default %default]"),
    make_option("--save-annotated-rds", action = "store_true", default = FALSE,
                help = "Write annotated Seurat object to output directory")
  )
}

stop_if_missing <- function(value, flag) {
  if (is.null(value) || identical(value, "")) {
    stop(sprintf("Missing required argument: %s", flag), call. = FALSE)
  }
}

read_marker_table <- function(path) {
  extension <- tolower(tools::file_ext(path))
  sep <- if (extension == "csv") "," else "\t"
  markers <- read.table(
    path,
    header = TRUE,
    sep = sep,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = ""
  )

  required_columns <- c("cell_type", "gene")
  missing_columns <- setdiff(required_columns, colnames(markers))
  if (length(missing_columns) > 0) {
    stop(
      sprintf(
        "Marker table is missing required columns: %s",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!"weight" %in% colnames(markers)) {
    markers$weight <- 1
  }
  if (!"direction" %in% colnames(markers)) {
    markers$direction <- "positive"
  }

  markers$cell_type <- trimws(markers$cell_type)
  markers$gene <- trimws(markers$gene)
  markers$direction <- tolower(trimws(markers$direction))
  markers$weight <- as.numeric(markers$weight)
  markers$weight[is.na(markers$weight)] <- 1

  invalid_direction <- !markers$direction %in% c("positive", "negative")
  if (any(invalid_direction)) {
    stop("Marker `direction` values must be 'positive' or 'negative'.", call. = FALSE)
  }

  markers <- markers[markers$cell_type != "" & markers$gene != "", , drop = FALSE]
  if (nrow(markers) == 0) {
    stop("Marker table is empty after cleaning.", call. = FALSE)
  }

  markers
}

safe_scale_rows <- function(mat) {
  scaled <- t(apply(mat, 1, function(x) {
    if (all(is.na(x)) || stats::sd(x) == 0) {
      rep(0, length(x))
    } else {
      as.numeric(scale(x))
    }
  }))
  rownames(scaled) <- rownames(mat)
  colnames(scaled) <- colnames(mat)
  scaled
}

compute_cluster_detection <- function(expr_mat, clusters, threshold) {
  cluster_ids <- unique(clusters)
  detection <- sapply(cluster_ids, function(cluster_id) {
    cluster_cells <- names(clusters)[clusters == cluster_id]
    if (length(cluster_cells) == 0) {
      return(rep(0, nrow(expr_mat)))
    }
    rowMeans(expr_mat[, cluster_cells, drop = FALSE] > threshold)
  })

  if (is.null(dim(detection))) {
    detection <- matrix(detection, ncol = 1)
    colnames(detection) <- cluster_ids
  }

  rownames(detection) <- rownames(expr_mat)
  detection
}

score_celltype_for_cluster <- function(marker_subset, cluster_id, expr_z, detect_mat) {
  genes_present <- intersect(marker_subset$gene, rownames(expr_z))
  if (length(genes_present) == 0) {
    return(list(
      score = NA_real_,
      n_markers_total = nrow(marker_subset),
      n_markers_found = 0
    ))
  }

  rows <- match(genes_present, marker_subset$gene)
  matched_markers <- marker_subset[rows, , drop = FALSE]

  expr_component <- expr_z[genes_present, cluster_id]
  detect_component <- detect_mat[genes_present, cluster_id]

  sign_multiplier <- ifelse(matched_markers$direction == "negative", -1, 1)
  weighted_score <- sign_multiplier * matched_markers$weight *
    (opt$score_weight_expression * expr_component +
       opt$score_weight_detection * detect_component)

  list(
    score = sum(weighted_score) / sum(matched_markers$weight),
    n_markers_total = nrow(marker_subset),
    n_markers_found = length(genes_present)
  )
}

assign_cluster_labels <- function(score_table) {
  split_by_cluster <- split(score_table, score_table$cluster)
  assignments <- lapply(split_by_cluster, function(df) {
    ranked <- df[order(df$score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    top_row <- ranked[1, , drop = FALSE]
    second_score <- if (nrow(ranked) > 1) ranked$score[2] else NA_real_
    top_row$score_margin <- top_row$score - second_score
    top_row
  })
  do.call(rbind, assignments)
}

option_parser <- OptionParser(option_list = build_option_list())
opt <- parse_args(option_parser)

stop_if_missing(opt$`seurat-rds`, "--seurat-rds")
stop_if_missing(opt$`marker-file`, "--marker-file")

dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)

message("Loading Seurat object...")
seurat_obj <- readRDS(opt$`seurat-rds`)

if (!inherits(seurat_obj, "Seurat")) {
  stop("Input file is not a Seurat object.", call. = FALSE)
}

if (!(opt$`cluster-column` %in% colnames(seurat_obj@meta.data))) {
  stop(
    sprintf("Cluster column '%s' was not found in Seurat metadata.", opt$`cluster-column`),
    call. = FALSE
  )
}

assay_to_use <- if (is.null(opt$assay)) DefaultAssay(seurat_obj) else opt$assay
if (!(assay_to_use %in% names(seurat_obj@assays))) {
  stop(sprintf("Assay '%s' not found in object.", assay_to_use), call. = FALSE)
}

DefaultAssay(seurat_obj) <- assay_to_use

message("Reading marker table...")
markers <- read_marker_table(opt$`marker-file`)

message("Extracting expression data...")
expr_mat <- GetAssayData(seurat_obj, assay = assay_to_use, slot = opt$slot)
expr_mat <- as.matrix(expr_mat)

clusters <- seurat_obj@meta.data[[opt$`cluster-column`]]
names(clusters) <- colnames(seurat_obj)
clusters <- as.character(clusters)
names(clusters) <- colnames(seurat_obj)
cluster_ids <- unique(clusters)

expr_cluster_avg <- AverageExpression(
  seurat_obj,
  assays = assay_to_use,
  group.by = opt$`cluster-column`,
  slot = opt$slot,
  verbose = FALSE
)[[assay_to_use]]
expr_cluster_avg <- as.matrix(expr_cluster_avg)
expr_cluster_avg <- expr_cluster_avg[, cluster_ids, drop = FALSE]

message("Computing cluster-level marker summaries...")
expr_z <- safe_scale_rows(expr_cluster_avg)
detect_mat <- compute_cluster_detection(expr_mat, clusters, opt$`min-pct-expressing`)
detect_mat <- detect_mat[, cluster_ids, drop = FALSE]

marker_coverage <- unique(markers[, c("cell_type", "gene", "direction", "weight")])
marker_coverage$present_in_object <- marker_coverage$gene %in% rownames(expr_mat)

cell_types <- unique(markers$cell_type)
score_rows <- lapply(cluster_ids, function(cluster_id) {
  lapply(cell_types, function(cell_type) {
    marker_subset <- markers[markers$cell_type == cell_type, , drop = FALSE]
    score_result <- score_celltype_for_cluster(marker_subset, cluster_id, expr_z, detect_mat)
    data.frame(
      cluster = cluster_id,
      cell_type = cell_type,
      score = score_result$score,
      n_markers_total = score_result$n_markers_total,
      n_markers_found = score_result$n_markers_found,
      marker_fraction_found = score_result$n_markers_found / score_result$n_markers_total,
      stringsAsFactors = FALSE
    )
  })
})

score_table <- do.call(rbind, unlist(score_rows, recursive = FALSE))
assignment_table <- assign_cluster_labels(score_table)
colnames(assignment_table)[colnames(assignment_table) == "cell_type"] <- "predicted_cell_type"
assignment_table$confidence <- with(
  assignment_table,
  ifelse(is.na(score_margin), NA_real_, plogis(score_margin))
)

cluster_annotation_map <- setNames(
  assignment_table$predicted_cell_type,
  assignment_table$cluster
)

seurat_obj$predicted_cell_type <- unname(cluster_annotation_map[clusters])
seurat_obj$predicted_cell_type_source_cluster <- clusters

write.table(
  score_table,
  file = file.path(opt$`output-dir`, "cluster_annotation_scores.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  assignment_table,
  file = file.path(opt$`output-dir`, "cluster_annotations.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  marker_coverage,
  file = file.path(opt$`output-dir`, "marker_coverage.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

if (isTRUE(opt$`save-annotated-rds`)) {
  saveRDS(seurat_obj, file = file.path(opt$`output-dir`, "annotated_seurat.rds"))
}

message("Cluster annotation complete.")
