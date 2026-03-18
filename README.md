# Cluster annotation pipeline for scRNA-seq

This workspace contains an R-based pipeline to annotate scRNA-seq clusters from a processed Seurat object using a target marker list.

## Files

- `scripts/annotate_clusters_from_markers.R`: main pipeline
- `markers/marker_template.tsv`: example marker table

## Expected inputs

1. A processed Seurat object saved as `.rds`
2. A marker table in `.tsv` or `.csv` format with at least:
   - `cell_type`
   - `gene`

Optional marker columns:

- `weight`: numeric importance weight for the marker. Defaults to `1`.
- `direction`: `positive` or `negative`. Defaults to `positive`.

## What the pipeline does

- Loads a processed Seurat object
- Uses a chosen cluster column from metadata
- Pulls normalized expression from the selected assay and slot
- Scores each cluster against the supplied marker sets
- Assigns the best-matching cell type to each cluster
- Writes summary tables and optionally a new annotated Seurat object

## Example usage

```bash
Rscript scripts/annotate_clusters_from_markers.R \
  --seurat-rds data/processed_seurat.rds \
  --marker-file markers/marker_template.tsv \
  --cluster-column seurat_clusters \
  --output-dir results/cluster_annotation
```

To also save an annotated Seurat object:

```bash
Rscript scripts/annotate_clusters_from_markers.R \
  --seurat-rds data/processed_seurat.rds \
  --marker-file markers/marker_template.tsv \
  --cluster-column seurat_clusters \
  --output-dir results/cluster_annotation \
  --save-annotated-rds
```

## Outputs

- `cluster_annotation_scores.tsv`: score matrix for every cluster and candidate cell type
- `cluster_annotations.tsv`: best label per cluster with confidence metrics
- `marker_coverage.tsv`: which markers were found in the object
- `annotated_seurat.rds`: optional updated object with cluster and cell-level annotations

## Notes

- The script expects the Seurat object to already be normalized.
- Cluster scoring combines average marker expression and percent of cells expressing each marker.
- If your cluster metadata column is not `seurat_clusters`, pass the correct column name with `--cluster-column`.
