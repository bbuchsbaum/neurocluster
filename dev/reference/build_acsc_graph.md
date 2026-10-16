# Build ACSC adjacency graph

Build ACSC adjacency graph

## Usage

``` r
build_acsc_graph(
  block_summary,
  ann_k,
  alpha,
  spatial_weighting,
  block_size,
  correlation_metric = c("pearson", "spearman", "robust"),
  knn_proj_dim = 0L,
  knn_proj_method = c("none", "dct", "rp"),
  knn_proj_seed = 1L
)
```
