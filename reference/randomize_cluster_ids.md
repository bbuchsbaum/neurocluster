# Randomize Cluster IDs

Randomly permute cluster IDs to break any spatial ordering. This is
useful when visualizing clusters with continuous color scales, as
spatially-ordered cluster IDs create artificial gradients that can be
misleading.

## Usage

``` r
randomize_cluster_ids(cluster_result, seed = NULL)
```

## Arguments

- cluster_result:

  A clustering result object (e.g., from cluster4d, supervoxels_flash3d,
  slice_msf, etc.) containing at minimum:

  - `cluster`: vector of cluster assignments

  - `clusvol`: ClusteredNeuroVol object

  - Optionally: `centers`, `coord_centers`, etc.

- seed:

  Optional random seed for reproducibility. If NULL, uses R's current
  random state. Default is NULL.

## Value

A modified cluster_result object with randomized cluster IDs. All
components (cluster vector, clusvol, centers, coord_centers, etc.) are
updated consistently.

## Details

Many clustering algorithms assign cluster IDs sequentially as clusters
are encountered during processing. When voxels are processed in memory
order (which is typically z-major in 3D arrays), this creates systematic
ordering where cluster ID correlates with z-position (inferior to
superior in brain imaging).

This spatial ordering creates two problems:

1.  **Visualization artifacts**: When viewed with continuous color
    scales (e.g., viridis, rainbow), ordered cluster IDs create
    artificial gradients that suggest spatial structure where none
    exists. Randomized IDs show each cluster as a distinct color patch.

2.  **Interpretation bias**: Ordered IDs may mislead users into thinking
    cluster numbering has anatomical significance (e.g., "cluster 1 is
    always inferior, cluster 100 is always superior").

This function breaks the ordering by randomly permuting cluster IDs
while maintaining all cluster memberships and properties.

## Examples

``` r
if (FALSE) { # \dontrun{
# Standard workflow with ordered IDs
result <- cluster4d(vec, mask, K = 100, method = "flash3d")

# Randomize for better visualization
result_random <- randomize_cluster_ids(result, seed = 42)

# Compare visualizations
# Ordered IDs show gradient from inferior to superior
plot(result$clusvol)

# Randomized IDs show distinct color patches
plot(result_random$clusvol)

# Works with any clustering method
result_msf <- cluster4d(vec, mask, K = 100, method = "slice_msf")
result_msf_random <- randomize_cluster_ids(result_msf)
} # }
```
