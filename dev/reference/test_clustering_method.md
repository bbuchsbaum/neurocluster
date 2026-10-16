# Test Clustering Method on Synthetic Data

Runs a clustering method on synthetic ground truth data and reports
accuracy metrics.

## Usage

``` r
test_clustering_method(
  method,
  noise_sd = 0,
  n_clusters = 27,
  seed = 42,
  verbose = TRUE,
  ...
)
```

## Arguments

- method:

  Clustering method name (e.g., "supervoxels", "snic", "flash3d")

- noise_sd:

  Noise level for synthetic data. Default 0 (perfect separability).

- n_clusters:

  Number of ground truth clusters (must be perfect cube). Default 27.

- seed:

  Random seed for reproducibility.

- verbose:

  Print results. Default TRUE.

- ...:

  Additional arguments passed to cluster4d()

## Value

A list with:

- metrics:

  Output from clustering_accuracy()

- result:

  Full clustering result object

- data:

  Synthetic data used

- time:

  Elapsed time in seconds

## Examples

``` r
if (FALSE) { # \dontrun{
# Test FLASH-3D with 27 cubic clusters
test_clustering_method("flash3d")

# Test with 8 clusters
test_clustering_method("slic", n_clusters = 8)

# Test with noise
test_clustering_method("snic", noise_sd = 0.5)
} # }
```
