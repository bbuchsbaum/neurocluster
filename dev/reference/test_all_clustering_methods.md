# Test All Clustering Methods on Synthetic Data

Comprehensive test of all available clustering methods on synthetic
ground truth data with cubic clusters.

## Usage

``` r
test_all_clustering_methods(
  noise_levels = c(0, 0.1, 0.5),
  methods = c("flash3d", "corr_slic", "snic", "slic", "slice_msf", "supervoxels", "g3s",
    "rena"),
  n_clusters = 27,
  seed = 42
)
```

## Arguments

- noise_levels:

  Vector of noise standard deviations to test. Default c(0, 0.1, 0.5)
  tests perfect, low noise, and moderate noise.

- methods:

  Vector of method names to test. Default tests all main methods.

- n_clusters:

  Number of ground truth clusters (must be perfect cube). Default 27.

- seed:

  Random seed for reproducibility.

## Value

A data.frame with columns: method, noise_sd, ari, nmi, accuracy,
n_clusters, time, perfect

## Examples

``` r
if (FALSE) { # \dontrun{
# Test all methods with 27 clusters
results <- test_all_clustering_methods()
print(results)

# Test with 8 clusters (simpler)
results <- test_all_clustering_methods(n_clusters = 8)

# Test specific methods
results <- test_all_clustering_methods(
  methods = c("flash3d", "corr_slic", "snic"),
  noise_levels = c(0, 0.2)
)
} # }
```
