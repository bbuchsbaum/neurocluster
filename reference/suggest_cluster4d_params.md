# Suggest cluster4d parameters based on data characteristics

Provides parameter recommendations based on data size and user
priorities. Experimental methods are reported separately and are never
returned as the recommended method.

## Usage

``` r
suggest_cluster4d_params(
  n_voxels,
  n_timepoints,
  priority = c("balanced", "speed", "quality", "memory")
)
```

## Arguments

- n_voxels:

  Number of voxels in mask

- n_timepoints:

  Number of time points

- priority:

  What to optimize for: "speed", "quality", "memory", or "balanced"

## Value

A list with suggested parameters for supported methods plus
`experimental_methods` and an `experimental_note`. Experimental methods
are not recommended or assigned suggested parameter sets.

## Examples

``` r
# Get parameter suggestions for a typical fMRI dataset
params <- suggest_cluster4d_params(50000, 200, priority = "balanced")
print(params$recommended_method)
#> [1] "flash3d"
print(params$n_clusters)
#> [1] 200

# Speed-optimized parameters for large dataset
speed_params <- suggest_cluster4d_params(100000, 300, priority = "speed")
print(speed_params$recommended_method)
#> [1] "flash3d"

# Quality-optimized parameters for smaller dataset
quality_params <- suggest_cluster4d_params(10000, 150, priority = "quality")
print(quality_params$n_clusters)
#> [1] 48

# Memory-efficient parameters
memory_params <- suggest_cluster4d_params(200000, 100, priority = "memory")
print(memory_params$recommended_method)
#> [1] "snic"
```
