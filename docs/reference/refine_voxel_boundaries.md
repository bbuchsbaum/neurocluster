# Refine voxel boundaries using C++ or R implementation

For each boundary voxel, compare correlation with each neighboring
cluster's cached centroid. This approach is much faster than comparing
against all voxel time-series.

## Usage

``` r
refine_voxel_boundaries(
  voxel_labels,
  feature_mat,
  coords,
  max_iter,
  use_cpp = TRUE
)
```

## Arguments

- use_cpp:

  Logical; if TRUE, use C++ accelerated version (default)
