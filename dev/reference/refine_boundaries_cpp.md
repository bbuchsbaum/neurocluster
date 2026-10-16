# Refine boundary voxels using C++ acceleration

Refine boundary voxels using C++ acceleration

## Usage

``` r
refine_boundaries_cpp(
  voxel_labels,
  feature_mat_normalized,
  neighbor_indices,
  boundary_voxels,
  max_iter = 5L
)
```

## Arguments

- voxel_labels:

  Integer vector of current cluster assignments (0-based internally)

- feature_mat_normalized:

  Numeric matrix (voxels x time) with NORMALIZED features

- neighbor_indices:

  Integer matrix (voxels x K) of nearest neighbor indices (1-based from
  R)

- boundary_voxels:

  Integer vector of boundary voxel indices (1-based from R)

- max_iter:

  Maximum number of refinement iterations

## Value

List containing updated labels and number of iterations performed
