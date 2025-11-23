# Cluster4d using slice_msf method

Wrapper for slice_msf algorithm with standardized interface.

## Usage

``` r
cluster4d_slice_msf(
  vec,
  mask,
  n_clusters = 100,
  spatial_weight = 0.5,
  connectivity = 8,
  parallel = TRUE,
  verbose = FALSE,
  num_runs = 3,
  consensus = TRUE,
  stitch_z = TRUE,
  theta_link = 0.85,
  min_contact = 1,
  r = 12,
  gamma = 1.5,
  ...
)
```

## Arguments

- vec:

  A `NeuroVec` instance supplying the 4D data (x, y, z, time) to cluster

- mask:

  A `NeuroVol` mask defining the voxels to include in clustering. If
  numeric, nonzero values define included voxels. If logical, TRUE
  values define included voxels.

- n_clusters:

  Target number of clusters (default 100). Note that some methods may
  produce slightly different numbers of clusters due to algorithmic
  constraints.

- spatial_weight:

  Balance between spatial and feature similarity (0-1). Higher values
  emphasize spatial compactness. Default 0.5. Maps to method-specific
  parameters:

  - supervoxels: `alpha = 1 - spatial_weight` (0 = all spatial, 1 = all
    feature)

  - snic/slic: `compactness = spatial_weight * 20` (typical range 1-20)

  - slice_msf: `compactness = spatial_weight * 10` (typical range 1-10)

  - flash3d: `lambda_s = spatial_weight` (direct mapping)

- connectivity:

  Neighborhood connectivity (6 or 26). Default 26. 6 = face neighbors
  only, 26 = face + edge + corner neighbors.

- parallel:

  Enable parallel processing where supported. Default TRUE.

- verbose:

  Print progress messages. Default FALSE.

- num_runs:

  Number of independent runs

- consensus:

  Use consensus fusion

- stitch_z:

  Stitch clusters across z-slices

- ...:

  Additional parameters passed to slice_msf

## Value

A cluster4d_result object
