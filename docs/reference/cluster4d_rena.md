# Cluster4d using ReNA method

Wrapper for ReNA algorithm with standardized cluster4d interface.

## Usage

``` r
cluster4d_rena(
  vec,
  mask,
  n_clusters = 100,
  spatial_weight = 0.5,
  max_iterations = 50,
  connectivity = 26,
  verbose = FALSE,
  exact_k = TRUE,
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

- max_iterations:

  Maximum iterations for iterative methods. Default 10. Maps to:
  `iterations` (supervoxels), `max_iter` (snic/slic), `rounds`
  (flash3d).

- connectivity:

  Neighborhood connectivity (6 or 26). Default 26. 6 = face neighbors
  only, 26 = face + edge + corner neighbors.

- verbose:

  Print progress messages. Default FALSE.

- exact_k:

  Logical; if TRUE, apply edge pruning to target exactly `n_clusters`.

- ...:

  Additional parameters (currently unused for ReNA)

## Value

A cluster4d_result object
