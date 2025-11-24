# ReNA++: Edge-aware Reciprocal Multi-Level ReNA with Ward refinement

Two-stage pipeline:

1.  Edge-aware reciprocal-nearest-neighbor ReNA coarsening to K' = r \*
    K super-voxels.

2.  Spatially-constrained Ward refinement on the super-graph to reach
    exactly K clusters.

## Usage

``` r
rena_plus(
  bvec,
  mask,
  K = 100,
  r = 5,
  lambda = 1,
  grad_img = NULL,
  connectivity = 26,
  max_iterations = 50,
  verbose = FALSE,
  ...
)
```

## Arguments

- bvec:

  A `NeuroVec` providing voxel time-series.

- mask:

  A `NeuroVol` mask; non-zero voxels are included.

- K:

  Target number of clusters.

- r:

  Over-clustering factor; coarsening stops at K' = ceiling(r \* K).
  Default 5.

- lambda:

  Gradient penalty strength (\>=0). 0 disables edge weighting.

- grad_img:

  Optional numeric vector: gradient/intensity per voxel. Either length =
  prod(dim(mask)) (in which case values are subset by mask) or length =
  number of masked voxels. If NULL, no gradient.

- connectivity:

  Neighborhood connectivity (6, 18, or 26). Default 26.

- max_iterations:

  Max iterations for coarsening stage.

- verbose:

  Logical for progress messages.

- ...:

  Reserved for future options.

## Value

A `cluster4d_result` with additional class `rena_plus_cluster_result`.
