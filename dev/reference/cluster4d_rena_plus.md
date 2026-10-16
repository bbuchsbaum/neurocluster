# Cluster4d wrapper for ReNA++

Cluster4d wrapper for ReNA++

## Usage

``` r
cluster4d_rena_plus(
  vec,
  mask,
  n_clusters = 100,
  r = 2,
  lambda = 1,
  grad_img = NULL,
  connectivity = 26,
  max_iterations = 50,
  verbose = FALSE,
  ...
)
```

## Arguments

- vec:

  A `NeuroVec` instance supplying the 4D data (x, y, z, time) to cluster

- mask:

  A `NeuroVol` mask defining the voxels to include in clustering. Values
  must be finite. Voxels are included exactly when the mask value is
  strictly greater than zero; zero and negative values are excluded.

- n_clusters:

  Exact finite integer target number of clusters, from one through the
  number of included mask voxels.

- r:

  Over-clustering factor for coarsening stage.

- lambda:

  Gradient penalty weight.

- grad_img:

  Optional gradient/intensity vector (see `rena_plus`).

- connectivity:

  Integer neighborhood connectivity. Supported values are
  method-specific: supervoxels accepts 6, 18, 26, or 27; slic,
  corr_slic, brs_slic, and slice_msf accept 6 or 26; G3S, ReNA variants,
  and MCL accept 6, 18, or 26. It is inactive for snic, flash3d, acsc,
  and commute.

- max_iterations:

  Positive finite integer iteration limit for iterative methods. It is
  inactive for `"snic"`, `"slice_msf"`, and `"commute"`; supplying it
  explicitly for those methods is an error.

- verbose:

  Print progress messages. Default FALSE.

- ...:

  Additional arguments passed to
  [`rena_plus()`](https://bbuchsbaum.github.io/neurocluster/dev/reference/rena_plus.md).
  ReNA++ has no spatial/feature trade-off, so `spatial_weight` is not a
  formal argument; supplying it through `...` is an error.

## Value

A `cluster4d_result`.
