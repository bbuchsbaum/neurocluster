# Cluster4d using ReNA method

Wrapper for ReNA algorithm with standardized cluster4d interface.

## Usage

``` r
cluster4d_rena(
  vec,
  mask,
  n_clusters = 100,
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

  A `NeuroVol` mask defining the voxels to include in clustering. Values
  must be finite. Voxels are included exactly when the mask value is
  strictly greater than zero; zero and negative values are excluded.

- n_clusters:

  Finite integer target number of clusters (default 100), from one
  through the number of included mask voxels. Note that some methods may
  produce slightly different numbers of clusters due to algorithmic
  constraints.

- max_iterations:

  Positive finite integer iteration limit for iterative methods. It is
  inactive for `"snic"`, `"slice_msf"`, and `"commute"`; supplying it
  explicitly for those methods is an error.

- connectivity:

  Integer neighborhood connectivity. Supported values are
  method-specific: supervoxels accepts 6, 18, 26, or 27; slic,
  corr_slic, brs_slic, and slice_msf accept 6 or 26; G3S, ReNA variants,
  and MCL accept 6, 18, or 26. It is inactive for snic, flash3d, acsc,
  and commute.

- verbose:

  Print progress messages. Default FALSE.

- exact_k:

  Logical; if TRUE, apply edge pruning to target exactly `n_clusters`.

- ...:

  Additional parameters (currently unused for ReNA). ReNA has no
  spatial/feature trade-off, so `spatial_weight` is not a formal
  argument; supplying it through `...` is an error.

## Value

A cluster4d_result object
