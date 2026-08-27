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

  A `NeuroVol` mask defining the voxels to include in clustering. Values
  must be finite. Voxels are included exactly when the mask value is
  strictly greater than zero; zero and negative values are excluded.

- n_clusters:

  Finite integer target number of clusters (default 100), from one
  through the number of included mask voxels. Note that some methods may
  produce slightly different numbers of clusters due to algorithmic
  constraints.

- spatial_weight:

  Balance between spatial and feature similarity (0-1). Both endpoints
  are supported: 0 disables the spatial term and 1 disables the feature
  term. Higher values emphasize spatial compactness. Default 0.5. This
  parameter is inactive for `"rena"` and `"rena_plus"`; supplying it
  explicitly for either method is an error. Maps to method-specific
  parameters:

  - supervoxels: `alpha = 1 - spatial_weight` (0 = all spatial, 1 = all
    feature)

  - snic: `compactness = spatial_weight * 10` (range 0-10)

  - slic: `compactness = spatial_weight * 20` (typical range 1-20)

  - corr_slic/brs_slic: direct convex blend between correlation and
    scaled spatial distance

  - slice_msf: `compactness = spatial_weight * 10` (typical range 1-10)

  - flash3d: `lambda_s = spatial_weight` (direct mapping)

  - mcl: direct blend between feature and spatial edge similarities

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

  Additional parameters (currently unused for ReNA)

## Value

A cluster4d_result object
