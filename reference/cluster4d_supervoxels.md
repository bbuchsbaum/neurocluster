# Cluster4d using supervoxels method

Wrapper for supervoxels algorithm with standardized interface.

## Usage

``` r
cluster4d_supervoxels(
  vec,
  mask,
  n_clusters = 100,
  spatial_weight = 0.5,
  max_iterations = 50,
  connectivity = 27,
  parallel = TRUE,
  verbose = FALSE,
  sigma1 = 1,
  sigma2 = 2.5,
  use_gradient = TRUE,
  converge_thresh = 0.001,
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

  Method-specific spatial control in `[0, 1]`. For methods that define a
  convex feature/spatial blend, zero disables the spatial term and one
  disables the feature term. Other methods use the value as a bounded
  mapping to a native scale; their endpoints do not have that blend
  meaning. Default 0.5. This parameter is inactive for `"rena"` and
  `"rena_plus"`; supplying it explicitly for either method is an error.
  Maps to method-specific parameters:

  - supervoxels: `alpha = 1 - spatial_weight` (0 = all spatial, 1 = all
    feature)

  - snic: `compactness = spatial_weight * 10` (range 0-10)

  - slic: `compactness = spatial_weight * 20` (typical range 1-20)

  - corr_slic/brs_slic: direct convex blend between correlation and
    scaled spatial distance

  - slice_msf: `compactness = spatial_weight * 10`, followed by FH
    component scale `2 / (compactness + 1)`. This changes the
    graph-segmentation scale; it is not a convex spatial/feature blend.

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

- parallel:

  Enable parallel processing for supervoxels and slic. For all other
  methods, an explicit `TRUE` is rejected; omitted or `FALSE` selects
  the serial implementation. Default TRUE for supported methods.

- verbose:

  Print progress messages. Default FALSE.

- sigma1:

  Bandwidth of heat kernel for features (supervoxels-specific)

- sigma2:

  Bandwidth of heat kernel for coordinates (supervoxels-specific)

- use_gradient:

  Use gradient-based initialization

- converge_thresh:

  Convergence threshold

- ...:

  Additional parameters passed to supervoxels

## Value

A cluster4d_result object
