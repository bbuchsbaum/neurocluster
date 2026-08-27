# Cluster4d using slice_msf method

Wrapper for slice_msf algorithm with standardized interface.

## Usage

``` r
cluster4d_slice_msf(
  vec,
  mask,
  n_clusters = 100,
  spatial_weight = 0.2,
  connectivity = 26,
  parallel = TRUE,
  verbose = FALSE,
  num_runs = 1,
  consensus = FALSE,
  stitch_z = TRUE,
  r = 12,
  gamma = 1,
  z_mult = 0.5,
  min_size = NULL,
  seed = 1L,
  ensemble_fraction = 0.8,
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

- num_runs:

  Number of independent runs

- consensus:

  Use consensus fusion

- stitch_z:

  Whether graph edges may cross axial slices.

- r:

  Number of non-DC DCT modes.

- gamma:

  Reliability exponent used in feature edge distances.

- z_mult:

  Optional axial sketch-smoothing fraction in `[0, 1]`.

- min_size:

  Minimum component size after connectivity enforcement. Default NULL
  (auto).

- seed:

  Integer seed for multi-run DCT subspaces.

- ensemble_fraction:

  Candidate-pool fraction for multi-run DCT subspaces.

- ...:

  Additional parameters passed to slice_msf

## Value

A cluster4d_result object
