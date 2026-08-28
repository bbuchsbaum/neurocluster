# Cluster4d using slice_msf method

Wrapper for slice_msf algorithm with standardized interface. Slice-MSF
is experimental and is deliberately excluded from
[`suggest_cluster4d_params()`](https://bbuchsbaum.github.io/neurocluster/reference/suggest_cluster4d_params.md)
recommendations pending release recertification.

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
  gamma = 0,
  z_mult = 0,
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

  Reserved compatibility parameter; must be zero.

- z_mult:

  Axial DCT-sketch smoothing fraction in `[0, 1]`; zero leaves sketches
  unsmoothed.

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

A cluster4d_result object with natural/exact-K and effective-control
provenance in `metadata$exact_k_repair`; see
[`slice_msf()`](https://bbuchsbaum.github.io/neurocluster/reference/slice_msf.md).

## Details

Unlike direct
[`slice_msf()`](https://bbuchsbaum.github.io/neurocluster/reference/slice_msf.md),
this wrapper always requests exactly `n_clusters`, defaults to one
non-consensus run, and chooses a conservative automatic `min_size`.
`spatial_weight` maps to `compactness = spatial_weight * 10` and then to
FH scale `2 / (compactness + 1)`; it is not a convex spatial/feature
blend. Positive `gamma` values fail with class
`slice_msf_unsupported_gamma`.
