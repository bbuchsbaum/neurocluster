# Cluster4d using correlation-embedded SLIC method

Wrapper for a correlation-first SLIC variant that builds a compact
random projection embedding of voxel time series, then runs local 3D
SLIC updates.

## Usage

``` r
cluster4d_corrslic(
  vec,
  mask,
  n_clusters = 100,
  spatial_weight = 0.5,
  max_iterations = 6,
  connectivity = 6,
  parallel = FALSE,
  verbose = FALSE,
  embedding_dim = 64L,
  adaptive_embedding = FALSE,
  embedding_basis = c("hash", "dct"),
  embedding_whiten = FALSE,
  sketch_repeats = 1L,
  alpha = NULL,
  assign_stride = 1L,
  quantize_assign = FALSE,
  refine_exact_iters = 0L,
  refine_boundary_only = TRUE,
  refine_stride = 1L,
  refine_alpha = NULL,
  seed = 1L,
  min_size = NULL,
  ...,
  .input_contract = NULL
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

- parallel:

  Enable parallel processing for supervoxels and slic. For all other
  methods, an explicit `TRUE` is rejected; omitted or `FALSE` selects
  the serial implementation. Default TRUE for supported methods.

- verbose:

  Print progress messages. Default FALSE.

- embedding_dim:

  Embedding dimension used to approximate correlations. Use `"auto"` (or
  set `adaptive_embedding = TRUE`) for data-adaptive selection.

- adaptive_embedding:

  Logical; if TRUE, choose embedding dimension from data size and time
  length.

- embedding_basis:

  Embedding basis, either `"hash"` (CountSketch-like) or `"dct"`
  (demeaned DCT basis).

- embedding_whiten:

  Logical; if TRUE, whiten embedding dimensions across voxels.

- sketch_repeats:

  Number of independent hash sketches combined per voxel.

- alpha:

  Optional backward-compatible alias for the spatial blend weight. It
  must be in `[0, 1]`: 0 is correlation-only and 1 is spatial-only.

- assign_stride:

  Assignment subsampling stride for coarse updates. 1 disables
  subsampling; values \> 1 process one z-slice phase per iteration
  followed by a final full assignment pass.

- quantize_assign:

  Logical; if TRUE, quantize embeddings and centroid features to int8
  during assignment for faster dot-product distance checks.

- refine_exact_iters:

  Number of exact-correlation refinement passes after coarse
  embedding-SLIC iterations.

- refine_boundary_only:

  Logical; if TRUE, exact refinement is restricted to boundary voxels
  each pass.

- refine_stride:

  Time-axis subsampling stride for exact refinement. Values \> 1 speed
  up refinement by computing correlations on every `refine_stride`-th
  time point (approximate).

- refine_alpha:

  Optional spatial blend weight for exact refinement. If `NULL`, uses
  `alpha`; both endpoints are supported.

- seed:

  Random seed for embedding hash and seed initialization.

- min_size:

  Minimum component size for connectivity enforcement.

- ...:

  Additional arguments (currently unused; reserved for forward
  compatibility).

- .input_contract:

  Internal prevalidated input receipt used by
  [`cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md);
  direct callers should leave this as `NULL`.

## Value

A cluster4d_result object

## Details

`spatial_weight` (or its alias `alpha`) is a convex blend:
`(1 - w) * correlation_distance + w * scaled_spatial_distance`. Coarse
correlation is approximated by a hash or DCT sketch. Assignment may
additionally use int8 quantization or z-plane subsampling; exact
refinement may use a temporal stride. These choices and their stage are
recorded under `metadata$approximation`. The public `centers` are always
recomputed from the final labels in the original time-series space;
sketch centers are separately typed metadata and are never substituted
for them.
