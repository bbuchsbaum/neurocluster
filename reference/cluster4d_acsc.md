# Cluster4d using ACSC method

Wrapper for Adaptive Correlation Superclustering algorithm with
standardized interface. ACSC uses graph-based clustering with Louvain
community detection and optional boundary refinement.

## Usage

``` r
cluster4d_acsc(
  vec,
  mask,
  n_clusters = 100,
  spatial_weight = 0.5,
  max_iterations = 5,
  verbose = FALSE,
  block_size = 2,
  correlation_metric = c("pearson", "spearman", "robust"),
  spatial_weighting = c("gaussian", "binary"),
  refine = TRUE,
  max_refine_iter = NULL,
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

- verbose:

  Print progress messages. Default FALSE.

- block_size:

  Approximate side length of initial blocks. Default 2.

- correlation_metric:

  Correlation definition used consistently for graph construction and
  refinement: `"pearson"`, `"spearman"`, or `"robust"`.

- spatial_weighting:

  Spatial block-edge weights: `"gaussian"` or `"binary"`.

- refine:

  Logical; whether to refine boundaries. Default TRUE.

- max_refine_iter:

  Maximum iterations for boundary refinement. Default 5.

- ...:

  Additional parameters passed to acsc

## Value

A cluster4d_result object
