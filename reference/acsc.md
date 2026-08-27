# Adaptive Correlation Superclustering (ACSC)

Clusters fMRI voxels into spatially-coherent groups based on temporal
correlation and spatial proximity. Includes optional refinement for
boundary corrections. The algorithm supports parallel processing via the
future package for improved performance.

## Usage

``` r
acsc(
  bvec,
  mask,
  block_size = 2,
  ann_k = 10,
  alpha = 0.5,
  correlation_metric = c("pearson", "spearman", "robust"),
  spatial_weighting = c("gaussian", "binary"),
  refine = TRUE,
  max_refine_iter = 5,
  K = NULL,
  knn_proj_dim = 0L,
  knn_proj_method = c("none", "dct", "rp"),
  knn_proj_seed = 1L
)
```

## Arguments

- bvec:

  A NeuroVec-like object containing 4D fMRI data.

- mask:

  A NeuroVol-like object (logical or numeric mask).

- block_size:

  Approximate side length of blocks (e.g., 2 or 3). Must be \> 0.

- ann_k:

  Number of approximate (or exact) nearest neighbors per block. Must be
  \>= 1.

- alpha:

  Feature-similarity weight in `[0, 1]`. Zero is spatial-only and one is
  feature-only; both endpoints are supported.

- correlation_metric:

  Similarity definition: Pearson correlation, Spearman correlation of
  average ranks, or robust biweight midcorrelation with tuning constant
  9.

- spatial_weighting:

  Spatial adjacency weighting ("gaussian", "binary").

- refine:

  Logical; whether to refine boundaries.

- max_refine_iter:

  Maximum iterations for boundary refinement. Must be \>= 0.

- K:

  (Optional) Desired number of clusters.

- knn_proj_dim:

  Optional integer; if \> 0, project block summaries to this dimension
  before kNN graph construction (speeds up
  [`FNN::get.knn()`](https://rdrr.io/pkg/FNN/man/get.knn.html) for long
  time-series).

- knn_proj_method:

  Projection method for `knn_proj_dim`: one of `"none"`, `"dct"`, or
  `"rp"` (random projection).

- knn_proj_seed:

  Integer seed used when `knn_proj_method="rp"` (RNG state is restored
  after graph construction).

## Value

A standardized `cluster4d_result` (also classed `acsc_result`) with
final contiguous labels, original-space centers, physical coordinate
centers, provenance, and these ACSC-specific elements:

- cluster_map:

  3D array with cluster labels per voxel.

- graph:

  An `igraph` object used for clustering.

- init_block_label:

  Initial coarse partition (3D array) matching `mask` dimensions.

## Details

ACSC builds feature-neighbor edges from the declared correlation metric
and unions them with six-neighbor block-grid edges. Edge correlations
are mapped to non-negative Louvain similarities as
`(correlation + 1) / 2`. At `alpha = 0`, feature-neighbor selection is
bypassed entirely, so the graph topology and weights depend only on
spatial adjacency. Optional DCT or random projection changes
feature-neighbor candidate selection, but edge weights are always
recomputed using the declared metric. Louvain calls use a fixed internal
seed and restore the caller's RNG state, so identical ACSC inputs
produce the same partition independently of the surrounding random
stream.

Boundary refinement uses the same metric representation as graph
building. Pearson uses centered unit vectors, Spearman uses centered
average ranks, and robust mode uses median-centered Tukey-bisquare
residuals. The C++ path accelerates centroid dot products and recomputes
the active boundary after every synchronous pass; failures fall back to
the equivalent R path.

After refinement, labels are flood-filled into connected components. If
`K` is supplied, the shared adjacency-preserving exact-K engine merges
or splits those components and fails closed when the target is
topologically infeasible.
