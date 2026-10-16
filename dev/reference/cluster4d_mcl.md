# Cluster4d using sparse Markov Clustering (MCL)

Builds a sparse voxel graph, runs an MCL flow process, and maps
attractors to cluster labels. The implementation is optimized for sparse
neuroimaging graphs.

## Usage

``` r
cluster4d_mcl(
  vec,
  mask,
  n_clusters = 100,
  spatial_weight = 0.2,
  max_iterations = 8,
  connectivity = 6,
  verbose = FALSE,
  inflation = 1.6,
  expansion = 2L,
  loop_weight = 1,
  prune_k = NULL,
  prune_threshold = 1e-06,
  tol = 1e-04,
  feature_metric = c("correlation", "euclidean"),
  feature_sigma = NULL,
  spatial_sigma = NULL,
  exact_k = FALSE
)
```

## Arguments

- vec:

  A NeuroVec, SparseNeuroVec, or NeuroVol containing the features.

- mask:

  A NeuroVol mask; finite values strictly greater than zero are
  included.

- n_clusters:

  Requested cluster count. With `exact_k = FALSE`, this is the target
  used to select the closest natural MCL partition across a fixed,
  deterministic inflation search. With `exact_k = TRUE`, it is an exact
  adjacency-preserving postcondition.

- spatial_weight:

  Convex weight from zero through one for physical edge similarity.

- max_iterations:

  Positive MCL iteration limit for each inflation candidate.

- connectivity:

  Exact masked grid connectivity: 6, 18, or 26.

- verbose:

  Logical; report iteration diagnostics.

- inflation:

  Nominal MCL inflation parameter (\>1). In natural target mode, the
  deterministic search evaluates this value and fixed multipliers around
  it.

- expansion:

  MCL expansion power (integer \>= 2). Default 2.

- loop_weight:

  Finite non-negative self-loop weight added before normalization.

- prune_k:

  Maximum nonzero entries kept per column during MCL iterations. Must be
  a positive integer; the default is 64 and is independent of the target
  K.

- prune_threshold:

  Finite threshold in \[0, 1). The strongest entry in every column is
  retained even when all entries fall below the threshold.

- tol:

  Positive convergence tolerance on maximum absolute matrix delta.

- feature_metric:

  Feature similarity metric: `"correlation"` or `"euclidean"`.

- feature_sigma:

  Optional bandwidth for euclidean heat-kernel feature similarity.

- spatial_sigma:

  Optional bandwidth for spatial heat-kernel similarity.

- exact_k:

  If TRUE, use the shared adjacency-preserving exact-K engine.

## Value

A `cluster4d_result`. Metadata records natural and final K, convergence,
selected inflation, candidate search results, flow sparsity, final
normalization error, and exact-K policy.

## Details

Every expansion, inflation, and pruning step is checked to remain
finite, non-negative, and column stochastic. The `parallel` argument was
removed because this sparse Matrix implementation has no parallel
backend. Supplying it is therefore an ordinary unused-argument error
rather than a silently sequential execution.
