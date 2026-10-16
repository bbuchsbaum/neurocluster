# Commute-time clustering on a physical k-nearest-neighbor graph

Builds an explicitly declared symmetric k-nearest-neighbor graph in
physical coordinate space, computes a Laplacian commute-time embedding,
and applies deterministically initialized k-means.

## Usage

``` r
commute_cluster(
  bvec,
  mask,
  K = 100,
  ncomp = ceiling(sqrt(K * 2)),
  alpha = 0.5,
  sigma1 = 0.73,
  sigma2 = 5,
  connectivity = 6L,
  weight_mode = c("binary", "heat"),
  bad_voxel_policy = c("error", "zero_constant"),
  variance_tol = 1e-09,
  eigen_tol = 1e-10,
  verbose = FALSE
)
```

## Arguments

- bvec:

  A NeuroVec or SparseNeuroVec supplying voxel features.

- mask:

  A NeuroVol mask. Finite values strictly greater than zero are
  included.

- K:

  Finite integer cluster count from one through the included voxel
  count.

- ncomp:

  Finite integer embedding dimension from one through N - 1.

- alpha:

  Feature-similarity weight from zero through one.

- sigma1:

  Positive physical-distance heat-kernel bandwidth.

- sigma2:

  Positive feature-distance heat-kernel bandwidth.

- connectivity:

  Number of nearest physical neighbors per directed search. The final
  undirected graph is the symmetric union of those directed edges; this
  is not masked grid adjacency and may bridge mask holes.

- weight_mode:

  `"binary"` or `"heat"` edge weights.

- bad_voxel_policy:

  `"error"` (default), or `"zero_constant"` to map finite zero-variance
  voxel series to deterministic zero standardized vectors. NA and Inf
  always fail closed. No noise is injected.

- variance_tol:

  Non-negative threshold defining a constant voxel series.

- eigen_tol:

  Positive tolerance for the Laplacian zero eigenspace.

- verbose:

  Logical; print progress.

## Value

A `commute_time_cluster_result` with contiguous labels, K-by-T
original-space feature means, K-by-3 physical coordinate means, the
N-by-ncomp embedding, graph and spectral provenance, and explicit
data-policy provenance.

## Details

The implementation is deterministic. Any RNG initialization performed
inside [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) is
scoped to a fixed local seed, after which the caller's prior RNG state
(including absence of `.Random.seed`) is restored exactly. It forms a
dense Laplacian eigendecomposition, requiring quadratic memory and cubic
time; use only for small regions.
