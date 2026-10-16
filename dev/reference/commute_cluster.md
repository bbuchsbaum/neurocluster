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
  eigen_solver = c("auto", "dense", "sparse"),
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

- eigen_solver:

  `"auto"` (default), `"dense"`, or `"sparse"`. The dense solver runs a
  full [`eigen()`](https://rdrr.io/r/base/eigen.html) on the N-by-N
  Laplacian: O(N^2) memory and O(N^3) time, practical up to a few
  thousand voxels. The sparse solver computes only the `ncomp + 1`
  smallest eigenpairs by shift-invert Lanczos
  ([`RSpectra::eigs_sym()`](https://rdrr.io/pkg/RSpectra/man/eigs.html),
  a suggested package) and scales to whole-brain masks. `"auto"` uses
  the dense solver for at most 3,000 masked voxels and the sparse solver
  above that. Both span the same eigenspace; embeddings agree to solver
  tolerance, up to the sign of each eigenvector (and the basis inside a
  repeated eigenvalue), which commute-time distances and hence k-means
  do not depend on.

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
(including absence of `.Random.seed`) is restored exactly. With the
dense solver it forms a full Laplacian eigendecomposition, requiring
quadratic memory and cubic time, so `eigen_solver = "auto"` switches to
the sparse solver above 3,000 voxels. With the sparse solver the run
time is dominated by k-means and is a few seconds for about 30,000
voxels and K = 200.
