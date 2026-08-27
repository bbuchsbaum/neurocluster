# Supervoxel Clustering (3D volumes)

Cluster a `NeuroVec` instance into a set of spatially constrained
clusters.

## Usage

``` r
supervoxels(
  bvec,
  mask,
  K = 500,
  sigma1 = NULL,
  sigma2 = 2.5,
  iterations = 50,
  connectivity = 27,
  use_medoid = FALSE,
  use_gradient = TRUE,
  alpha = 0.5,
  parallel = TRUE,
  grain_size = 100,
  num_threads = NULL,
  verbose = FALSE,
  converge_thresh = 0.001
)
```

## Arguments

- bvec:

  A
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html)
  instance supplying the data to cluster. Can also be a 3D
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.html)
  for structural image segmentation, which will be automatically
  converted to a single-timepoint NeuroVec internally.

- mask:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.html)
  mask defining the voxels to include. If numeric, nonzero = included.

- K:

  The number of clusters to find (default 500).

- sigma1:

  The bandwidth of the heat kernel for the data vectors.

- sigma2:

  The bandwidth of the heat kernel for the coordinate vectors.

- iterations:

  The maximum number of cluster iterations.

- connectivity:

  The number of nearest neighbors defining the neighborhood.

- use_medoid:

  Logical; whether to use medoids rather than means for cluster centers.

- use_gradient:

  Logical; use the image gradient to initialize clusters if possible.

- alpha:

  The relative weighting of data similarity vs spatial similarity;
  `alpha=1` = all data weighting, `alpha=0` = purely spatial weighting.

- parallel:

  Logical; whether to use parallel processing for cluster assignment
  updates. Default is TRUE. Parallel processing is automatically
  disabled for small datasets (\<1000 voxels).

- grain_size:

  Integer; the minimum number of voxels to process per parallel task.
  Default is 100. Smaller values provide better load balancing but
  increase overhead.

- num_threads:

  Optional integer to override the number of threads used by
  RcppParallel (defaults to package/global setting). Ignored if
  `parallel = FALSE`.

- verbose:

  Logical; whether to print detailed progress messages including
  convergence metrics. Default is FALSE.

- converge_thresh:

  Numeric; convergence threshold as proportion of voxels switching
  clusters. Algorithm stops when switch ratio falls below this value.
  Default is 0.001 (0.1% of voxels).

## Value

A `list` (of class `cluster_result`) with elements:

- clusvol:

  `ClusteredNeuroVol` containing the final clustering.

- cluster:

  Integer vector of cluster assignments for each voxel.

- centers:

  K-by-T matrix of final cluster means in the original feature space.

- coord_centers:

  K-by-3 matrix of final spatial means in physical units.

- metadata:

  Algorithm diagnostics, including whether parallel assignment was
  requested and used and any deterministic empty-label reseeds.

## Details

The algorithm:

1.  Scale input data (`bvec`) so each feature dimension is centered and
    scaled.

2.  If `use_gradient = TRUE`, initialize cluster seeds using
    gradient-based heuristics.

3.  Run an iterative, spatially-constrained clustering that updates
    voxel assignments based on both feature similarity (bandwidth
    `sigma1`) and spatial proximity (bandwidth `sigma2`), weighted by
    `alpha`.

4.  Return the final clusters, plus the feature-space and
    coordinate-space centers.

### Exact-K and parallel contracts

Every assignment center must have a positive membership count. If an
update empties one or more labels, distinct voxels are deterministically
moved from non-singleton clusters in descending combined feature/spatial
loss order. Empty all-zero centers are therefore never passed to an
assignment kernel. Since `K <= N` is required, this reseed policy is
feasible and the returned labels, center rows, and requested K agree on
every convergence path.

With `parallel = TRUE`, voxel assignment uses RcppParallel when there
are more than 1000 included voxels; smaller inputs use the same
sequential scoring kernel. Parallel and sequential kernels implement the
same scoring and deterministic tie rule. Centroid reduction may also run
in parallel for sufficiently large K. `metadata$algorithm$parallel_used`
reports whether the parallel assignment path actually ran. No
architecture-specific fallback is applied.

## Note

Consider using
[`cluster4d`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md)
with `method = "supervoxels"` for a standardized interface across all
clustering methods.

## Examples

``` r
if (FALSE) { # \dontrun{
mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
bvec <- replicate(10,
                  NeuroVol(array(runif(20*20*20), c(20,20,20)),
                           NeuroSpace(c(20,20,20))),
                  simplify=FALSE)
bvec <- do.call(concat, bvec)
cres1 <- supervoxels(bvec, mask, K=100, sigma1=1, sigma2=3)
} # }
```
