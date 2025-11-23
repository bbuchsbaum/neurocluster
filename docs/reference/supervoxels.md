# Supervoxel Clustering (3D volumes)

Cluster a `NeuroVec` instance into a set of spatially constrained
clusters.

## Usage

``` r
supervoxels(
  bvec,
  mask,
  K = 500,
  sigma1 = 1,
  sigma2 = 2.5,
  iterations = 50,
  connectivity = 27,
  use_medoid = FALSE,
  use_gradient = TRUE,
  alpha = 0.5,
  parallel = TRUE,
  grain_size = 100,
  verbose = FALSE,
  converge_thresh = 0.001
)
```

## Arguments

- bvec:

  A
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html)
  instance supplying the data to cluster.

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

  Matrix of cluster centers in feature space.

- coord_centers:

  Matrix of cluster spatial centroids.

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

### Parallelization Status

**NOW PARALLELIZED with RcppParallel!** The supervoxels algorithm can
now run cluster assignment updates in parallel across multiple CPU
cores.

#### Parallel Operations:

1.  **Heat Kernel Computation**: Parallel across voxels using
    RcppParallel

    - Each voxel's best cluster assignment computed independently

    - Automatic load balancing with configurable grain size

    - Scales linearly with number of CPU cores

2.  **Sequential Operations** (still):

    - Initialization (K-means or gradient-based seed selection)

    - Centroid updates after each iteration

    - Convergence checking between iterations

#### Performance Characteristics:

- **Speedup**: Typically 2-8x faster on multicore systems

- **Automatic optimization**: Disabled for small datasets (\<1000
  voxels)

- **Memory overhead**: Minimal - uses shared memory via RcppParallel

- **Computational complexity**: Still O(N × K × iterations) but
  parallelized over N

#### Parallel Configuration:

- **parallel**: Set to FALSE to force sequential execution

- **grain_size**: Controls work distribution (default 100)

  - Smaller values = better load balancing but more overhead

  - Larger values = less overhead but potential imbalance

- **Thread control**: Set threads via
  [`RcppParallel::setThreadOptions()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html)

#### When Parallelization Helps Most:

- Large numbers of voxels (N \> 10,000)

- Many clusters (K \> 100)

- Multiple iterations needed for convergence

- Systems with 4+ CPU cores

#### Performance Tips:

- **Set threads**: `RcppParallel::setThreadOptions(numThreads = 4)`

- **Tune grain_size**: Start with nvoxels/nthreads/10

- **Monitor CPU usage**: Should see near 100% on all cores during
  updates

- **Memory considerations**: Parallel version uses slightly more RAM

- **Disable for debugging**: Set `parallel = FALSE` for reproducible
  debugging

## Note

Consider using [`cluster4d`](cluster4d.md) with `method = "supervoxels"`
for a standardized interface across all clustering methods.

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
