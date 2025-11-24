# Commute Time Clustering

Performs spatially constrained clustering on a `NeuroVec` instance using
commute time distance (spectral embedding) and K-means clustering.

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
  connectivity = 27,
  weight_mode = c("binary", "heat"),
  noise_seed = NULL,
  verbose = TRUE
)
```

## Arguments

- bvec:

  A `NeuroVec` instance supplying the data to cluster.

- mask:

  A `NeuroVol` mask defining the voxels to include in the clustering
  result. If the mask contains `numeric` data, nonzero values will
  define the included voxels. If the mask is a
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.html),
  then `TRUE` will define the set of included voxels.

- K:

  The number of clusters to find. Default is 100.

- ncomp:

  The number of components to use for the commute time embedding.
  Default is the ceiling of `sqrt(K*2)`.

- alpha:

  A numeric value controlling the balance between spatial and feature
  similarity. Default is 0.5 (balanced). Range: 0 (spatial only) to 1
  (feature only).

- sigma1:

  A numeric value controlling the spatial weighting function. Default is
  0.73.

- sigma2:

  A numeric value controlling the feature weighting function. Default is
  5.

- connectivity:

  An integer representing the number of nearest neighbors to consider
  when constructing the similarity graph. Default is 27 (full 3D
  neighborhood).

- weight_mode:

  A character string indicating the type of weight function for the
  similarity graph. Options are "binary" and "heat". Default is "heat".

- noise_seed:

  Optional integer seed for reproducible noise injection when handling
  zero-variance voxels. If NULL (default), uses non-deterministic noise.
  **Note**: This only controls noise injection, not k-means
  initialization. For full reproducibility, wrap the entire call in
  [`set.seed()`](https://rdrr.io/r/base/Random.html).

- verbose:

  Logical. If TRUE, print progress messages. Default is TRUE.

## Value

A `list` of class `commute_time_cluster_result` (inheriting from
`cluster_result`) with the following elements:

- clusvol:

  An instance of type
  [ClusteredNeuroVol](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.html).

- cluster:

  A vector of cluster indices equal to the number of voxels in the mask.

- centers:

  A matrix of cluster centers (K x T) where T is the number of
  timepoints.

- coord_centers:

  A matrix of spatial coordinates (K x 3) with each row corresponding to
  a cluster centroid.

- embedding:

  The spectral embedding coordinates (N x ncomp) from commute time
  distance.

- n_clusters:

  The number of clusters (same as K).

- method:

  Character string "commute_time".

## Details

### Algorithm Overview

Commute time clustering uses spectral graph theory to embed voxels into
a lower-dimensional space where geodesic distances on the graph
approximate commute times (expected random walk return times). This
embedding respects both spatial proximity and feature similarity.

The algorithm has three main steps:

1.  **Graph Construction**: Build a weighted adjacency matrix combining
    spatial and feature similarity using
    [`neighborweights::weighted_spatial_adjacency()`](https://rdrr.io/pkg/neighborweights/man/weighted_spatial_adjacency.html).

2.  **Spectral Embedding**: Compute commute time distances via
    eigendecomposition of the graph Laplacian using
    [`neighborweights::commute_time_distance()`](https://rdrr.io/pkg/neighborweights/man/commute_time_distance.html).

3.  **Clustering**: Apply k-means to the embedded coordinates.

### Scalability Warning

**This method is computationally expensive and NOT recommended for
whole-brain clustering.**

- **Complexity**: O(N³) for eigendecomposition where N = number of
  voxels

- **Memory**: O(N²) for adjacency matrix storage

- **Practical limit**: ~10,000 voxels (ROI-based analysis)

- **Whole-brain**: ~100,000+ voxels will likely crash or take hours

For large-scale clustering, consider:

- [`slice_msf()`](slice_msf.md): Slice-based minimum spanning forests

- [`acsc()`](acsc.md): Adaptive correlation superclustering

- [`supervoxels()`](supervoxels.md) or [`snic()`](snic.md): Iterative
  spatial methods

### Parallelization Status

**Currently NOT explicitly parallelized.** The algorithm runs
sequentially, but matrix operations may use multi-threaded BLAS/LAPACK
libraries.

#### Why Not Parallelized:

- **External dependencies**: Uses `neighborweights` package functions

- **Eigendecomposition**: Difficult to parallelize efficiently in R

- **Already optimized**: BLAS/LAPACK typically use multiple threads
  automatically

- **Bottleneck**: Eigendecomposition dominates runtime regardless of
  parallelization

#### Performance Tips:

- **Use optimized BLAS**: OpenBLAS, Intel MKL, or Apple Accelerate

- **Reduce connectivity**: Smaller neighborhoods = sparser matrices
  (e.g., 6 or 18)

- **Increase alpha**: Higher values emphasize features over space,
  reducing graph density

- **Use fewer components**: Set `ncomp` lower (e.g., `ncomp = K/2`) for
  faster embedding

- **Pre-filter voxels**: Remove low-variance voxels before clustering

- **ROI analysis**: Apply to small regions of interest rather than whole
  brain

#### Common Issues:

**Eigenvalue Errors**: Often due to singular or near-singular weight
matrices.

Causes:

- Duplicate or perfectly correlated time series

- Disconnected graph components

- Insufficient connectivity parameter

Solutions:

- Increase `connectivity` (e.g., 27 instead of 6)

- Adjust `alpha` to balance spatial/feature weights

- Use `noise_seed` for reproducible noise injection

- Check for and remove constant voxels beforehand

**Memory Errors**: Adjacency matrix requires O(N²) memory.

Solutions:

- Reduce number of voxels (subsample or use smaller ROI)

- Use alternative method for large N

## See also

[`snic`](snic.md) for a faster non-iterative method
[`slice_msf`](slice_msf.md) for scalable slice-based clustering
[`acsc`](acsc.md) for large-scale correlation-based clustering

## Examples

``` r
if (FALSE) { # \dontrun{
# Small example with synthetic data
library(neuroim2)
mask <- NeuroVol(array(1, c(20, 20, 20)), NeuroSpace(c(20, 20, 20)))
vec <- replicate(10, NeuroVol(array(runif(20*20*20), c(20, 20, 20)),
  NeuroSpace(c(20, 20, 20))), simplify = FALSE)
vec <- do.call(concat, vec)

# Run clustering (8000 voxels - feasible for this method)
commute_res <- commute_cluster(vec, mask, K = 50, verbose = TRUE)

# Access results
print(commute_res$n_clusters)
plot(commute_res$clusvol)
} # }

if (FALSE) { # \dontrun{
# With reproducible noise injection (for zero-variance voxels)
commute_res <- commute_cluster(vec, mask, K = 50, noise_seed = 42)

# For full reproducibility (including k-means), use set.seed() wrapper
set.seed(123)
commute_res <- commute_cluster(vec, mask, K = 50, noise_seed = 42)

# ROI-based analysis (recommended workflow)
roi_mask <- mask  # In practice, use a smaller ROI
roi_mask[1:10, , ] <- 0  # Reduce voxels
commute_res <- commute_cluster(vec, roi_mask, K = 30)
} # }
```
