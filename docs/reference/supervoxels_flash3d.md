# FLASH-3D: Fast 3D Superclustering for fMRI

Performs spatially-constrained clustering using Fast Low-rank
Approximate Superclusters for Hemodynamics (FLASH-3D). This algorithm
uses DCT-based temporal hashing and 3D jump-flood propagation for
efficient clustering.

## Usage

``` r
supervoxels_flash3d(
  vec,
  mask,
  K,
  lambda_s = 0.6,
  lambda_t = 1,
  lambda_g = 0,
  rounds = 2L,
  bits = 64L,
  dctM = 12L,
  vox_scale = NULL,
  barrier = NULL,
  verbose = FALSE
)
```

## Arguments

- vec:

  A `NeuroVec` instance supplying the 4D data to cluster

- mask:

  A `NeuroVol` mask defining the voxels to include in the clustering
  result. If the mask contains `numeric` data, nonzero values will
  define the included voxels. If the mask is a
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.html),
  then `TRUE` will define the set of included voxels.

- K:

  The number of clusters to find

- lambda_s:

  Spatial weight for distance penalty (default 0.6). Will be annealed
  upward over rounds.

- lambda_t:

  Temporal weight for Hamming distance of time-series hashes (default
  1.0)

- lambda_g:

  Optional barrier weight (default 0.0). If \>0, supply `barrier`
  volume.

- rounds:

  Number of outer rounds (seed→flood→recenter). Default 2, typically 2-3
  is sufficient.

- bits:

  Length of temporal hash (64 or 128). Default 64.

- dctM:

  Number of low DCT coefficients to rank-hash (default 12, range 4-32)

- vox_scale:

  Voxel size scaling for spatial distance, e.g., c(dx,dy,dz). Default
  c(1,1,1)

- barrier:

  Optional 3D numeric array same dimensions as mask; higher values
  resist region growth

- verbose:

  Logical indicating whether to print progress messages

## Value

A `list` of class `cluster_result` with the following elements:

- clusvol:

  An instance of type
  [ClusteredNeuroVol](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.html)

- cluster:

  An integer vector of cluster assignments for each voxel in mask

- centers:

  A matrix of cluster centers in feature space (K x T)

- coord_centers:

  A matrix of spatial coordinates of cluster centers (K x 3)

- K:

  The number of clusters

- method:

  A character string indicating the method used ("FLASH-3D")

## Details

FLASH-3D uses a novel approach combining:

- DCT-based temporal feature hashing for fast similarity computation

- Jump-flood algorithm for efficient spatial propagation

- Blue-noise seeding for optimal initial cluster placement

- Annealing of spatial weights to encourage compact clusters

The algorithm is particularly efficient for large-scale fMRI data,
offering significant speed improvements over iterative methods while
maintaining clustering quality.

## Note

Consider using [`cluster4d`](cluster4d.md) with `method = "flash3d"` for
a standardized interface across all clustering methods.

## References

FLASH-3D algorithm for fast superclustering of fMRI data (2025)

## See also

[`supervoxels`](supervoxels.md), [`snic`](snic.md),
[`slic4d_supervoxels`](slic4d_supervoxels.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  # Basic usage
  result <- supervoxels_flash3d(vec, mask, K = 100)
  
  # With custom parameters
  result <- supervoxels_flash3d(vec, mask, K = 100, 
                                lambda_s = 0.8, lambda_t = 1.2,
                                bits = 128, dctM = 16)
  
  # With barrier for anatomy-aware clustering
  barrier_vol <- create_anatomical_barrier(mask)
  result <- supervoxels_flash3d(vec, mask, K = 100,
                                lambda_g = 0.5, barrier = barrier_vol)
} # }
```
