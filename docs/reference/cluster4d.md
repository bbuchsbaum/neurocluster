# Unified 4D Clustering for Neuroimaging Data

Performs spatially-constrained clustering on 4D neuroimaging data using
various algorithms. This is the main entry point for all clustering
methods in the neurocluster package.

## Usage

``` r
cluster4d(
  vec,
  mask,
  n_clusters = 100,
  method = c("supervoxels", "snic", "slic", "slice_msf", "flash3d", "g3s", "rena",
    "rena_plus"),
  spatial_weight = 0.5,
  max_iterations = 10,
  connectivity = 26,
  parallel = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- vec:

  A `NeuroVec` instance supplying the 4D data (x, y, z, time) to cluster

- mask:

  A `NeuroVol` mask defining the voxels to include in clustering. If
  numeric, nonzero values define included voxels. If logical, TRUE
  values define included voxels.

- n_clusters:

  Target number of clusters (default 100). Note that some methods may
  produce slightly different numbers of clusters due to algorithmic
  constraints.

- method:

  Clustering algorithm to use. Options:

  - `"supervoxels"`: Iterative heat kernel-based clustering (default)

  - `"snic"`: Simple Non-Iterative Clustering

  - `"slic"`: SLIC superpixels extended to 4D

  - `"slice_msf"`: Slice-wise Minimum Spanning Forest (fast but may show
    z-artifacts)

  - `"flash3d"`: Fast Low-rank Approximate Superclusters

  - `"g3s"`: Gradient-Guided Geodesic Supervoxels (NEW - recommended for
    best quality/speed)

  - `"rena"`: Recursive Nearest Agglomeration (fast, balanced,
    topology-aware)

- spatial_weight:

  Balance between spatial and feature similarity (0-1). Higher values
  emphasize spatial compactness. Default 0.5. Maps to method-specific
  parameters:

  - supervoxels: `alpha = 1 - spatial_weight` (0 = all spatial, 1 = all
    feature)

  - snic/slic: `compactness = spatial_weight * 20` (typical range 1-20)

  - slice_msf: `compactness = spatial_weight * 10` (typical range 1-10)

  - flash3d: `lambda_s = spatial_weight` (direct mapping)

- max_iterations:

  Maximum iterations for iterative methods. Default 10. Maps to:
  `iterations` (supervoxels), `max_iter` (snic/slic), `rounds`
  (flash3d).

- connectivity:

  Neighborhood connectivity (6 or 26). Default 26. 6 = face neighbors
  only, 26 = face + edge + corner neighbors.

- parallel:

  Enable parallel processing where supported. Default TRUE.

- verbose:

  Print progress messages. Default FALSE.

- ...:

  Additional method-specific parameters. See method documentation for
  details.

## Value

A `cluster4d_result` object (also inherits from `cluster_result`)
containing:

- clusvol:

  A `ClusteredNeuroVol` with cluster assignments

- cluster:

  Integer vector of cluster assignments for each masked voxel

- centers:

  Matrix of cluster centers in feature space (n_clusters x timepoints)

- coord_centers:

  Matrix of cluster spatial centers (n_clusters x 3)

- n_clusters:

  Actual number of clusters produced

- method:

  Clustering method used

- parameters:

  List of all parameters used

- metadata:

  Method-specific additional information

## Algorithm Comparison

|  |  |  |  |  |  |  |
|----|----|----|----|----|----|----|
| **Method** | **Speed** | **3D Continuity** | **Memory** | **Parallel** | **Best For** | supervoxels |
| Slow | Excellent | High | Yes | Small-medium data, smooth parcels | snic | Fast |
| Good | Low | No | Large data, non-iterative | slic | Fast | Good |
| Medium | Yes | Balanced speed/quality | slice_msf | Very Fast | Moderate | Low |
| Yes | High-res data, accept z-artifacts | flash3d | Fast | Good | Medium | Partial |
| Large data, hash-based | rena | Fast | Excellent | Low | No | Balanced clusters, topology-aware |

## Parameter Guidelines

**For whole-brain parcellation:**

- n_clusters: 100-1000 depending on desired granularity

- spatial_weight: 0.4-0.6 for balanced clustering

- connectivity: 26 for smoother boundaries

**For ROI analysis:**

- n_clusters: 10-100 depending on ROI size

- spatial_weight: 0.2-0.4 to emphasize functional similarity

- connectivity: 6 for more discrete parcels

**For high-resolution data (\< 2mm):**

- method: "slice_msf" or "flash3d" for speed

- n_clusters: Scale with voxel count (roughly n_voxels/200)

## See also

Method-specific functions:
[`cluster4d_supervoxels`](cluster4d_supervoxels.md),
[`cluster4d_snic`](cluster4d_snic.md),
[`cluster4d_slic`](cluster4d_slic.md),
[`cluster4d_slice_msf`](cluster4d_slice_msf.md),
[`cluster4d_flash3d`](cluster4d_flash3d.md)

Legacy functions (deprecated): [`supervoxels`](supervoxels.md),
[`snic`](snic.md), [`slic4d_supervoxels`](slic4d_supervoxels.md),
[`slice_msf`](slice_msf.md),
[`supervoxels_flash3d`](supervoxels_flash3d.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simple synthetic example (runs quickly for testing)
library(neuroim2)
mask <- NeuroVol(array(1, c(4,4,4)), NeuroSpace(c(4,4,4)))
vec <- NeuroVec(array(rnorm(4*4*4*10), c(4,4,4,10)), 
                NeuroSpace(c(4,4,4,10)))
result <- cluster4d(vec, mask, n_clusters = 3, method = "g3s", 
                   max_iterations = 1)
print(result$n_clusters)
} # }

if (FALSE) { # \dontrun{
# More realistic examples with larger data
mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
vec <- replicate(50, NeuroVol(array(runif(20*20*20), c(20,20,20)),
                              NeuroSpace(c(20,20,20))), simplify=FALSE)
vec <- do.call(concat, vec)

# Basic usage with default supervoxels method
result <- cluster4d(vec, mask, n_clusters = 100)

# Fast clustering with FLASH-3D (hash-based)
result <- cluster4d(vec, mask, n_clusters = 100, method = "flash3d")

# Emphasize spatial compactness
result <- cluster4d(vec, mask, n_clusters = 100, spatial_weight = 0.8)

# Use specific method with custom parameters
result <- cluster4d(vec, mask, n_clusters = 100, 
                   method = "slice_msf",
                   num_runs = 3,  # slice_msf-specific parameter
                   consensus = TRUE)

# Get parameter suggestions for your data
n_vox <- sum(mask > 0)
n_time <- dim(vec)[4]
params <- suggest_cluster4d_params(n_vox, n_time, priority = "quality")
result <- cluster4d(vec, mask, 
                   n_clusters = params$n_clusters,
                   method = params$recommended_method)
} # }
```
