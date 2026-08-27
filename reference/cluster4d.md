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
  method = c("supervoxels", "snic", "slic", "corr_slic", "brs_slic", "slice_msf",
    "flash3d", "g3s", "rena", "rena_plus", "mcl", "acsc", "commute"),
  spatial_weight = 0.5,
  max_iterations = NULL,
  connectivity = NULL,
  parallel = TRUE,
  verbose = FALSE,
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

- method:

  Clustering algorithm to use. Options:

  - `"supervoxels"`: Iterative heat kernel-based clustering (default)

  - `"snic"`: Simple Non-Iterative Clustering

  - `"slic"`: SLIC superpixels extended to 4D

  - `"brs_slic"`: Boundary-Refined Sketch SLIC (coarse sketch + boundary
    exact-correlation refinement)

  - `"slice_msf"`: Slice-wise Minimum Spanning Forest (fast but may show
    z-artifacts)

  - `"flash3d"`: Fast Low-rank Approximate Superclusters

  - `"g3s"`: Gradient-Guided Geodesic Supervoxels (NEW - recommended for
    best quality/speed)

  - `"rena"`: Recursive Nearest Agglomeration (fast, balanced,
    topology-aware)

  - `"mcl"`: Sparse Markov Cluster Algorithm on a weighted voxel graph

  - `"acsc"`: Adaptive Correlation Superclustering (graph-based with
    boundary refinement)

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

- connectivity:

  Integer neighborhood connectivity. Supported values are
  method-specific: supervoxels accepts 6, 18, 26, or 27; slic,
  corr_slic, brs_slic, and slice_msf accept 6 or 26; G3S, ReNA variants,
  and MCL accept 6, 18, or 26. It is inactive for snic, flash3d, acsc,
  and commute.

- parallel:

  Enable parallel processing for supervoxels and slic. For all other
  methods, an explicit `TRUE` is rejected; omitted or `FALSE` selects
  the serial implementation. Default TRUE for supported methods.

- verbose:

  Print progress messages. Default FALSE.

- ...:

  Additional method-specific parameters. See method documentation for
  details.

## Value

A `cluster4d_result` object (also inherits from `cluster_result`)
containing:

- labels:

  Contiguous positive integer final labels in mask voxel order

- cluster:

  Backward-compatible alias that is identical to `labels`

- clusvol:

  A `ClusteredNeuroVol` with cluster assignments

- centers:

  Actual final-label means in the original feature space, with shape
  actual K by timepoints

- coord_centers:

  Actual final-label coordinate means in physical mm, with shape actual
  K by 3

- actual_k:

  Actual number of clusters produced

- n_clusters:

  Backward-compatible alias of `actual_k`

- label_ids:

  Explicit mapping from center rows to contiguous label IDs

- method:

  Clustering method used

- parameters:

  List of all parameters used

- provenance:

  Typed label, original-feature, physical-coordinate, geometry, and mask
  provenance

- metadata:

  Method-specific additional information

## Algorithm Comparison

|  |  |  |  |  |  |  |
|----|----|----|----|----|----|----|
| **Method** | **Speed** | **3D Continuity** | **Memory** | **Parallel** | **Best For** | supervoxels |
| Slow | Excellent | High | Yes | Small-medium data, smooth parcels | snic | Fast |
| Good | Low | No | Large data, non-iterative | slic | Fast | Good |
| Medium | Yes | Balanced speed/quality | slice_msf | Very Fast | Moderate | Low |
| No | High-res data, accept z-artifacts | flash3d | Fast | Good | Medium | No |
| Large data, hash-based | rena | Fast | Excellent | Low | No | Balanced clusters, topology-aware |
| mcl | Fast | Good | Medium | No | Sparse graph clustering with tunable granularity |  |

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
[`cluster4d_supervoxels`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_supervoxels.md),
[`cluster4d_snic`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_snic.md),
[`cluster4d_slic`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_slic.md),
[`cluster4d_slice_msf`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_slice_msf.md),
[`cluster4d_flash3d`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_flash3d.md),
[`cluster4d_mcl`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_mcl.md),
[`cluster4d_commute`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_commute.md)

Legacy functions (deprecated):
[`supervoxels`](https://bbuchsbaum.github.io/neurocluster/reference/supervoxels.md),
[`snic`](https://bbuchsbaum.github.io/neurocluster/reference/snic.md),
[`slic4d_supervoxels`](https://bbuchsbaum.github.io/neurocluster/reference/slic4d_supervoxels.md),
[`slice_msf`](https://bbuchsbaum.github.io/neurocluster/reference/slice_msf.md),
[`supervoxels_flash3d`](https://bbuchsbaum.github.io/neurocluster/reference/supervoxels_flash3d.md)

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
