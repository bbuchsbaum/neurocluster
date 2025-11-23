# SLiCE-MSF: Slice-wise, Low-rank, Minimum-Spanning Forest Clustering

Performs spatially constrained clustering on neuroimaging time series
data using a slice-based approach with optional 3D stitching. The
algorithm applies DCT sketching for temporal compression, reliability
weighting, and graph-based segmentation. While computationally
efficient, the slice-based approach may create visible boundaries
between z-slices (see Details for mitigation strategies).

## Usage

``` r
slice_msf(
  vec,
  mask,
  target_k_global = -1,
  target_k_per_slice = -1,
  r = 12,
  compactness = 5,
  min_size = 80,
  num_runs = 3,
  consensus = TRUE,
  stitch_z = TRUE,
  theta_link = 0.85,
  min_contact = 1,
  nbhd = 8,
  gamma = 1.5,
  k_fuse = NULL,
  min_size_fuse = NULL,
  use_features = FALSE,
  lambda = 0.7,
  z_mult = 0
)
```

## Arguments

- vec:

  A `NeuroVec` or `SparseNeuroVec` instance supplying the time series
  data to cluster.

- mask:

  A `NeuroVol` mask defining the voxels to include in the clustering
  result. If the mask contains `numeric` data, nonzero values will
  define the included voxels. If the mask is a
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.html),
  then `TRUE` will define the set of included voxels.

- target_k_global:

  Target number of clusters across entire volume. When positive, uses
  region adjacency graph (RAG) agglomeration to achieve exactly K
  clusters. Default is -1 (no target, uses natural
  Felzenszwalb-Huttenlocher clustering).

- target_k_per_slice:

  Target number of clusters per slice. Only used when positive and
  stitch_z=FALSE. Useful for consistent parcellation across slices.
  Default is -1.

- r:

  DCT sketch rank (number of basis functions, excluding DC component).
  Higher values preserve more temporal detail but increase computation.
  Range: 4-20, Default: 12.

- compactness:

  Controls spatial compactness vs feature similarity balance (1-10).
  Lower values (1-3): Feature-driven, may create irregular shapes but
  better cross-slice alignment. Medium (4-6): Balanced. Higher (7-10):
  Spatially compact but may show more z-artifacts. Default is 5.

- min_size:

  Minimum cluster size in voxels. Smaller clusters are merged with
  nearest neighbors. Affects granularity of parcellation. Default is 80.

- num_runs:

  Number of independent segmentation runs. Single run (1) is faster but
  less stable. Multiple runs (3-5) with consensus fusion reduce
  variability and can smooth z-transitions. More than 5 has diminishing
  returns. Default is 3.

- consensus:

  Logical. If TRUE and num_runs \> 1, applies consensus fusion across
  runs, improving stability and potentially reducing z-artifacts.
  Default is TRUE.

- stitch_z:

  Logical. If TRUE, attempts to stitch 2D slice clusters into coherent
  3D clusters by merging across z-boundaries. Essential for 3D
  continuity. Default is TRUE.

- theta_link:

  Centroid correlation threshold for cross-slice stitching (0-1). Lower
  values (0.70-0.80): Aggressive stitching, reduces z-plane artifacts
  but may over-merge. Default (0.85): Balanced. Higher (0.90-0.95):
  Conservative, preserves boundaries but more z-artifacts. Critical
  parameter for z-continuity.

- min_contact:

  Minimum number of touching voxels between slices required for
  stitching attempt. Lower (1-2): More connections, better continuity.
  Higher (3-5): Stricter requirement, prevents spurious bridges. Default
  is 1.

- nbhd:

  Neighborhood connectivity for within-slice clustering. Options: 4 (von
  Neumann), 6 (includes z but mapped to 8), 8 (Moore). Higher
  connectivity can improve within-slice coherence. Default is 8.

- gamma:

  Reliability weighting exponent for split-half correlation. Higher
  values (\>1.5) emphasize high-reliability voxels, useful for noisy
  data. Lower values (\<1) treat all voxels more equally. Default is
  1.5.

- k_fuse:

  Scale parameter for consensus fusion graph. If NULL, uses same as
  compactness-derived scale. Lower values create more clusters in
  fusion. Default is NULL.

- min_size_fuse:

  Minimum cluster size during consensus fusion. If NULL, uses min_size.
  Can be set lower to preserve small consistent regions. Default is
  NULL.

- use_features:

  Include feature similarity in consensus fusion. When TRUE, uses both
  label agreement and temporal similarity, improving cross-slice
  consistency. Recommended when targeting exact K. Default is FALSE.

- lambda:

  Mixing parameter for consensus (0-1). Controls balance between label
  agreement and feature similarity when use_features=TRUE. Higher values
  weight label agreement more. Default is 0.7.

- z_mult:

  Z-smoothing factor (0-1). Values \> 0 softly blend DCT sketches
  between adjacent slices before clustering, reducing visible z-plane
  seams. 0 preserves the legacy per-slice behavior. Recommended range
  0.1-0.4. Default is 0.0.

## Value

A `list` of class `slice_msf_cluster_result` with the following
elements:

- clusvol:

  An instance of type
  [ClusteredNeuroVol](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.html).

- cluster:

  A vector of cluster indices equal to the number of voxels in the mask.

- centers:

  A matrix of cluster centers with each column representing the feature
  vector for a cluster.

- coord_centers:

  A matrix of spatial coordinates with each row corresponding to a
  cluster.

- runs:

  If num_runs \> 1, a list of individual run results.

## Details

### Algorithm Overview

SLiCE-MSF (Slice-wise, Low-rank, Minimum-Spanning Forest) is designed
for fast clustering of high-resolution fMRI data. The algorithm proceeds
in stages:

1.  **Temporal Sketching**: Each voxel's time series is compressed using
    Discrete Cosine Transform (DCT) basis functions, reducing from T
    timepoints to r coefficients.

2.  **Reliability Weighting**: Split-half correlations are computed to
    identify reliable voxels, which receive higher weight in clustering
    decisions.

3.  **Slice-wise Clustering**: Each axial slice is clustered
    independently using Felzenszwalb-Huttenlocher graph segmentation,
    which efficiently finds regions with high internal similarity and
    low external similarity.

4.  **Z-Stitching** (optional): Clusters from adjacent slices are merged
    based on spatial contact and centroid similarity to create 3D
    parcels.

5.  **Consensus Fusion** (optional): Multiple runs are combined using
    co-association matrices to improve stability.

### Z-Plane Artifacts and Mitigation

Because clustering is performed slice-by-slice, the algorithm can
produce visible horizontal lines when viewed in sagittal or coronal
planes. These artifacts are inherent to the slice-based approach but can
be minimized:

#### Artifact Reduction Strategies

**Mild artifacts** (slight discontinuities):

- Reduce theta_link to 0.75-0.80 for more aggressive stitching

- Set min_contact to 2-3 for balanced connectivity

- Use lower compactness (2-4) for more flexible cluster shapes

- Set z_mult between 0.1-0.3 to softly bleed information across slices

**Moderate artifacts** (visible lines):

- Use multiple runs (num_runs = 3-5) with consensus = TRUE

- Enable use_features = TRUE for feature-based consensus

- Adjust lambda to 0.5-0.6 to weight features more

- Combine with z_mult smoothing for additional continuity without heavy
  fusion

**Severe artifacts** (strong discontinuities):

- Consider alternative algorithms like supervoxels() or snic() for true
  3D clustering

- These provide smoother 3D parcels at the cost of increased computation
  time

### Parameter Selection Guidelines

#### For Whole-Brain Parcellation

- r = 10-15 (balance detail and speed)

- compactness = 4-6 (balanced spatial/feature weighting)

- min_size = 80-150 (appropriate for ~3mm resolution)

- num_runs = 3-5 with consensus = TRUE

#### For ROI Analysis

- r = 15-20 (preserve more temporal detail)

- compactness = 2-4 (feature-driven clustering)

- min_size = 20-50 (allow smaller parcels)

- theta_link = 0.75 (aggressive stitching within ROI)

#### For High-Resolution Data (\< 2mm)

- Increase min_size proportionally (e.g., 200-300 voxels)

- Use target_k_global to control final parcel count

- Consider target_k_per_slice for consistent slice-wise parcellation

### Performance Considerations

- **Memory**: Scales with mask size × r × num_runs

- **Speed**: Much faster than full 3D methods, especially for
  high-resolution data

- **Trade-off**: Speed vs. z-continuity - for smooth 3D parcels, use
  supervoxels()

### Parallelization Strategy

SLiCE-MSF uses **automatic parallelization via RcppParallel** for key
operations:

#### Parallel Operations:

1.  **DCT Sketching** (`SliceSketchWorker`): Each z-slice is processed
    in parallel

    - Detrending and z-scoring of time series

    - Split-half reliability computation

    - DCT coefficient calculation

    - Threads automatically assigned by RcppParallel based on available
      cores

2.  **Consensus Fusion** (`FuseSliceWorker`): When num_runs \> 1

    - Co-association matrix computation parallelized across slices

    - Edge weight calculations done in parallel

    - Graph segmentation remains sequential within each slice

#### Performance Characteristics:

- **Scaling**: Near-linear speedup with cores for the sketching phase

- **Optimal for**: High-resolution data (many slices to process)

- **Automatic**: No user configuration needed - uses all available cores

- **Memory-efficient**: Each thread processes independent slices

- **Bottleneck**: Graph segmentation phase is sequential per slice

#### Controlling Parallelization:

RcppParallel automatically determines the number of threads. To control:

    # Set number of threads globally
    RcppParallel::setThreadOptions(numThreads = 4)

    # Reset to automatic
    RcppParallel::setThreadOptions(numThreads = "auto")

## Troubleshooting Z-Plane Artifacts

Common issues and solutions:

**Horizontal lines in sagittal/coronal views:**

1.  First, verify artifacts aren't anatomically meaningful (e.g., at
    gray/white boundaries)

2.  Try progressive adjustments:

    - Step 1: theta_link = 0.75 (from default 0.85)

    - Step 2: Add num_runs = 3, consensus = TRUE

    - Step 3: Set use_features = TRUE, lambda = 0.5

    - Step 4: Reduce compactness to 2-3

**Over-merging after reducing theta_link:**

- Increase min_contact to 3-5 to require more evidence for merging

- Increase min_size to prevent small bridging clusters

**Inconsistent cluster sizes across slices:**

- Use target_k_per_slice (with stitch_z = FALSE) for consistent slice
  parcellation

- Or use target_k_global for consistent total cluster count

**Poor performance in low-SNR regions:**

- Increase gamma to 2.0-2.5 to emphasize reliable voxels

- Consider masking out low-SNR regions before clustering

## References

Felzenszwalb, P. F., & Huttenlocher, D. P. (2004). Efficient graph-based
image segmentation. International journal of computer vision, 59(2),
167-181.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example data
mask <- NeuroVol(array(1, c(64,64,32)), NeuroSpace(c(64,64,32)))
vec <- replicate(100, NeuroVol(array(rnorm(64*64*32), c(64,64,32)),
                 NeuroSpace(c(64,64,32))), simplify=FALSE)
vec <- do.call(concat, vec)

# Example 1: Basic usage (may show z-plane artifacts)
result_basic <- slice_msf(vec, mask, 
                          num_runs = 1,
                          compactness = 5)

# Example 2: Reduced z-artifacts with aggressive stitching
# Recommended for minimizing slice boundaries
result_smooth <- slice_msf(vec, mask, 
                           theta_link = 0.75,      # More aggressive stitching
                           min_contact = 2,        # Moderate contact requirement
                           compactness = 3,        # Lower compactness
                           num_runs = 3,           # Multiple runs
                           consensus = TRUE,       # Enable consensus
                           use_features = TRUE,    # Feature-based consensus
                           lambda = 0.6)           # Balance features/labels

# Example 3: Conservative approach for anatomical boundaries
# Preserves natural boundaries but may show more z-artifacts
result_conservative <- slice_msf(vec, mask,
                                 theta_link = 0.90,    # Conservative stitching
                                 min_contact = 5,      # Strict contact requirement
                                 compactness = 7,      # Compact clusters
                                 min_size = 120)       # Larger minimum size

# Example 4: Exact K targeting with 100 clusters
result_exact <- slice_msf(vec, mask, 
                         target_k_global = 100,   # Exactly 100 clusters
                         use_features = TRUE,     # Required for exact K
                         num_runs = 3,
                         consensus = TRUE)

# Example 5: Per-slice consistency (useful for group studies)
result_per_slice <- slice_msf(vec, mask,
                              target_k_per_slice = 50,  # 50 clusters per slice
                              stitch_z = FALSE,         # No z-stitching
                              compactness = 5)

# Example 6: High-resolution data optimization
# For data with voxel size < 2mm
result_highres <- slice_msf(vec, mask,
                            min_size = 250,         # Larger clusters for high-res
                            r = 15,                 # More DCT components
                            compactness = 4,
                            theta_link = 0.78,
                            num_runs = 5,
                            consensus = TRUE)

# Example 7: Comparison with true 3D algorithm
# If z-artifacts are unacceptable, use supervoxels instead
result_3d <- supervoxels(vec, mask, n_supvox = 500, alpha = 0.5)
# supervoxels provides smooth 3D parcels without slice artifacts
} # }
```
