#' SLiCE-MSF: Slice-wise, Low-rank, Minimum-Spanning Forest Clustering
#'
#' Performs spatially constrained clustering on neuroimaging time series data using a 
#' slice-based approach with optional 3D stitching. The algorithm applies DCT sketching
#' for temporal compression, reliability weighting, and graph-based segmentation. While
#' computationally efficient, the slice-based approach may create visible boundaries
#' between z-slices (see Details for mitigation strategies).
#'
#' @param vec A \code{NeuroVec} or \code{SparseNeuroVec} instance supplying the time series data to cluster.
#' @param mask A \code{NeuroVol} mask defining the voxels to include in the clustering result.
#' If the mask contains \code{numeric} data, nonzero values will define the included voxels.
#' If the mask is a \code{\linkS4class{LogicalNeuroVol}}, then \code{TRUE} will define the set
#' of included voxels.
#' @param target_k_global Target number of clusters across entire volume. When positive,
#'   uses region adjacency graph (RAG) agglomeration to achieve exactly K clusters.
#'   Default is -1 (no target, uses natural Felzenszwalb-Huttenlocher clustering).
#' @param target_k_per_slice Target number of clusters per slice. Only used when positive
#'   and stitch_z=FALSE. Useful for consistent parcellation across slices. Default is -1.
#' @param r DCT sketch rank (number of basis functions, excluding DC component). Higher
#'   values preserve more temporal detail but increase computation. Range: 4-20, Default: 12.
#' @param compactness Controls spatial compactness vs feature similarity balance (1-10).
#'   Lower values (1-3): Feature-driven, may create irregular shapes but better cross-slice
#'   alignment. Medium (4-6): Balanced. Higher (7-10): Spatially compact but may show more
#'   z-artifacts. Default is 5.
#' @param min_size Minimum cluster size in voxels. Smaller clusters are merged with
#'   nearest neighbors. Affects granularity of parcellation. Default is 80.
#' @param num_runs Number of independent segmentation runs. Single run (1) is faster but
#'   less stable. Multiple runs (3-5) with consensus fusion reduce variability and can
#'   smooth z-transitions. More than 5 has diminishing returns. Default is 3.
#' @param consensus Logical. If TRUE and num_runs > 1, applies consensus fusion across
#'   runs, improving stability and potentially reducing z-artifacts. Default is TRUE.
#' @param stitch_z Logical. If TRUE, attempts to stitch 2D slice clusters into coherent
#'   3D clusters by merging across z-boundaries. Essential for 3D continuity. Default is TRUE.
#' @param theta_link Centroid correlation threshold for cross-slice stitching (0-1).
#'   Lower values (0.70-0.80): Aggressive stitching, reduces z-plane artifacts but may
#'   over-merge. Default (0.85): Balanced. Higher (0.90-0.95): Conservative, preserves
#'   boundaries but more z-artifacts. Critical parameter for z-continuity.
#' @param min_contact Minimum number of touching voxels between slices required for
#'   stitching attempt. Lower (1-2): More connections, better continuity. Higher (3-5):
#'   Stricter requirement, prevents spurious bridges. Default is 1.
#' @param nbhd Neighborhood connectivity for within-slice clustering. Options: 4 (von
#'   Neumann), 6 (includes z but mapped to 8), 8 (Moore). Higher connectivity can
#'   improve within-slice coherence. Default is 8.
#' @param gamma Reliability weighting exponent for split-half correlation. Higher values
#'   (>1.5) emphasize high-reliability voxels, useful for noisy data. Lower values (<1)
#'   treat all voxels more equally. Default is 1.5.
#' @param k_fuse Scale parameter for consensus fusion graph. If NULL, uses same as
#'   compactness-derived scale. Lower values create more clusters in fusion. Default is NULL.
#' @param min_size_fuse Minimum cluster size during consensus fusion. If NULL, uses
#'   min_size. Can be set lower to preserve small consistent regions. Default is NULL.
#' @param use_features Include feature similarity in consensus fusion. When TRUE, uses
#'   both label agreement and temporal similarity, improving cross-slice consistency.
#'   Recommended when targeting exact K. Default is FALSE.
#' @param lambda Mixing parameter for consensus (0-1). Controls balance between label
#'   agreement and feature similarity when use_features=TRUE. Higher values weight
#'   label agreement more. Default is 0.7.
#' @param z_mult Z-smoothing factor (0-1). Values > 0 softly blend DCT sketches between
#'   adjacent slices before clustering, reducing visible z-plane seams. 0 preserves the
#'   legacy per-slice behavior. Recommended range 0.1-0.4. Default is 0.0.
#'
#' @return A \code{list} of class \code{slice_msf_cluster_result} with the following elements:
#' \describe{
#' \item{clusvol}{An instance of type \linkS4class{ClusteredNeuroVol}.}
#' \item{cluster}{A vector of cluster indices equal to the number of voxels in the mask.}
#' \item{centers}{A matrix of cluster centers with each column representing the feature vector for a cluster.}
#' \item{coord_centers}{A matrix of spatial coordinates with each row corresponding to a cluster.}
#' \item{runs}{If num_runs > 1, a list of individual run results.}
#' }
#' 
#' @details
#' ## Algorithm Overview
#' 
#' SLiCE-MSF (Slice-wise, Low-rank, Minimum-Spanning Forest) is designed for fast
#' clustering of high-resolution fMRI data. The algorithm proceeds in stages:
#' 
#' 1. **Temporal Sketching**: Each voxel's time series is compressed using Discrete 
#'    Cosine Transform (DCT) basis functions, reducing from T timepoints to r coefficients.
#' 
#' 2. **Reliability Weighting**: Split-half correlations are computed to identify
#'    reliable voxels, which receive higher weight in clustering decisions.
#' 
#' 3. **Slice-wise Clustering**: Each axial slice is clustered independently using
#'    Felzenszwalb-Huttenlocher graph segmentation, which efficiently finds regions
#'    with high internal similarity and low external similarity.
#' 
#' 4. **Z-Stitching** (optional): Clusters from adjacent slices are merged based on
#'    spatial contact and centroid similarity to create 3D parcels.
#' 
#' 5. **Consensus Fusion** (optional): Multiple runs are combined using co-association
#'    matrices to improve stability.
#' 
#' ## Z-Plane Artifacts and Mitigation
#' 
#' Because clustering is performed slice-by-slice, the algorithm can produce visible
#' horizontal lines when viewed in sagittal or coronal planes. These artifacts are
#' inherent to the slice-based approach but can be minimized:
#' 
#' ### Artifact Reduction Strategies
#' 
#' **Mild artifacts** (slight discontinuities):
#' - Reduce theta_link to 0.75-0.80 for more aggressive stitching
#' - Set min_contact to 2-3 for balanced connectivity
#' - Use lower compactness (2-4) for more flexible cluster shapes
#' - Set z_mult between 0.1-0.3 to softly bleed information across slices
#' 
#' **Moderate artifacts** (visible lines):
#' - Use multiple runs (num_runs = 3-5) with consensus = TRUE
#' - Enable use_features = TRUE for feature-based consensus
#' - Adjust lambda to 0.5-0.6 to weight features more
#' - Combine with z_mult smoothing for additional continuity without heavy fusion
#' 
#' **Severe artifacts** (strong discontinuities):
#' - Consider alternative algorithms like supervoxels() or snic() for true 3D clustering
#' - These provide smoother 3D parcels at the cost of increased computation time
#' 
#' ## Parameter Selection Guidelines
#' 
#' ### For Whole-Brain Parcellation
#' - r = 10-15 (balance detail and speed)
#' - compactness = 4-6 (balanced spatial/feature weighting)
#' - min_size = 80-150 (appropriate for ~3mm resolution)
#' - num_runs = 3-5 with consensus = TRUE
#' 
#' ### For ROI Analysis
#' - r = 15-20 (preserve more temporal detail)
#' - compactness = 2-4 (feature-driven clustering)
#' - min_size = 20-50 (allow smaller parcels)
#' - theta_link = 0.75 (aggressive stitching within ROI)
#' 
#' ### For High-Resolution Data (< 2mm)
#' - Increase min_size proportionally (e.g., 200-300 voxels)
#' - Use target_k_global to control final parcel count
#' - Consider target_k_per_slice for consistent slice-wise parcellation
#' 
#' ## Performance Considerations
#' 
#' - **Memory**: Scales with mask size × r × num_runs
#' - **Speed**: Much faster than full 3D methods, especially for high-resolution data
#' - **Trade-off**: Speed vs. z-continuity - for smooth 3D parcels, use supervoxels()
#' 
#' ## Parallelization Strategy
#' 
#' SLiCE-MSF uses **automatic parallelization via RcppParallel** for key operations:
#' 
#' ### Parallel Operations:
#' 1. **DCT Sketching** (`SliceSketchWorker`): Each z-slice is processed in parallel
#'    - Detrending and z-scoring of time series
#'    - Split-half reliability computation
#'    - DCT coefficient calculation
#'    - Threads automatically assigned by RcppParallel based on available cores
#' 
#' 2. **Consensus Fusion** (`FuseSliceWorker`): When num_runs > 1
#'    - Co-association matrix computation parallelized across slices
#'    - Edge weight calculations done in parallel
#'    - Graph segmentation remains sequential within each slice
#' 
#' ### Performance Characteristics:
#' - **Scaling**: Near-linear speedup with cores for the sketching phase
#' - **Optimal for**: High-resolution data (many slices to process)
#' - **Automatic**: No user configuration needed - uses all available cores
#' - **Memory-efficient**: Each thread processes independent slices
#' - **Bottleneck**: Graph segmentation phase is sequential per slice
#' 
#' ### Controlling Parallelization:
#' RcppParallel automatically determines the number of threads. To control:
#' ```r
#' # Set number of threads globally
#' RcppParallel::setThreadOptions(numThreads = 4)
#' 
#' # Reset to automatic
#' RcppParallel::setThreadOptions(numThreads = "auto")
#' ```
#' 
#' @section Troubleshooting Z-Plane Artifacts:
#' 
#' Common issues and solutions:
#' 
#' **Horizontal lines in sagittal/coronal views:**
#' 1. First, verify artifacts aren't anatomically meaningful (e.g., at gray/white boundaries)
#' 2. Try progressive adjustments:
#'    - Step 1: theta_link = 0.75 (from default 0.85)
#'    - Step 2: Add num_runs = 3, consensus = TRUE
#'    - Step 3: Set use_features = TRUE, lambda = 0.5
#'    - Step 4: Reduce compactness to 2-3
#' 
#' **Over-merging after reducing theta_link:**
#' - Increase min_contact to 3-5 to require more evidence for merging
#' - Increase min_size to prevent small bridging clusters
#' 
#' **Inconsistent cluster sizes across slices:**
#' - Use target_k_per_slice (with stitch_z = FALSE) for consistent slice parcellation
#' - Or use target_k_global for consistent total cluster count
#' 
#' **Poor performance in low-SNR regions:**
#' - Increase gamma to 2.0-2.5 to emphasize reliable voxels
#' - Consider masking out low-SNR regions before clustering
#'
#' @examples
#' \dontrun{
#' # Load example data
#' mask <- NeuroVol(array(1, c(64,64,32)), NeuroSpace(c(64,64,32)))
#' vec <- replicate(100, NeuroVol(array(rnorm(64*64*32), c(64,64,32)),
#'                  NeuroSpace(c(64,64,32))), simplify=FALSE)
#' vec <- do.call(concat, vec)
#' 
#' # Example 1: Basic usage (may show z-plane artifacts)
#' result_basic <- slice_msf(vec, mask, 
#'                           num_runs = 1,
#'                           compactness = 5)
#' 
#' # Example 2: Reduced z-artifacts with aggressive stitching
#' # Recommended for minimizing slice boundaries
#' result_smooth <- slice_msf(vec, mask, 
#'                            theta_link = 0.75,      # More aggressive stitching
#'                            min_contact = 2,        # Moderate contact requirement
#'                            compactness = 3,        # Lower compactness
#'                            num_runs = 3,           # Multiple runs
#'                            consensus = TRUE,       # Enable consensus
#'                            use_features = TRUE,    # Feature-based consensus
#'                            lambda = 0.6)           # Balance features/labels
#'
#' # Example 3: Conservative approach for anatomical boundaries
#' # Preserves natural boundaries but may show more z-artifacts
#' result_conservative <- slice_msf(vec, mask,
#'                                  theta_link = 0.90,    # Conservative stitching
#'                                  min_contact = 5,      # Strict contact requirement
#'                                  compactness = 7,      # Compact clusters
#'                                  min_size = 120)       # Larger minimum size
#'
#' # Example 4: Exact K targeting with 100 clusters
#' result_exact <- slice_msf(vec, mask, 
#'                          target_k_global = 100,   # Exactly 100 clusters
#'                          use_features = TRUE,     # Required for exact K
#'                          num_runs = 3,
#'                          consensus = TRUE)
#'
#' # Example 5: Per-slice consistency (useful for group studies)
#' result_per_slice <- slice_msf(vec, mask,
#'                               target_k_per_slice = 50,  # 50 clusters per slice
#'                               stitch_z = FALSE,         # No z-stitching
#'                               compactness = 5)
#'
#' # Example 6: High-resolution data optimization
#' # For data with voxel size < 2mm
#' result_highres <- slice_msf(vec, mask,
#'                             min_size = 250,         # Larger clusters for high-res
#'                             r = 15,                 # More DCT components
#'                             compactness = 4,
#'                             theta_link = 0.78,
#'                             num_runs = 5,
#'                             consensus = TRUE)
#'
#' # Example 7: Comparison with true 3D algorithm
#' # If z-artifacts are unacceptable, use supervoxels instead
#' result_3d <- supervoxels(vec, mask, n_supvox = 500, alpha = 0.5)
#' # supervoxels provides smooth 3D parcels without slice artifacts
#' }
#'
#' @references
#' Felzenszwalb, P. F., & Huttenlocher, D. P. (2004). Efficient graph-based image segmentation.
#' International journal of computer vision, 59(2), 167-181.
#'
#' @importFrom neuroim2 NeuroVec NeuroVol series spacing index_to_coord ClusteredNeuroVol concat
#' @importFrom Matrix t
#' @import assertthat
#' @export
slice_msf <- function(vec, mask, 
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
                      z_mult = 0.0) {
  
  # Basic parameter validation
  assert_that(nbhd %in% c(4, 6, 8))
  assert_that(r >= 1)
  assert_that(num_runs >= 1)
  assert_that(is.numeric(z_mult), length(z_mult) == 1, z_mult >= 0, z_mult <= 1)
  
  # Check spatial dimension compatibility
  vec_dims <- dim(vec)[1:3]  # First 3 dimensions are spatial (x, y, z)
  mask_dims <- dim(mask)
  if (!identical(vec_dims, mask_dims)) {
    stop("NeuroVec and mask must have identical spatial dimensions. ",
         "NeuroVec dimensions: ", paste(vec_dims, collapse="x"), 
         ", mask dimensions: ", paste(mask_dims, collapse="x"))
  }
  
  # Map nbhd=6 to nbhd=8 (C++ implementation only supports 4 or 8)
  if (nbhd == 6) {
    message("Note: nbhd=6 is mapped to nbhd=8 for compatibility")
    nbhd <- 8
  }
  
  # Extract data
  mask.idx <- which(mask > 0)
  if (length(mask.idx) == 0) {
    stop("No valid voxels found in mask. Mask must contain at least one positive value.")
  }
  effective_k <- if (target_k_global > 0) target_k_global else min(100, length(mask.idx))
  validate_cluster4d_inputs(vec, mask, effective_k, "slice_msf")
  
  # Get time series matrix for all voxels
  # The C++ function expects data for the full volume
  all_idx <- seq_len(prod(dim(mask)))
  # series() returns T x N when vec is a NeuroVec
  TS <- series(vec, all_idx)
  
  # Prepare mask and dimensions
  mask_flat <- as.integer(mask@.Data)
  vol_dim <- dim(mask)
  
  # Map compactness to fh_scale parameter (inverse relationship)
  # Higher compactness -> smaller fh_scale -> more compact clusters
  fh_scale <- 2.0 / (compactness + 1.0)  # Maps compactness 5 -> fh_scale ≈ 0.33
  
  # Get voxel dimensions
  voxel_dim <- spacing(mask)
  
  # Set fusion parameters
  if (is.null(k_fuse)) k_fuse <- fh_scale
  if (is.null(min_size_fuse)) min_size_fuse <- min_size
  
  # Single run or multiple runs
  if (num_runs == 1 || !consensus) {
    # Single run
    result <- slice_msf_runwise(
      TS = TS,
      mask = mask_flat,
      vol_dim = vol_dim,
      r = r,
      fh_scale = fh_scale,
      min_size = min_size,
      nbhd = nbhd,
      stitch_z = stitch_z,
      theta_link = theta_link,
      min_contact = min_contact,
      rows_are_time = TRUE,
      gamma = gamma,
      voxel_dim = voxel_dim,
      spatial_beta = 0.0,
      target_k_global = target_k_global,
      target_k_per_slice = target_k_per_slice,
      z_mult = z_mult
    )
    
    labels <- result$labels
    
  } else {
    # Multiple runs with consensus
    runs <- vector("list", num_runs)
    
    for (i in seq_len(num_runs)) {
      runs[[i]] <- slice_msf_runwise(
        TS = TS,
        mask = mask_flat,
        vol_dim = vol_dim,
        r = r,
        fh_scale = fh_scale,
        min_size = min_size,
        nbhd = nbhd,
        stitch_z = stitch_z,
        theta_link = theta_link,
        min_contact = min_contact,
        rows_are_time = TRUE,
        gamma = gamma,
        voxel_dim = voxel_dim,
        spatial_beta = 0.0,
        target_k_global = target_k_global,
        target_k_per_slice = target_k_per_slice,
        z_mult = z_mult
      )
    }
    
    # Consensus fusion
    consensus_result <- slice_fuse_consensus(
      run_results = runs,
      vol_dim = vol_dim,
      nbhd = nbhd,
      fh_scale = k_fuse,
      min_size = min_size_fuse,
      use_features = use_features,
      lambda = lambda,
      voxel_dim = voxel_dim,
      spatial_beta = 0.0,
      target_k_global = target_k_global,
      target_k_per_slice = target_k_per_slice,
      stitch_z = stitch_z
    )
    
    labels <- consensus_result$labels
    result <- list(labels = labels, runs = runs)
  }
  
  # Extract labels for mask voxels only
  cluster_ids <- labels[mask.idx]
  
  # Compute cluster centers and spatial centroids
  centers_info <- compute_slice_cluster_centers(vec, mask, cluster_ids, mask.idx)
  
  # Create ClusteredNeuroVol with consistent logical mask (only positive values are TRUE)
  logical_mask <- mask > 0
  kvol <- ClusteredNeuroVol(logical_mask, clusters = cluster_ids)
  
  # Prepare return structure
  ret <- structure(
    list(
      clusvol = kvol,
      cluster = cluster_ids,
      centers = centers_info$centers,
      coord_centers = centers_info$coord_centers
    ),
    class = c("slice_msf_cluster_result", "cluster_result", "list")
  )
  
  # Add run information if multiple runs
  if (num_runs > 1 && consensus) {
    ret$runs <- runs
  }
  
  ret
}

#' Single Run SLiCE-MSF Segmentation
#'
#' Lower-level function that performs a single run of SLiCE-MSF segmentation without
#' consensus fusion. This function is useful for testing parameters or when you need
#' direct access to the DCT sketch and reliability weights. Most users should use 
#' \code{slice_msf} instead for better stability through consensus.
#'
#' @param vec A \code{NeuroVec} instance supplying the time series data.
#' @param mask A \code{NeuroVol} mask defining the voxels to include.
#' @param r DCT sketch rank (number of basis functions). Higher values preserve more
#'   temporal detail but increase computation. Range: 4-20. Default is 12.
#' @param k Scale parameter for Felzenszwalb-Huttenlocher segmentation (0-2). Controls
#'   the trade-off between under- and over-segmentation. Smaller values (0.1-0.3) create
#'   more clusters, larger values (0.5-1.0) create fewer, larger clusters. Default is 0.32.
#' @param min_size Minimum cluster size in voxels. Clusters smaller than this are merged
#'   with neighbors. Default is 80.
#' @param nbhd Neighborhood connectivity for within-slice edges (4 or 8). Note that 6 is
#'   automatically mapped to 8. Default is 8.
#' @param stitch_z Enable cross-slice stitching to create 3D clusters. When FALSE, each
#'   slice is treated independently. Default is TRUE.
#' @param theta_link Centroid correlation threshold for z-stitching (0-1). Lower values
#'   allow more aggressive merging across slices. Only used if stitch_z=TRUE. Default is 0.85.
#' @param min_contact Minimum number of vertically adjacent voxels required to consider
#'   stitching two clusters. Only used if stitch_z=TRUE. Default is 1.
#' @param gamma Reliability weighting exponent based on split-half correlations. Higher
#'   values give more weight to reliable voxels. Default is 1.5.
#' @param z_mult Z-smoothing factor (0-1). Values > 0 blend slice sketches prior to
#'   clustering to reduce boundary artifacts. Default is 0.0.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{labels}{Integer vector of cluster labels for each voxel}
#'   \item{weights}{Numeric vector of reliability weights for each voxel}
#'   \item{sketch}{Matrix of DCT coefficients (r × n_voxels)}
#'   \item{params}{List of parameters used in the segmentation}
#' }
#' 
#' @details
#' This function exposes the core SLiCE-MSF algorithm without consensus fusion. The
#' returned sketch matrix contains the DCT coefficients that summarize each voxel's
#' time series, and the weights indicate the reliability of each voxel based on
#' split-half correlations. These can be useful for diagnostic purposes or custom
#' post-processing.
#' 
#' The k parameter is particularly important: it directly controls the scale of
#' segmentation. For typical fMRI data, values between 0.2-0.5 work well. Lower
#' values are needed for fine-grained parcellation, while higher values produce
#' coarser segmentation.
#' @export
slice_msf_single <- function(vec, mask, 
                            r = 12,
                            k = 0.32,
                            min_size = 80,
                            nbhd = 8,
                            stitch_z = TRUE,
                            theta_link = 0.85,
                            min_contact = 1,
                            gamma = 1.5,
                            z_mult = 0.0) {
  
  # Validate inputs
  assert_that(inherits(vec, "NeuroVec") || inherits(vec, "SparseNeuroVec"))
  assert_that(inherits(mask, "NeuroVol"))
  
  # Check spatial dimension compatibility
  vec_dims <- dim(vec)[1:3]  # First 3 dimensions are spatial (x, y, z)
  mask_dims <- dim(mask)
  if (!identical(vec_dims, mask_dims)) {
    stop("NeuroVec and mask must have identical spatial dimensions. ",
         "NeuroVec dimensions: ", paste(vec_dims, collapse="x"), 
         ", mask dimensions: ", paste(mask_dims, collapse="x"))
  }
  
  # Map nbhd=6 to nbhd=8 (C++ implementation only supports 4 or 8)
  if (nbhd == 6) {
    message("Note: nbhd=6 is mapped to nbhd=8 for compatibility")
    nbhd <- 8
  }
  
  # Extract data
  mask.idx <- which(mask > 0)
  all_idx <- seq_len(prod(dim(mask)))
  # series() returns T x N when vec is a NeuroVec
  TS <- series(vec, all_idx)
  
  mask_flat <- as.integer(mask@.Data)
  vol_dim <- dim(mask)
  voxel_dim <- spacing(mask)
  assert_that(is.numeric(z_mult), length(z_mult) == 1, z_mult >= 0, z_mult <= 1)
  
  # Call C++ function
  result <- slice_msf_runwise(
    TS = TS,
    mask = mask_flat,
    vol_dim = vol_dim,
    r = r,
    fh_scale = k,
    min_size = min_size,
    nbhd = nbhd,
    stitch_z = stitch_z,
    theta_link = theta_link,
    min_contact = min_contact,
    rows_are_time = TRUE,
    gamma = gamma,
    voxel_dim = voxel_dim,
    spatial_beta = 0.0,
    z_mult = z_mult
  )
  
  result
}

#' Consensus Fusion for SLiCE-MSF
#'
#' Combines multiple SLiCE-MSF segmentation runs using consensus clustering to improve
#' stability and reduce variability. This function can help reduce z-plane artifacts
#' by averaging cluster assignments across multiple independent runs.
#'
#' @param run_results List of results from \code{slice_msf_single}, each containing
#'   labels, weights, and sketch matrices.
#' @param mask A \code{NeuroVol} mask used in the original segmentation.
#' @param nbhd Neighborhood connectivity (4, 6, or 8). Should match the value used
#'   in the original segmentation. Default is 8.
#' @param k_fuse Scale parameter for consensus fusion graph (0-2). Lower values create
#'   more refined clusters during fusion. If too low, may fragment consensus clusters.
#'   Default is 0.30.
#' @param min_size_fuse Minimum cluster size during fusion. Can be set lower than
#'   the original min_size to preserve small but consistent clusters. Default is 80.
#' @param use_features Include feature similarity in consensus fusion. When TRUE,
#'   uses both label co-occurrence and temporal similarity from DCT sketches. This
#'   can improve cross-slice consistency and is recommended when targeting exact K.
#'   Default is FALSE.
#' @param lambda Mixing parameter (0-1) controlling the balance between label agreement
#'   and feature similarity when use_features=TRUE. Higher values (0.7-0.9) emphasize
#'   label agreement, lower values (0.3-0.6) emphasize features. Default is 0.7.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{labels}{Integer vector of consensus cluster labels}
#'   \item{params}{List of parameters used in consensus fusion}
#' }
#' 
#' @details
#' Consensus fusion works by building a co-association matrix that counts how often
#' each pair of voxels is assigned to the same cluster across runs. This matrix is
#' then used as edge weights in a new graph segmentation. When use_features=TRUE,
#' the edge weights also incorporate the similarity of DCT sketches, providing
#' additional stability.
#' 
#' For reducing z-artifacts, using use_features=TRUE with lambda around 0.5-0.6 can
#' help ensure that clusters are consistent both in terms of their boundaries and
#' their temporal characteristics.
#' @export
slice_msf_consensus <- function(run_results, mask,
                               nbhd = 8,
                               k_fuse = 0.30,
                               min_size_fuse = 80,
                               use_features = FALSE,
                               lambda = 0.7) {
  
  vol_dim <- dim(mask)
  voxel_dim <- spacing(mask)
  
  # Call C++ consensus function
  result <- slice_fuse_consensus(
    run_results = run_results,
    vol_dim = vol_dim,
    nbhd = nbhd,
    k_fuse = k_fuse,
    min_size_fuse = min_size_fuse,
    use_features = use_features,
    lambda = lambda,
    voxel_dim = voxel_dim,
    spatial_beta = 0.0
  )
  
  result
}

# Helper function to compute cluster centers
compute_slice_cluster_centers <- function(vec, mask, cluster_ids, mask.idx) {
  unique_clusters <- sort(unique(cluster_ids[cluster_ids > 0]))
  n_clusters <- length(unique_clusters)
  
  if (n_clusters == 0) {
    return(list(centers = NULL, coord_centers = NULL))
  }
  
  # Get time series data
  vecmat <- series(vec, mask.idx)
  
  # Compute data centers (mean time series per cluster)
  centers <- matrix(0, nrow = nrow(vecmat), ncol = n_clusters)
  coord_centers <- matrix(0, nrow = n_clusters, ncol = 3)
  
  coords <- index_to_coord(mask, mask.idx)
  
  for (i in seq_along(unique_clusters)) {
    clust_id <- unique_clusters[i]
    idx <- which(cluster_ids == clust_id)
    
    if (length(idx) > 0) {
      # Mean time series
      centers[, i] <- rowMeans(vecmat[, idx, drop = FALSE])
      
      # Spatial centroid
      coord_centers[i, ] <- colMeans(coords[idx, , drop = FALSE])
    }
  }
  
  list(centers = centers, coord_centers = coord_centers)
}
