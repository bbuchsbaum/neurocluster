#' SLiCE-MSF: Slice-wise, Low-rank, Minimum-Spanning Forest Clustering
#'
#' The slice_msf function performs spatially constrained clustering on a \code{NeuroVec} instance
#' using the SLiCE-MSF algorithm. This method uses temporal sketching via DCT basis functions,
#' split-half reliability weighting, and Felzenszwalb-Huttenlocher graph segmentation.
#'
#' @param vec A \code{NeuroVec} or \code{SparseNeuroVec} instance supplying the time series data to cluster.
#' @param mask A \code{NeuroVol} mask defining the voxels to include in the clustering result.
#' If the mask contains \code{numeric} data, nonzero values will define the included voxels.
#' If the mask is a \code{\linkS4class{LogicalNeuroVol}}, then \code{TRUE} will define the set
#' of included voxels.
#' @param target_k_global Target number of clusters across entire volume (exact if positive, uses RAG agglomeration).
#' Default is -1 (no target, uses natural FH clustering).
#' @param target_k_per_slice Target number of clusters per slice (exact if positive, ignored if stitch_z=TRUE).
#' Default is -1 (no target).
#' @param r DCT sketch rank (number of basis functions, excluding DC). Default is 12.
#' @param compactness A numeric value controlling the compactness of the clusters, with larger values resulting
#' in more compact clusters. Internally mapped to the scale parameter k. Default is 5.
#' @param min_size Minimum cluster size in voxels. Smaller clusters are merged. Default is 80.
#' @param num_runs Number of independent segmentation runs. If > 1, consensus fusion is applied. Default is 3.
#' @param consensus Logical. If TRUE and num_runs > 1, apply consensus fusion. Default is TRUE.
#' @param stitch_z Logical. If TRUE, stitch 2D slice clusters into 3D clusters. Default is TRUE.
#' @param theta_link Centroid similarity threshold for cross-slice stitching (0-1). Default is 0.85.
#' @param min_contact Minimum touching voxels between slices to attempt stitching. Default is 1.
#' @param nbhd Neighborhood connectivity (4, 6, or 8). Default is 8.
#' @param gamma Reliability weighting exponent. Higher values emphasize reliable voxels. Default is 1.5.
#' @param k_fuse Scale parameter for consensus fusion. If NULL, uses same as k. Default is NULL.
#' @param min_size_fuse Minimum cluster size for consensus. If NULL, uses min_size. Default is NULL.
#' @param use_features Use feature similarity in consensus fusion. Default is FALSE.
#' @param lambda Mix parameter for consensus (0-1). Higher values weight label agreement. Default is 0.7.
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
#' @examples
#' \dontrun{
#' mask <- NeuroVol(array(1, c(64,64,32)), NeuroSpace(c(64,64,32)))
#' vec <- replicate(100, NeuroVol(array(rnorm(64*64*32), c(64,64,32)),
#'                  NeuroSpace(c(64,64,32))), simplify=FALSE)
#' vec <- do.call(concat, vec)
#' 
#' # Single run
#' result <- slice_msf(vec, mask, K=200, num_runs=1)
#' 
#' # Multi-run consensus
#' result <- slice_msf(vec, mask, K=200, num_runs=3, consensus=TRUE)
#' }
#' 
#' # Using exact K targeting
#' result <- slice_msf(vec, mask, target_k_global=100, use_features=TRUE)
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
                      lambda = 0.7) {
  
  # Validate inputs
  assert_that(inherits(vec, "NeuroVec") || inherits(vec, "SparseNeuroVec"))
  assert_that(inherits(mask, "NeuroVol"))
  assert_that(nbhd %in% c(4, 6, 8))
  assert_that(r >= 1)
  assert_that(num_runs >= 1)
  
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
      target_k_per_slice = target_k_per_slice
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
        target_k_per_slice = target_k_per_slice
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
  centers_info <- compute_cluster_centers(vec, mask, cluster_ids, mask.idx)
  
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
#' Lower-level function that performs a single run of SLiCE-MSF segmentation.
#' Most users should use \code{slice_msf} instead.
#'
#' @param vec A \code{NeuroVec} instance supplying the time series data.
#' @param mask A \code{NeuroVol} mask defining the voxels to include.
#' @param r DCT sketch rank. Default is 12.
#' @param k Scale parameter (0-2). Smaller values create more clusters. Default is 0.32.
#' @param min_size Minimum cluster size. Default is 80.
#' @param nbhd Neighborhood connectivity (4, 6, or 8). Default is 8.
#' @param stitch_z Enable cross-slice stitching. Default is TRUE.
#' @param theta_link Centroid similarity threshold for stitching. Default is 0.85.
#' @param min_contact Minimum contact voxels for stitching. Default is 1.
#' @param gamma Reliability weighting exponent. Default is 1.5.
#'
#' @return A list with labels, weights, and sketch matrices.
#' @export
slice_msf_single <- function(vec, mask, 
                            r = 12,
                            k = 0.32,
                            min_size = 80,
                            nbhd = 8,
                            stitch_z = TRUE,
                            theta_link = 0.85,
                            min_contact = 1,
                            gamma = 1.5) {
  
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
    spatial_beta = 0.0
  )
  
  result
}

#' Consensus Fusion for SLiCE-MSF
#'
#' Combines multiple SLiCE-MSF segmentation runs using consensus clustering.
#'
#' @param run_results List of results from \code{slice_msf_single}.
#' @param mask A \code{NeuroVol} mask used in the original segmentation.
#' @param nbhd Neighborhood connectivity (4, 6, or 8). Default is 8.
#' @param k_fuse Scale parameter for fusion. Default is 0.30.
#' @param min_size_fuse Minimum cluster size. Default is 80.
#' @param use_features Use feature similarity in fusion. Default is FALSE.
#' @param lambda Mix parameter (0-1). Default is 0.7.
#'
#' @return A list with fused labels.
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
compute_cluster_centers <- function(vec, mask, cluster_ids, mask.idx) {
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