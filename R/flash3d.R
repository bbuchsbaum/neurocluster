#' FLASH-3D: Fast 3D Superclustering for fMRI
#'
#' Performs spatially-constrained clustering using Fast Low-rank Approximate 
#' Superclusters for Hemodynamics (FLASH-3D). This algorithm uses DCT-based 
#' temporal hashing and 3D jump-flood propagation for efficient clustering.
#' 
#' @note Consider using \code{\link{cluster4d}} with \code{method = "flash3d"} for a 
#' standardized interface across all clustering methods.
#'
#' @param vec A \code{NeuroVec} instance supplying the 4D data to cluster
#' @param mask A \code{NeuroVol} mask defining the voxels to include in the clustering result.
#'   If the mask contains \code{numeric} data, finite values strictly greater
#'   than zero define the included voxels.
#'   If the mask is a \code{\link[neuroim2:LogicalNeuroVol-class]{LogicalNeuroVol}}, then \code{TRUE} will define the set
#'   of included voxels.
#' @param K The number of clusters to find
#' @param lambda_s Finite non-negative spatial weight (default 0.6). It is
#'   annealed upward over rounds.
#' @param lambda_t Finite non-negative temporal weight for Hamming distance of
#'   time-series hashes (default 1.0).
#' @param lambda_g Finite non-negative boundary-field weight (default 0.0). If
#'   positive, \code{barrier} is required.
#' @param rounds Number of outer rounds (seed→flood→recenter). Default 2, typically 2-3 is sufficient.
#' @param bits Length of temporal hash (64 or 128). Default 64.
#' @param dctM Number of low DCT coefficients to rank-hash (default 12, range 4-32)
#' @param vox_scale Voxel size scaling for spatial distance, e.g., c(dx,dy,dz). Default c(1,1,1)
#' @param barrier Optional non-negative finite 3D boundary field with the same
#'   geometry as the mask. A candidate site is penalized when its seed and the
#'   target voxel have different boundary-field values.
#' @param verbose Logical indicating whether to print progress messages
#'
#' @return A \code{list} of class \code{cluster_result} with the following elements:
#' \describe{
#'   \item{clusvol}{An instance of type \link[neuroim2:ClusteredNeuroVol-class]{ClusteredNeuroVol}}
#'   \item{cluster}{An integer vector of cluster assignments for each voxel in mask}
#'   \item{centers}{A matrix of cluster centers in feature space (K x T)}
#'   \item{coord_centers}{A matrix of spatial coordinates of cluster centers (K x 3)}
#'   \item{K}{The number of clusters}
#'   \item{method}{A character string indicating the method used ("FLASH-3D")}
#' }
#'
#' @details
#' FLASH-3D uses a novel approach combining:
#' \itemize{
#'   \item DCT-based temporal feature hashing for fast similarity computation
#'   \item Jump-flood algorithm for efficient spatial propagation
#'   \item Blue-noise seeding for optimal initial cluster placement
#'   \item Annealing of spatial weights to encourage compact clusters
#' }
#'
#' The algorithm is particularly efficient for large-scale fMRI data, offering
#' significant speed improvements over iterative methods while maintaining
#' clustering quality.
#'
#' @examples
#' \dontrun{
#'   # Basic usage
#'   result <- supervoxels_flash3d(vec, mask, K = 100)
#'   
#'   # With custom parameters
#'   result <- supervoxels_flash3d(vec, mask, K = 100, 
#'                                 lambda_s = 0.8, lambda_t = 1.2,
#'                                 bits = 128, dctM = 16)
#'   
#'   # With barrier for anatomy-aware clustering
#'   barrier_vol <- create_anatomical_barrier(mask)
#'   result <- supervoxels_flash3d(vec, mask, K = 100,
#'                                 lambda_g = 0.5, barrier = barrier_vol)
#' }
#'
#' @references
#' FLASH-3D algorithm for fast superclustering of fMRI data (2025)
#'
#' @seealso \code{\link{supervoxels}}, \code{\link{snic}}, \code{\link{slic4d_supervoxels}}
#' @importFrom neuroim2 NeuroVec NeuroVol ClusteredNeuroVol series index_to_coord spacing
#' @export
supervoxels_flash3d <- function(vec, mask, K,
                               lambda_s = 0.6,
                               lambda_t = 1.0,
                               lambda_g = 0.0,
                               rounds = 2L,
                               bits = 64L,
                               dctM = 12L,
                               vox_scale = NULL,
                               barrier = NULL,
                               verbose = FALSE) {
  vec <- ensure_neurovec(vec)

  input <- validate_cluster4d_inputs(vec, mask, K, "supervoxels_flash3d")
  K <- input$n_clusters
  lambda_s <- .cluster4d_scalar_number(lambda_s, "lambda_s", "supervoxels_flash3d")
  lambda_t <- .cluster4d_scalar_number(lambda_t, "lambda_t", "supervoxels_flash3d")
  lambda_g <- .cluster4d_scalar_number(lambda_g, "lambda_g", "supervoxels_flash3d")
  if (any(c(lambda_s, lambda_t, lambda_g) < 0)) {
    stop("supervoxels_flash3d: lambda weights must be non-negative", call. = FALSE)
  }
  rounds <- .cluster4d_scalar_number(
    rounds, "rounds", "supervoxels_flash3d", integer = TRUE
  )
  bits <- .cluster4d_scalar_number(
    bits, "bits", "supervoxels_flash3d", integer = TRUE
  )
  dctM <- .cluster4d_scalar_number(
    dctM, "dctM", "supervoxels_flash3d", integer = TRUE
  )
  if (!bits %in% c(64L, 128L)) {
    stop("supervoxels_flash3d: bits must be 64 or 128", call. = FALSE)
  }
  if (dctM < 4L || dctM > 32L) {
    stop("supervoxels_flash3d: dctM must be between 4 and 32", call. = FALSE)
  }
  if (rounds < 1L) {
    stop("supervoxels_flash3d: rounds must be positive", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("supervoxels_flash3d: verbose must be TRUE or FALSE", call. = FALSE)
  }
  
  # Get dimensions
  dims <- dim(vec)
  nx <- dims[1]
  ny <- dims[2] 
  nz <- dims[3]
  ntime <- dims[4]
  
  # Handle voxel scaling
  if (is.null(vox_scale)) {
    vox_scale <- spacing(mask)
    if (length(vox_scale) != 3) {
      vox_scale <- c(1, 1, 1)
    }
  }
  if (!is.numeric(vox_scale) || length(vox_scale) != 3L ||
      any(!is.finite(vox_scale)) || any(vox_scale <= 0)) {
    stop(
      "supervoxels_flash3d: vox_scale must contain three finite positive values",
      call. = FALSE
    )
  }
  
  # Get mask indices and extract time series
  mask_idx <- input$mask_idx
  nmask <- length(mask_idx)

  # Validate the boundary field even for trivial K=N output so no supplied
  # parameter becomes silently inactive on a special-case path.
  barrier_vec <- NULL
  if (!is.null(barrier)) {
    if (is.array(barrier)) {
      if (!identical(as.integer(dim(barrier)), as.integer(c(nx, ny, nz)))) {
        stop(sprintf(
          "barrier dimensions (%s) must match mask dimensions (%d, %d, %d)",
          paste(dim(barrier), collapse = ", "), nx, ny, nz
        ), call. = FALSE)
      }
      barrier_vec <- as.numeric(barrier)
    } else if (inherits(barrier, "NeuroVol")) {
      barrier_geometry <- .cluster4d_geometry(barrier)
      if (!identical(barrier_geometry$dimensions, input$geometry$dimensions) ||
          !.cluster4d_same_numeric(barrier_geometry$spacing, input$geometry$spacing) ||
          !.cluster4d_same_numeric(barrier_geometry$affine, input$geometry$affine)) {
        stop("supervoxels_flash3d: barrier NeuroVol geometry must match mask", call. = FALSE)
      }
      barrier_vec <- as.numeric(as.array(barrier))
    } else {
      stop("supervoxels_flash3d: barrier must be an array or NeuroVol", call. = FALSE)
    }
    if (any(!is.finite(barrier_vec)) || any(barrier_vec < 0)) {
      stop(
        "supervoxels_flash3d: barrier values must be finite and non-negative",
        call. = FALSE
      )
    }
  } else if (lambda_g > 0) {
    stop("supervoxels_flash3d: positive lambda_g requires barrier", call. = FALSE)
  }
  
  if (K == nmask) {
    warning("K equals number of voxels, returning trivial clustering")
    labels <- rep(1:K, length.out = nmask)
    clusvol <- suppressWarnings(ClusteredNeuroVol(mask, labels))
    coords <- index_to_coord(mask, mask_idx)
    ts_matrix <- as.matrix(series(vec, mask_idx))
    if (ncol(ts_matrix) != nmask) ts_matrix <- matrix(ts_matrix, ncol = nmask)
    
    # Compute proper centers for each cluster
    centers <- matrix(0, nrow = K, ncol = nrow(ts_matrix))
    coord_centers <- matrix(0, nrow = K, ncol = 3)
    
    for (k in 1:K) {
      cluster_mask <- labels == k
      if (sum(cluster_mask) == 1) {
        centers[k, ] <- ts_matrix[, cluster_mask]
        coord_centers[k, ] <- coords[cluster_mask, ]
      } else if (sum(cluster_mask) > 1) {
        centers[k, ] <- rowMeans(ts_matrix[, cluster_mask, drop = FALSE])
        coord_centers[k, ] <- colMeans(coords[cluster_mask, , drop = FALSE])
      }
    }
    
    result <- structure(
      list(
        clusvol = clusvol,
        cluster = labels,
        centers = centers,
        coord_centers = coord_centers,
        K = K,
        n_clusters = K,
        method = "flash3d",
        parameters = list(
          n_clusters_requested = K,
          actual_K = K,
          lambda_s = lambda_s,
          lambda_t = lambda_t,
          lambda_g = lambda_g,
          rounds = rounds,
          bits = bits,
          dctM = dctM,
          dctM_effective = min(dctM, nrow(ts_matrix)),
          vox_scale = as.numeric(vox_scale)
        ),
        metadata = list(reseed_count = 0L)
      ),
      class = c("flash3d_result", "cluster_result", "list")
    )
    return(finalize_cluster4d_result(
      result, vec, mask, "flash3d", result$parameters
    ))
  }
  
  # Extract time series matrix - series returns T x Nmask which is what we need
  ts_matrix <- series(vec, mask_idx)
  if (!is.matrix(ts_matrix)) {
    ts_matrix <- matrix(as.numeric(ts_matrix), nrow = 1L, ncol = nmask)
  }
  dctM_eff <- min(dctM, nrow(ts_matrix), 32L)
  
  # Debug: check dimensions
  if (verbose) {
    cat("ts_matrix dimensions:", dim(ts_matrix), "\n")
    cat("mask_idx length:", length(mask_idx), "\n")
  }
  
  # Call C++ implementation (now returns List with labels, centers, coords)
  cpp_result <- flash3d_supervoxels_cpp(
    ts = ts_matrix,
    mask_lin0 = as.integer(mask_idx),  # R 1-based indices (C++ will convert)
    dims = as.integer(c(nx, ny, nz)),
    K = as.integer(K),
    lambda = c(lambda_s, lambda_t, lambda_g),
    rounds = as.integer(rounds),
    bits = as.integer(bits),
    dctM = as.integer(dctM),
    vox_scale = as.numeric(vox_scale),
    barrier_opt = barrier_vec,
    verbose = verbose
  )

  # Extract results from C++
  labels_mask <- cpp_result$labels
  centers <- cpp_result$centers      # K x T matrix (already computed in C++)
  coord_centers <- cpp_result$coords # K x 3 matrix (already computed in C++)
  
  used <- sort(unique(labels_mask[!is.na(labels_mask)]))
  if (!identical(as.integer(used), seq_len(K))) {
    stop(sprintf(
      "supervoxels_flash3d: native core returned %d clusters; exactly %d were requested",
      length(used), K
    ), call. = FALSE)
  }
  labels_mask <- as.integer(labels_mask)
  n_clusters <- K

  clusvol <- suppressWarnings(ClusteredNeuroVol(mask, labels_mask))

  result <- structure(
    list(
      clusvol = clusvol,
      cluster = labels_mask,
      centers = centers,           # Already computed in C++ (K_eff x T)
      coord_centers = coord_centers, # Already computed in C++ (K_eff x 3)
      K = n_clusters,
      n_clusters = n_clusters,
      method = "flash3d",
      parameters = list(
        n_clusters_requested = K,
        requested_K = K,
        actual_K = n_clusters,
        lambda_s = lambda_s,
        lambda_t = lambda_t,
        lambda_g = lambda_g,
        rounds = rounds,
        bits = bits,
        dctM = dctM,
        dctM_effective = dctM_eff,
        vox_scale = vox_scale
      ),
      metadata = list(
        cpp_coord_centers_voxel = coord_centers,
        reseed_count = cpp_result$reseed_count,
        round_reseed_counts = cpp_result$round_reseed_counts,
        reseed_indices_grid = cpp_result$reseed_indices_grid
      )
    ),
    class = c("cluster_result", "flash3d_result", "list")
  )

  finalize_cluster4d_result(result, vec, mask, "flash3d", result$parameters)
}
