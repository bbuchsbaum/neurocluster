#' Fast 4D SLIC supervoxels (mask-aware, gradient relocation, preserve-K)
#'
#' Cluster a 4D \code{NeuroVec} (x,y,z,time) into compact, spatially contiguous
#' 3D supervoxels using an enhanced SLIC-style algorithm with mask-aware seeding,
#' gradient-based seed relocation, and exact K preservation.
#' 
#' @note Consider using \code{\link{cluster4d}} with \code{method = "slic"} for a
#' standardized interface across all clustering methods.
#'
#' @param bvec A \code{NeuroVec} with dims (X, Y, Z, T).
#'   Can also be a 3D \code{\linkS4class{NeuroVol}} for structural image segmentation,
#'   which will be automatically converted to a single-timepoint NeuroVec internally.
#' @param mask A 3D \code{NeuroVol} (or logical array) indicating voxels to include.
#' @param K Target number of supervoxels.
#' @param compactness Spatial vs feature tradeoff (like SLIC 'm'). Larger = more compact.
#' @param max_iter Maximum iterations (default 10).
#' @param n_threads Number of CPU threads to use (0 = auto).
#' @param step_mm Optional approximate spacing between seeds in millimeters; if NULL,
#'   computed as cubic root of bounding-box volume / K.
#' @param n_components If > 0, random-project each voxel's time series to this dimension
#'   for speed (Johnson-Lindenstrauss style). 0 = use full time series.
#' @param feature_norm One of "zscale", "l2", "none".
#' @param seed_method One of "mask_poisson" (Poisson disk in mask), "mask_grid" (grid in mask),
#'   "grid" (regular grid), "farthest" (farthest point sampling).
#' @param seed_relocate One of "correlation" (correlation gradient), "intensity" (mean intensity gradient),
#'   "spatial" (spatial gradient using neighborweights), "none" (no relocation).
#' @param seed_relocate_radius Search radius in voxels for gradient-based seed relocation (default 1).
#' @param connectivity Neighborhood connectivity: 6 (face neighbors) or 26 (all neighbors).
#' @param strict_connectivity Enforce exactly one connected component per label (default TRUE).
#' @param enforce_connectivity Alias for strict_connectivity (backward compatibility).
#' @param preserve_k Ensure exactly K non-empty labels by splitting largest clusters (default FALSE).
#' @param topup_iters Number of refinement iterations after splitting for preserve_k (default 2).
#' @param min_size Minimum component size (voxels) to keep before relabel (default 0).
#' @param verbose Logical.
#'
#' @return A list of class \code{cluster_result} with elements:
#' \itemize{
#'   \item \code{clusvol}: \code{ClusteredNeuroVol} with the final labels
#'   \item \code{cluster}: integer vector of length = #masked voxels
#'   \item \code{centers}: matrix (K x d_feat) center features
#'   \item \code{coord_centers}: matrix (K x 3) spatial centers in mm
#' }
#'
#' @examples
#' \dontrun{
#'   library(neuroim2)
#'   # Basic usage
#'   res <- slic4d_supervoxels(bvec, mask, K = 1000, compactness = 15)
#'   
#'   # With mask-aware seeding and gradient relocation
#'   res <- slic4d_supervoxels(bvec, mask, K = 1000, 
#'                            seed_method = "mask_poisson",
#'                            seed_relocate = "correlation",
#'                            preserve_k = TRUE)
#' }
#' @export
slic4d_supervoxels <- function(bvec, mask,
                              K,
                              compactness = 10,
                              max_iter = 10,
                              n_threads = 0,
                              step_mm = NULL,
                              n_components = 0,
                              feature_norm = c("zscale", "l2", "none"),
                              seed_method = c("mask_poisson", "mask_grid", "grid", "farthest"),
                              seed_relocate = c("none", "correlation", "intensity", "spatial"),
                              seed_relocate_radius = 1L,
                              connectivity = c(26L, 6L),
                              strict_connectivity = TRUE,
                              enforce_connectivity = TRUE,
                              preserve_k = FALSE,
                              topup_iters = 2L,
                              min_size = 0L,
                              verbose = FALSE) {
  
  # Allow single-volume NeuroVol by wrapping to NeuroVec
  bvec <- ensure_neurovec(bvec)

  feature_norm <- match.arg(feature_norm)
  seed_method <- match.arg(seed_method)
  seed_relocate <- match.arg(seed_relocate)
  connectivity <- as.integer(match.arg(as.character(connectivity), c("26", "6")))
  
  # Use common validation
  validate_cluster4d_inputs(bvec, mask, K, "slic4d_supervoxels")

  n_timepoints <- if (length(dim(bvec)) >= 4) dim(bvec)[4] else 1L
  if (n_timepoints <= 1 && seed_relocate == "correlation") {
    if (verbose) {
      message("slic4d_supervoxels: correlation relocation not applicable for single timepoint; using intensity gradient instead")
    }
    seed_relocate <- "intensity"
  }
  
  # Handle mask input
  if (inherits(mask, "NeuroVol")) {
    mask_arr <- as.logical(mask)
    sp <- neuroim2::space(mask)
  } else {
    mask_arr <- as.logical(mask)
    sp <- neuroim2::space(neuroim2::vol(mask_arr))
  }
  
  mask_idx <- which(mask_arr)
  if (length(mask_idx) == 0) stop("Mask is empty")
  
  # Get dimensions and voxel sizes
  dims <- dim(mask_arr)
  voxmm <- neuroim2::spacing(sp)
  
  # Spatial coordinates in mm for masked voxels
  coords <- neuroim2::index_to_coord(sp, mask_idx)
  coords <- as.matrix(coords)  # N x 3
  
  # Feature matrix: voxel time series (N x T)
  # series() returns T x N, so we need to transpose
  # Note: series() returns a vector for single-timepoint data
  feat <- neuroim2::series(bvec, mask_idx)
  if (!is.matrix(feat)) {
    feat <- matrix(feat, nrow = 1)  # 1 x N for single timepoint
  }
  feat <- t(as.matrix(feat))  # Now N x T
  
  # Optional feature normalization
  if (feature_norm == "zscale") {
    # Z-score each timepoint across voxels
    mu <- colMeans(feat)
    sdv <- apply(feat, 2, sd)
    sdv[sdv == 0] <- 1
    feat <- sweep(sweep(feat, 2, mu, "-"), 2, sdv, "/")
  } else if (feature_norm == "l2") {
    nrm <- sqrt(rowSums(feat * feat))
    nrm[nrm == 0] <- 1
    feat <- feat / nrm
  }
  
  # Optional random projection for speed
  if (n_components > 0 && n_components < ncol(feat)) {
    set.seed(123)
    R <- matrix(rnorm(ncol(feat) * n_components) / sqrt(n_components), 
                nrow = ncol(feat), ncol = n_components)
    feat <- feat %*% R  # N x n_components
  }
  
  # Compute default step from spatial extent / K if not supplied.
  # For thin volumes (e.g., single-slice 2D images), using full 3D volume can
  # make step_mm artificially tiny and collapse all voxels into one cluster.
  # Here we treat strongly anisotropic (quasi‑2D) data using an area‑based
  # heuristic and thicker volumes with the original 3D heuristic.
  if (is.null(step_mm)) {
    mins <- apply(coords, 2, min)
    maxs <- apply(coords, 2, max)
    extents <- pmax(1e-6, maxs - mins)
    # If z‑extent is much smaller than x/y, treat as 2D
    if (extents[3] < min(extents[1:2]) / 4) {
      area <- extents[1] * extents[2]
      step_mm <- sqrt(area / K)
    } else {
      vol <- prod(extents)
      step_mm <- (vol / K)^(1/3)
    }
  }
  
  # Compute gradient volume for seed relocation if needed
  grad_masked <- numeric(0)
  if (seed_relocate != "none") {
    if (seed_relocate == "correlation") {
      # Use correlation gradient for fMRI time series
      img4d <- as.array(bvec)
      # Convert mask to numeric while preserving dimensions
      mask_numeric <- array(as.numeric(mask_arr), dim = dims)
      grad3d <- correlation_gradient_cpp(img4d, mask_numeric)
      dim(grad3d) <- dims
      grad_masked <- grad3d[mask_idx]
    } else if (seed_relocate == "intensity") {
      # Use mean intensity gradient
      img4d <- as.array(bvec)
      mean3d <- apply(img4d, c(1,2,3), mean)
      grad3d <- .grad3d_fdiff(mean3d)
      grad_masked <- grad3d[mask_idx]
    } else if (seed_relocate == "spatial") {
      # Use spatial gradient via neighborweights
      if (!requireNamespace("neighborweights", quietly = TRUE)) {
        warning("neighborweights not installed; falling back to intensity gradient")
        img4d <- as.array(bvec)
        mean3d <- apply(img4d, c(1,2,3), mean)
        grad3d <- .grad3d_fdiff(mean3d)
      } else {
        grad3d <- spatial_gradient(
          neuroim2::vol(apply(as.array(bvec), c(1,2,3), mean), space = sp),
          neuroim2::vol(mask_arr, space = sp)
        )
        if (inherits(grad3d, "NeuroVol")) {
          grad3d <- as.array(grad3d)
        }
      }
      grad_masked <- grad3d[mask_idx]
    }
  }
  
  # Run core C++ implementation
  core <- slic4d_core(
    feat, coords,
    mask_lin_idx = as.integer(mask_idx) - 1L,  # Convert to 0-based
    dims = as.integer(dims),
    voxmm = as.numeric(voxmm),
    K = as.integer(K),
    compactness = as.numeric(compactness),
    max_iter = as.integer(max_iter),
    step_mm = as.numeric(step_mm),
    n_threads = as.integer(n_threads),
    seed_method = seed_method,
    enforce_connectivity = isTRUE(enforce_connectivity) || isTRUE(strict_connectivity),
    min_size = as.integer(min_size),
    connectivity = as.integer(connectivity),
    strict_connectivity = isTRUE(strict_connectivity),
    preserve_k = isTRUE(preserve_k),
    topup_iters = as.integer(topup_iters),
    grad_masked = as.numeric(grad_masked),
    seed_relocate_radius = as.integer(seed_relocate_radius),
    verbose = verbose
  )
  
  # Build ClusteredNeuroVol consistent with supervoxels.R
  kvol <- ClusteredNeuroVol(mask_arr, clusters = core$labels)
  
  # Prepare data for standardized result
  data_prep <- list(
    features = feat,
    coords = coords,
    mask_idx = mask_idx,
    n_voxels = length(mask_idx),
    dims = dims,
    spacing = voxmm
  )
  
  # Create standardized result
  result <- create_cluster4d_result(
    labels = core$labels,
    mask = mask,
    data_prep = data_prep,
    method = "slic",
    parameters = list(
      K = K,
      compactness = compactness,
      max_iter = max_iter,
      n_threads = n_threads,
      step_mm = step_mm,
      n_components = n_components,
      feature_norm = feature_norm,
      seed_method = seed_method,
      seed_relocate = seed_relocate,
      seed_relocate_radius = seed_relocate_radius,
      connectivity = connectivity,
      strict_connectivity = strict_connectivity,
      enforce_connectivity = enforce_connectivity,
      preserve_k = preserve_k,
      topup_iters = topup_iters,
      min_size = min_size
    ),
    metadata = list(
      centers = core$center_feats,
      coord_centers = core$center_coords
    ),
    compute_centers = FALSE  # Already computed by C++
  )
  
  # Ensure backward compatibility
  class(result) <- c("cluster_result", "list")
  result
}

# Helper function for finite difference gradient
.grad3d_fdiff <- function(arr) {
  dx <- arr * 0
  dy <- dx
  dz <- dx
  
  # Forward differences
  dx[-1,,] <- arr[-1,,] - arr[-nrow(arr),,]
  dy[,-1,] <- arr[,-1,] - arr[,-ncol(arr),]
  dz[,,-1] <- arr[,,-1] - arr[,,-dim(arr)[3]]
  
  # Gradient magnitude
  sqrt(dx*dx + dy*dy + dz*dz)
}

#' Visualize gradient used for seed relocation
#' 
#' Compute and return the gradient volume that would be used for seed relocation
#' in slic4d_supervoxels. Useful for debugging and visualization.
#' 
#' @param bvec A \code{NeuroVec} with dims (X, Y, Z, T).
#' @param mask A 3D \code{NeuroVol} (or logical array) indicating voxels to include.
#' @param method One of "correlation", "intensity", "spatial".
#' 
#' @return A 3D array containing the gradient values.
#' 
#' @export
slic4d_grad_summary <- function(bvec, mask, method = c("correlation", "intensity", "spatial")) {
  method <- match.arg(method)
  
  # Get dimensions and space
  dims <- if (inherits(mask, "NeuroVol")) dim(mask) else dim(mask)
  sp <- if (inherits(mask, "NeuroVol")) neuroim2::space(mask) else neuroim2::space(neuroim2::vol(mask))
  
  if (method == "correlation") {
    img4d <- as.array(bvec)
    # Convert mask to numeric array while preserving dimensions
    mask_numeric <- array(as.numeric(as.logical(mask)), dim = dims)
    grad3d <- correlation_gradient_cpp(img4d, mask_numeric)
    dim(grad3d) <- dims
    return(grad3d)
  } else if (method == "intensity") {
    img4d <- as.array(bvec)
    mean3d <- apply(img4d, c(1,2,3), mean)
    return(.grad3d_fdiff(mean3d))
  } else {
    if (!requireNamespace("neighborweights", quietly = TRUE)) {
      stop("neighborweights not installed for 'spatial' method.")
    }
    spatial_gradient(
      neuroim2::vol(apply(as.array(bvec), c(1,2,3), mean), space = sp),
      neuroim2::vol(as.logical(mask), space = sp)
    )
  }
}
