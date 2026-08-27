#' Fast 4D SLIC supervoxels (mask-aware, gradient relocation, preserve-K)
#'
#' Cluster a 4D \code{NeuroVec} (x,y,z,time) into compact, spatially contiguous
#' 3D supervoxels using an enhanced SLIC-style algorithm with mask-aware seeding,
#' gradient-based seed relocation, and topology-aware K preservation.
#'
#' Connectivity is the final partition constraint. With
#' \code{strict_connectivity = TRUE}, every returned label has exactly one
#' connected component under the requested neighborhood. If
#' \code{preserve_k = TRUE}, an adjacency-preserving repair reaches exactly K
#' whenever K is at least the number of connected mask components. When K is
#' below that topological minimum, strict connectivity takes precedence, the
#' minimum feasible number of labels is returned, and a warning is emitted.
#' Returned feature and coordinate centers are recomputed from the final public
#' labels in the original feature space and physical coordinate system.
#' 
#' @note Consider using \code{\link{cluster4d}} with \code{method = "slic"} for a
#' standardized interface across all clustering methods.
#'
#' @param bvec A \code{NeuroVec} with dims (X, Y, Z, T).
#'   Can also be a 3D \code{\link[neuroim2:NeuroVol-class]{NeuroVol}} for structural image segmentation,
#'   which will be automatically converted to a single-timepoint NeuroVec internally.
#' @param mask A 3D \code{NeuroVol} indicating voxels to include. Finite values
#'   greater than zero are included; zero and negative values are excluded.
#' @param K Finite integer target number of supervoxels, between one and the
#'   number of included voxels.
#' @param compactness Finite non-negative spatial vs feature tradeoff (like
#'   SLIC 'm'). Larger values produce more compact clusters.
#' @param max_iter Positive finite integer iteration limit (default 10).
#' @param n_threads Non-negative finite integer number of CPU threads. Zero
#'   requests RcppParallel's automatic thread count, one forces serial
#'   assignment, and values above one request that many threads. Assignment is
#'   intentionally serial for masks with fewer than 2,000 included voxels.
#' @param step_mm Optional approximate spacing between seeds in millimeters; if NULL,
#'   computed as cubic root of bounding-box volume / K.
#' @param n_components If > 0, random-project each voxel's time series to this dimension
#'   for speed (Johnson-Lindenstrauss style). 0 = use full time series.
#' @param feature_norm One of "zscale", "l2", "none".
#' @param seed_method One of "mask_poisson" (Poisson disk in mask), "mask_grid" (grid in mask),
#'   "grid" (regular grid), "farthest" (farthest point sampling).
#' @param seed_relocate One of "correlation" (correlation gradient), "intensity" (mean intensity gradient),
#'   "spatial" (spatial gradient using adjoin), "none" (no relocation).
#' @param seed_relocate_radius Search radius in voxels for gradient-based seed relocation (default 1).
#' @param connectivity Neighborhood connectivity: 6 (face neighbors) or 26 (all neighbors).
#' @param strict_connectivity Enforce exactly one connected component per label
#'   (default TRUE).
#' @param enforce_connectivity Deprecated alias for strict connectivity. When
#'   omitted, it follows \code{strict_connectivity}; if either supplied value is
#'   TRUE, connectivity is enforced.
#' @param preserve_k Return exactly K non-empty labels when compatible with
#'   strict connectivity (default FALSE). See Details.
#' @param topup_iters Non-negative number of native refinement iterations before
#'   the final connectivity and exact-K repair (default 2).
#' @param min_size Minimum component size (voxels) to keep before relabel (default 0).
#' @param verbose Logical.
#' @param .input_contract Internal prevalidated input receipt used by
#'   `cluster4d()`; direct callers should leave this as `NULL`.
#'
#' @return A \code{cluster4d_result} (also inheriting from
#' \code{cluster_result}) with elements:
#' \itemize{
#'   \item \code{clusvol}: \code{ClusteredNeuroVol} with the final labels
#'   \item \code{cluster}: integer vector of length = #masked voxels
#'   \item \code{centers}: matrix (actual_k x d_feat) center features
#'   \item \code{coord_centers}: matrix (actual_k x 3) spatial centers in mm
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
                              verbose = FALSE,
                              .input_contract = NULL) {
  
  # Allow single-volume NeuroVol by wrapping to NeuroVec
  bvec <- ensure_neurovec(bvec)

  method <- "slic4d_supervoxels"
  feature_norm <- match.arg(feature_norm)
  seed_method <- match.arg(seed_method)
  seed_relocate <- match.arg(seed_relocate)

  input <- .cluster4d_resolve_input_contract(
    .input_contract, bvec, mask, K, method
  )
  K <- input$n_clusters
  compactness <- .cluster4d_scalar_number(compactness, "compactness", method)
  max_iter <- .cluster4d_scalar_number(
    max_iter, "max_iter", method, integer = TRUE
  )
  n_threads <- .cluster4d_scalar_number(
    n_threads, "n_threads", method, integer = TRUE
  )
  n_components <- .cluster4d_scalar_number(
    n_components, "n_components", method, integer = TRUE
  )
  seed_relocate_radius <- .cluster4d_scalar_number(
    seed_relocate_radius, "seed_relocate_radius", method, integer = TRUE
  )
  topup_iters <- .cluster4d_scalar_number(
    topup_iters, "topup_iters", method, integer = TRUE
  )
  min_size <- .cluster4d_scalar_number(
    min_size, "min_size", method, integer = TRUE
  )
  strict_connectivity <- .cluster4d_scalar_logical(
    strict_connectivity, "strict_connectivity", method
  )
  enforce_connectivity <- if (missing(enforce_connectivity)) {
    strict_connectivity
  } else {
    .cluster4d_scalar_logical(
      enforce_connectivity, "enforce_connectivity", method
    )
  }
  preserve_k <- .cluster4d_scalar_logical(preserve_k, "preserve_k", method)
  verbose <- .cluster4d_scalar_logical(verbose, "verbose", method)
  connectivity <- if (missing(connectivity)) {
    26L
  } else {
    .cluster4d_scalar_number(
      connectivity, "connectivity", method, integer = TRUE
    )
  }
  if (!connectivity %in% c(6L, 26L)) {
    stop(method, ": connectivity must be 6 or 26", call. = FALSE)
  }
  if (compactness < 0) stop(method, ": compactness must be non-negative", call. = FALSE)
  if (max_iter < 1L) stop(method, ": max_iter must be positive", call. = FALSE)
  if (n_threads < 0L) stop(method, ": n_threads must be non-negative", call. = FALSE)
  if (n_components < 0L) stop(method, ": n_components must be non-negative", call. = FALSE)
  if (seed_relocate_radius < 0L) {
    stop(method, ": seed_relocate_radius must be non-negative", call. = FALSE)
  }
  if (topup_iters < 0L) stop(method, ": topup_iters must be non-negative", call. = FALSE)
  if (min_size < 0L) stop(method, ": min_size must be non-negative", call. = FALSE)
  if (!is.null(step_mm)) {
    step_mm <- .cluster4d_scalar_number(step_mm, "step_mm", method)
    if (step_mm <= 0) stop(method, ": step_mm must be positive", call. = FALSE)
  }

  n_timepoints <- nrow(input$features)
  if (n_timepoints <= 1 && seed_relocate == "correlation") {
    if (verbose) {
      message("slic4d_supervoxels: correlation relocation not applicable for single timepoint; using intensity gradient instead")
    }
    seed_relocate <- "intensity"
  }
  
  # Reuse the validated extraction and materialize its N x T view once.
  data_prep <- prepare_cluster4d_data(
    bvec, mask,
    scale_features = FALSE,
    scale_coords = FALSE,
    input_contract = input
  )
  mask_arr <- input$mask
  sp <- neuroim2::space(mask)

  mask_idx <- data_prep$mask_idx
  
  # Get dimensions and voxel sizes
  dims <- dim(mask_arr)
  voxmm <- neuroim2::spacing(sp)
  
  # Spatial coordinates in mm for masked voxels
  coords <- data_prep$coords
  feat <- data_prep$features

  # Thread selection is passed to each native parallelFor invocation. Never
  # mutate RcppParallel's process-wide environment from a clustering call.
  n_threads_requested <- as.integer(n_threads)
  parallel_requested <- n_threads_requested != 1L
  n_threads_eff <- n_threads_requested
  if (n_threads_eff == 0L && length(mask_idx) < 2000L) {
    n_threads_eff <- 1L
  }
  
  # Optional feature normalization
  if (feature_norm == "zscale") {
    # Z-score each timepoint across voxels
    mu <- colMeans(feat)
    sdv <- apply(feat, 2, sd)
    sdv[!is.finite(sdv) | sdv == 0] <- 1
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
      # Use spatial gradient via adjoin
      if (!requireNamespace("adjoin", quietly = TRUE)) {
        warning("adjoin not installed; falling back to intensity gradient")
        img4d <- as.array(bvec)
        mean3d <- apply(img4d, c(1,2,3), mean)
        grad3d <- .grad3d_fdiff(mean3d)
      } else {
        grad3d <- spatial_gradient(
          neuroim2::NeuroVol(apply(as.array(bvec), c(1,2,3), mean), space = sp),
          neuroim2::NeuroVol(mask_arr, space = sp)
        )
        if (inherits(grad3d, "NeuroVol")) {
          grad3d <- as.array(grad3d)
        }
      }
      grad_masked <- grad3d[mask_idx]
    }
  }
  if (length(grad_masked) && any(!is.finite(grad_masked))) {
    stop(method, ": relocation gradient must be finite", call. = FALSE)
  }

  connectivity_required <- strict_connectivity || enforce_connectivity
  
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
    n_threads = as.integer(n_threads_eff),
    seed_method = seed_method,
    enforce_connectivity = connectivity_required,
    min_size = as.integer(min_size),
    connectivity = as.integer(connectivity),
    strict_connectivity = isTRUE(strict_connectivity),
    preserve_k = isTRUE(preserve_k),
    topup_iters = as.integer(topup_iters),
    grad_masked = as.numeric(grad_masked),
    seed_relocate_radius = as.integer(seed_relocate_radius),
    verbose = verbose
  )

  preserve_k_feasible <- TRUE
  final_labels <- as.integer(core$labels)
  if (preserve_k) {
    if (connectivity_required) {
      graph_info <- .exact_k_graph(mask, connectivity)
      minimum_k <- length(unique(as.integer(graph_info$components)))
      repair_target <- max(K, minimum_k)
      preserve_k_feasible <- K >= minimum_k
      final_labels <- force_exact_k(
        final_labels, feat, repair_target, mask, connectivity,
        graph_info = graph_info, min_cluster_size = 1L
      )
      if (!preserve_k_feasible) {
        warning(
          method, ": strict connectivity takes precedence because K = ", K,
          " is below the mask's ", minimum_k, " connected components",
          call. = FALSE
        )
      }
    } else {
      final_labels <- .slic_force_exact_k_unconstrained(
        final_labels, feat, coords, K
      )
    }
  }
  
  # Finalize once in the original feature space. Native centers remain typed
  # metadata because they summarize the normalized/projected working space.
  result <- structure(
    list(
      labels = final_labels,
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
        native_centers = core$center_feats,
        native_coord_centers = core$center_coords,
        native_actual_k = core$actual_k,
        assignment_parallel_requested = parallel_requested,
        assignment_parallel_used = isTRUE(core$assignment_parallel_used),
        native_n_threads = as.integer(n_threads_eff),
        native_connectivity_changed = core$connectivity_changed,
        native_connectivity_elapsed_ms = core$connectivity_elapsed_ms,
        native_center_recompute_elapsed_ms = core$center_recompute_elapsed_ms,
        preserve_k_feasible = preserve_k_feasible
      )
    ),
    class = c("cluster4d_result", "cluster_result", "list")
  )

  finalize_cluster4d_result(
    result, bvec, mask, "slic", result$parameters,
    input_contract = input,
    data = data_prep
  )
}

.slic_force_exact_k_unconstrained <- function(labels, features, coords, K) {
  labels <- .exact_k_relabel(as.integer(labels))
  features <- as.matrix(features)
  coords <- as.matrix(coords)

  while (max(labels) > K) {
    summaries <- .exact_k_cluster_summaries(labels, cbind(features, coords))
    pairs <- utils::combn(seq_len(max(labels)), 2L)
    difference <- summaries$centers[pairs[1L, ], , drop = FALSE] -
      summaries$centers[pairs[2L, ], , drop = FALSE]
    costs <- summaries$counts[pairs[1L, ]] * summaries$counts[pairs[2L, ]] /
      (summaries$counts[pairs[1L, ]] + summaries$counts[pairs[2L, ]]) *
      rowSums(difference^2)
    pick <- order(costs, pairs[1L, ], pairs[2L, ])[[1L]]
    labels[labels == pairs[2L, pick]] <- pairs[1L, pick]
    labels <- .exact_k_relabel(labels)
  }

  while (max(labels) < K) {
    best <- NULL
    for (cluster in seq_len(max(labels))) {
      members <- which(labels == cluster)
      if (length(members) <= 1L) next
      values <- cbind(features[members, , drop = FALSE], coords[members, , drop = FALSE])
      center <- colMeans(values)
      gains <- rowSums((values - matrix(
        center, nrow(values), ncol(values), byrow = TRUE
      ))^2)
      pick <- order(-gains, members)[[1L]]
      candidate <- list(cluster = cluster, voxel = members[[pick]], gain = gains[[pick]])
      if (is.null(best) || candidate$gain > best$gain ||
          (candidate$gain == best$gain && candidate$cluster < best$cluster) ||
          (candidate$gain == best$gain && candidate$cluster == best$cluster &&
             candidate$voxel < best$voxel)) {
        best <- candidate
      }
    }
    if (is.null(best)) stop("slic4d_supervoxels: cannot create K non-empty labels", call. = FALSE)
    labels[[best$voxel]] <- max(labels) + 1L
  }
  .exact_k_relabel(labels)
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
#' @param mask A 3D \code{NeuroVol} indicating voxels to include.
#' @param method One of "correlation", "intensity", "spatial".
#' 
#' @return A 3D array containing the gradient values.
#' 
#' @export
slic4d_grad_summary <- function(bvec, mask, method = c("correlation", "intensity", "spatial")) {
  method <- match.arg(method)
  
  # Get dimensions and space
  mask_arr <- cluster4d_mask_array(mask, "slic4d_grad_summary")
  dims <- dim(mask_arr)
  sp <- neuroim2::space(mask)
  
  if (method == "correlation") {
    img4d <- as.array(bvec)
    # Convert mask to numeric array while preserving dimensions
    mask_numeric <- array(as.numeric(mask_arr), dim = dims)
    grad3d <- correlation_gradient_cpp(img4d, mask_numeric)
    dim(grad3d) <- dims
    return(grad3d)
  } else if (method == "intensity") {
    img4d <- as.array(bvec)
    mean3d <- apply(img4d, c(1,2,3), mean)
    return(.grad3d_fdiff(mean3d))
  } else {
    if (!requireNamespace("adjoin", quietly = TRUE)) {
      stop("adjoin not installed for 'spatial' method.")
    }
    spatial_gradient(
      neuroim2::NeuroVol(apply(as.array(bvec), c(1,2,3), mean), space = sp),
      neuroim2::NeuroVol(mask_arr, space = sp)
    )
  }
}
