#' ReNA++: Edge-aware Reciprocal Multi-Level ReNA with Ward refinement
#'
#' Two-stage pipeline:
#' 1) Edge-aware reciprocal-nearest-neighbor ReNA coarsening toward
#'    K' = ceiling(r * K) super-voxels.
#' 2) Spatially-constrained Ward refinement on the super-graph to reach exactly K clusters.
#'
#' A reciprocal-neighbor level is retained only when applying the complete
#' level would leave at least K' super-voxels. This prevents the coarsening
#' stage from overshooting below K and leaves the remaining one-at-a-time
#' adjacent merges to Ward refinement. Ward costs are recomputed from the
#' current size-weighted means before every accepted merge; equal-cost edges
#' are resolved by ascending super-voxel identifiers. Because the adjacency
#' constraint can expose a new edge after a merge, individual Ward deltas need
#' not be ordered, although every delta is non-negative and the cumulative
#' objective is non-decreasing. When Ward refinement is needed, its complete
#' trace and queue diagnostics are available in \code{result$metadata$ward}.
#'
#' @param bvec A \code{NeuroVec} providing voxel time-series.
#' @param mask A \code{NeuroVol} mask; finite positive voxels are included.
#' @param K Target number of clusters.
#' @param r Finite over-clustering factor at least 1. Coarsening approaches
#'   K' = ceiling(r * K) without applying a level that would overshoot below
#'   K'. Default 2.
#' @param lambda Gradient penalty strength (>=0). 0 disables edge weighting.
#' @param grad_img Optional numeric vector: gradient/intensity per voxel. Either length = prod(dim(mask))
#'   (in which case values are subset by mask) or length = number of masked voxels. If NULL, no gradient.
#' @param connectivity Neighborhood connectivity (6, 18, or 26). Default 26.
#' @param max_iterations Max iterations for coarsening stage.
#' @param verbose Logical for progress messages.
#' @param ... Reserved for future options.
#'
#' @return A \code{cluster4d_result} with additional class \code{rena_plus_cluster_result}.
#' @export
rena_plus <- function(bvec, mask,
                      K = 100,
                      r = 2,
                      lambda = 1,
                      grad_img = NULL,
                      connectivity = 26,
                      max_iterations = 50,
                      verbose = FALSE,
                      ...) {

  bvec <- ensure_neurovec(bvec)
  input <- validate_cluster4d_inputs(bvec, mask, K, "rena_plus")
  K <- input$n_clusters
  r <- .cluster4d_scalar_number(r, "r", "rena_plus")
  lambda <- .cluster4d_scalar_number(lambda, "lambda", "rena_plus")
  connectivity <- .cluster4d_scalar_number(
    connectivity, "connectivity", "rena_plus", integer = TRUE
  )
  max_iterations <- .cluster4d_scalar_number(
    max_iterations, "max_iterations", "rena_plus", integer = TRUE
  )
  verbose <- .cluster4d_scalar_logical(verbose, "verbose", "rena_plus")
  if (r < 1) stop("rena_plus: r must be at least 1", call. = FALSE)
  if (lambda < 0) stop("rena_plus: lambda must be non-negative", call. = FALSE)
  if (max_iterations < 1L) {
    stop("rena_plus: max_iterations must be positive", call. = FALSE)
  }

  if (!connectivity %in% c(6L, 18L, 26L)) {
    stop("rena_plus: connectivity must be one of 6, 18, or 26", call. = FALSE)
  }

  vox_idx <- input$mask_idx
  n_vox <- length(vox_idx)

  # Features: series returns T x N; scale then transpose to N x T
  feature_mat <- neuroim2::series(bvec, vox_idx)
  if (!is.matrix(feature_mat)) {
    feature_mat <- matrix(as.numeric(feature_mat), nrow = 1L, ncol = n_vox)
  }
  feature_mat <- base::scale(as.matrix(feature_mat), center = TRUE, scale = TRUE)
  feature_mat[is.na(feature_mat)] <- 0
  feature_mat <- t(feature_mat)  # voxels x time

  # Gradient vector
  total_vox <- prod(dim(mask))
  grad_vec <- if (is.null(grad_img)) {
    numeric(0)
  } else {
    if (length(grad_img) == total_vox) {
      as.numeric(grad_img[vox_idx])
    } else if (length(grad_img) == n_vox) {
      as.numeric(grad_img)
    } else {
      stop("rena_plus: grad_img must be NULL, length mask, or length masked voxels")
    }
  }
  if (length(grad_vec) && any(!is.finite(grad_vec))) {
    stop("rena_plus: grad_img values must be finite", call. = FALSE)
  }

  # Build grid adjacency (symmetric, binary)
  adj <- build_grid_adjacency(mask, vox_idx, connectivity)

  K_prime <- as.integer(min(n_vox, max(K, ceiling(K * r))))

  if (verbose) {
    message("ReNA++: coarsening to K' = ", K_prime, " (r = ", r, "), lambda = ", lambda)
  }

  coarse <- rena_rnn_coarse_cpp(
    X = feature_mat,
    G = adj,
    grad_img = grad_vec,
    stop_at = K_prime,
    lambda = lambda,
    max_iter = max_iterations
  )

  k_coarse <- nrow(coarse$X_coarse)
  k_target <- min(K, k_coarse)

  if (k_target < 1) {
    stop("rena_plus: coarsening produced zero clusters")
  }

  if (k_coarse > k_target) {
    ward <- ward_on_supervoxels_cpp(
      X_coarse = coarse$X_coarse,
      G_coarse = coarse$G_coarse,
      sizes = coarse$sizes_coarse,
      n_clusters = k_target
    )
    labels_super <- ward$labels_super
    final_k <- ward$n_clusters
  } else {
    # Already at or below target; keep current labels
    labels_super <- seq_len(k_coarse)
    final_k <- k_coarse
  }

  if (final_k != K) {
    stop(sprintf(
      "rena_plus: produced %d clusters but exactly %d were requested",
      final_k, K
    ), call. = FALSE)
  }

  final_labels <- labels_super[coarse$labels_coarse]

  data_prep <- list(
    features = feature_mat,  # voxels x time
    coords = index_to_coord(mask, vox_idx),
    mask_idx = vox_idx,
    n_voxels = n_vox,
    dims = dim(mask),
    spacing = spacing(mask)
  )

  result <- create_cluster4d_result(
    labels = final_labels,
    mask = mask,
    data_prep = data_prep,
    method = "rena_plus",
    parameters = list(
      K = K,
      requested_K = K,
      K_prime = K_prime,
      r = r,
      lambda = lambda,
      connectivity = connectivity,
      max_iterations = max_iterations
    ),
    metadata = list(
      coarse_n_clusters = k_coarse,
      coarse_iterations = coarse$iterations,
      coarse_stopped_before_overshoot = coarse$stopped_before_overshoot,
      final_n_clusters = final_k,
      ward = if (exists("ward", inherits = FALSE)) ward else NULL
    ),
    compute_centers = TRUE,
    center_method = "mean"
  )

  result <- finalize_cluster4d_result(
    result, bvec, mask, "rena_plus", result$parameters
  )
  class(result) <- c("rena_plus_cluster_result", class(result))
  result
}

#' Cluster4d wrapper for ReNA++
#'
#' @inheritParams cluster4d
#' @param n_clusters Exact finite integer target number of clusters, from one
#'   through the number of included mask voxels.
#' @param r Over-clustering factor for coarsening stage.
#' @param lambda Gradient penalty weight.
#' @param grad_img Optional gradient/intensity vector (see \code{rena_plus}).
#' @return A \code{cluster4d_result}.
#' @export
cluster4d_rena_plus <- function(vec, mask, n_clusters = 100,
                                spatial_weight = 0.5,
                                r = 2,
                                lambda = 1,
                                grad_img = NULL,
                                connectivity = 26,
                                max_iterations = 50,
                                verbose = FALSE,
                                ...) {

  vec <- ensure_neurovec(vec)
  if (!missing(spatial_weight)) {
    stop("cluster4d_rena_plus: spatial_weight is not supported", call. = FALSE)
  }
  validate_cluster4d_inputs(vec, mask, n_clusters, "cluster4d_rena_plus")

  result <- rena_plus(
    bvec = vec,
    mask = mask,
    K = n_clusters,
    r = r,
    lambda = lambda,
    grad_img = grad_img,
    connectivity = connectivity,
    max_iterations = max_iterations,
    verbose = verbose,
    ...
  )

  if (!"cluster4d_result" %in% class(result)) {
    class(result) <- c("cluster4d_result", class(result))
  }

  result$method <- "rena_plus"

  result$parameters <- modifyList(
    result$parameters,
    list(
      n_clusters_requested = n_clusters,
      connectivity = connectivity,
      r = r,
      lambda = lambda
    )
  )

  finalize_cluster4d_result(result, vec, mask, "rena_plus", result$parameters)
}

#' Internal: build 3D grid adjacency for masked voxels (fast C++ version)
#' @keywords internal
build_grid_adjacency <- function(mask, mask_idx, connectivity) {
  dims <- dim(mask)
  if (length(dims) != 3L) {
    stop("build_grid_adjacency: mask must be 3D")
  }
  build_grid_adjacency_cpp(
    mask_idx = as.integer(mask_idx),
    dims     = as.integer(dims),
    connectivity = as.integer(connectivity)
  )
}

#' Backwards-compatible alias for rena_plus
#'
#' @inheritParams rena_plus
#' @return Same as \code{\link{rena_plus}}.
#' @export
er_ml_rena <- function(bvec, mask, ...) {
  rena_plus(bvec, mask, ...)
}
