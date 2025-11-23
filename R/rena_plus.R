#' ReNA++: Edge-aware Reciprocal Multi-Level ReNA with Ward refinement
#'
#' Two-stage pipeline:
#' 1) Edge-aware reciprocal-nearest-neighbor ReNA coarsening to K' = r * K super-voxels.
#' 2) Spatially-constrained Ward refinement on the super-graph to reach exactly K clusters.
#'
#' @param bvec A \code{NeuroVec} providing voxel time-series.
#' @param mask A \code{NeuroVol} mask; non-zero voxels are included.
#' @param K Target number of clusters.
#' @param r Over-clustering factor; coarsening stops at K' = ceiling(r * K). Default 5.
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
                      r = 5,
                      lambda = 1,
                      grad_img = NULL,
                      connectivity = 26,
                      max_iterations = 50,
                      verbose = FALSE,
                      ...) {

  validate_cluster4d_inputs(bvec, mask, K, "rena_plus")

  if (!connectivity %in% c(6, 18, 26, 27)) {
    stop("rena_plus: connectivity must be one of 6, 18, 26, or 27")
  }

  vox_idx <- which(mask > 0)
  n_vox <- length(vox_idx)

  # Features: series returns T x N; scale then transpose to N x T
  feature_mat <- neuroim2::series(bvec, vox_idx)
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

  # Build grid adjacency (symmetric, binary)
  adj <- build_grid_adjacency(mask, vox_idx, ifelse(connectivity == 27, 26, connectivity))

  K_prime <- max(1L, as.integer(ceiling(K * r)))

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
      final_n_clusters = final_k
    ),
    compute_centers = TRUE,
    center_method = "mean"
  )

  class(result) <- c("rena_plus_cluster_result", class(result))
  result
}

#' Cluster4d wrapper for ReNA++
#'
#' @inheritParams cluster4d
#' @param r Over-clustering factor for coarsening stage.
#' @param lambda Gradient penalty weight.
#' @param grad_img Optional gradient/intensity vector (see \code{rena_plus}).
#' @return A \code{cluster4d_result}.
#' @export
cluster4d_rena_plus <- function(vec, mask, n_clusters = 100,
                                spatial_weight = 0.5,
                                r = 5,
                                lambda = 1,
                                grad_img = NULL,
                                connectivity = 26,
                                max_iterations = 50,
                                verbose = FALSE,
                                ...) {

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
      spatial_weight = spatial_weight,
      connectivity = connectivity,
      r = r,
      lambda = lambda
    )
  )

  result
}

#' Internal: build 3D grid adjacency for masked voxels
#' @keywords internal
build_grid_adjacency <- function(mask, mask_idx, connectivity) {
  dims <- dim(mask)
  n_vox <- length(mask_idx)

  idx_map <- array(-1L, dim = dims)
  idx_map[mask_idx] <- seq_len(n_vox)

  offsets <- switch(as.character(connectivity),
                    "6" = rbind(
                      c( 1, 0, 0), c(-1, 0, 0),
                      c( 0, 1, 0), c( 0,-1, 0),
                      c( 0, 0, 1), c( 0, 0,-1)
                    ),
                    "18" = {
                      base <- expand.grid(-1:1, -1:1, -1:1)
                      base <- as.matrix(base)
                      base <- base[rowSums(abs(base)) <= 2 & rowSums(abs(base)) > 0, ]
                      base
                    },
                    "26" = {
                      base <- expand.grid(-1:1, -1:1, -1:1)
                      base <- as.matrix(base)
                      base <- base[rowSums(abs(base)) > 0, ]
                      base
                    },
                    stop("Unsupported connectivity"))

  edges_i <- integer()
  edges_j <- integer()

  # iterate masked voxels
  coords <- index_to_grid(mask, mask_idx)
  for (v in seq_len(n_vox)) {
    xyz <- coords[v, ]
    for (k in seq_len(nrow(offsets))) {
      nb <- xyz + offsets[k, ]
      if (any(nb < 1) || nb[1] > dims[1] || nb[2] > dims[2] || nb[3] > dims[3]) {
        next
      }
      nb_idx <- idx_map[nb[1], nb[2], nb[3]]
      if (nb_idx > 0 && nb_idx > v) {  # upper-tri to dedup
        edges_i <- c(edges_i, v)
        edges_j <- c(edges_j, nb_idx)
      }
    }
  }

  if (length(edges_i) == 0) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(n_vox, n_vox)))
  }

  Matrix::sparseMatrix(
    i = c(edges_i, edges_j),
    j = c(edges_j, edges_i),
    x = 1,
    dims = c(n_vox, n_vox)
  )
}

# Backwards-compatible alias
#' @export
er_ml_rena <- function(bvec, mask, ...) {
  rena_plus(bvec, mask, ...)
}
