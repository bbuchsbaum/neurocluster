#' Gradient-guided geodesic clustering on a masked grid
#'
#' Compresses voxel features, selects low-gradient seeds, and assigns voxels by
#' deterministic multi-source Dijkstra propagation. Edges are exactly the
#' requested masked-grid neighbors; excluded gaps are never bridged. The edge
#' objective is `alpha * feature_distance + (1 - alpha) * physical_distance /
#' compactness`.
#'
#' @param vec A `NeuroVec` or `SparseNeuroVec` containing the input features.
#' @param mask A `NeuroVol`; exactly finite values greater than zero are included.
#' @param K Positive requested cluster count.
#' @param n_components Requested SVD feature count.
#' @param variance_threshold Minimum retained variance in `[0, 1]`.
#' @param alpha Feature-distance weight in `[0, 1]`; zero is spatial-only and
#'   one is feature-only.
#' @param compactness Positive finite physical length scale. Smaller values
#'   impose a stronger spatial penalty. The default uses masked voxel measure,
#'   K, and the effective dimension, so single-slice masks remain well-defined.
#' @param max_refinement_iter Non-negative number of boundary passes. Refinement
#'   uses the same feature/spatial mixture and final labels are repaired through
#'   the shared adjacency-preserving exact-K engine.
#' @param verbose Whether to report phase progress.
#' @param use_irlba Whether large SVDs may use `irlba`.
#' @param use_rsvd Whether randomized SVD may be used when available.
#' @param connectivity Exact grid connectivity: 6, 18, or 26.
#'
#' @return A typed `g3s_result` and `cluster4d_result` with K by T centers,
#'   physical coordinate centers, graph provenance, and physical-scale metadata.
#'
#' @references Dijkstra, E. W. (1959). A note on two problems in connexion
#'   with graphs. Numerische Mathematik, 1, 269-271.
#' @importFrom neuroim2 series spacing
#' @export
cluster4d_g3s <- function(vec, mask, K = 100,
                         n_components = 15,
                         variance_threshold = 0.95,
                         alpha = 0.5,
                         compactness = NULL,
                         max_refinement_iter = 3,
                         verbose = FALSE,
                         use_irlba = TRUE,
                         use_rsvd = TRUE,
                         connectivity = 26) {

  # Use the package-wide data, geometry, mask, and K contract.
  input_contract <- validate_cluster4d_inputs(
    vec, mask, K, "cluster4d_g3s"
  )
  K <- input_contract$n_clusters
  n_components <- .cluster4d_scalar_number(
    n_components, "n_components", "cluster4d_g3s", lower = 1, integer = TRUE
  )
  variance_threshold <- .cluster4d_scalar_number(
    variance_threshold, "variance_threshold", "cluster4d_g3s",
    lower = 0, upper = 1
  )
  verbose <- .cluster4d_scalar_logical(verbose, "verbose", "cluster4d_g3s")
  use_irlba <- .cluster4d_scalar_logical(
    use_irlba, "use_irlba", "cluster4d_g3s"
  )
  use_rsvd <- .cluster4d_scalar_logical(
    use_rsvd, "use_rsvd", "cluster4d_g3s"
  )

  mask.idx <- input_contract$mask_idx
  n_voxels <- length(mask.idx)

  alpha <- .cluster4d_scalar_number(
    alpha, "alpha", "cluster4d_g3s", lower = 0, upper = 1
  )
  max_refinement_iter <- .cluster4d_scalar_number(
    max_refinement_iter, "max_refinement_iter", "cluster4d_g3s",
    lower = 0, integer = TRUE
  )
  connectivity <- .cluster4d_scalar_number(
    connectivity, "connectivity", "cluster4d_g3s", integer = TRUE
  )
  if (!connectivity %in% c(6L, 18L, 26L)) {
    stop("cluster4d_g3s: connectivity must be 6, 18, or 26", call. = FALSE)
  }

  # The contract already pulled and validated the T x N series, including the
  # single-frame case that neuroim2 simplifies to a bare vector.
  feature_mat_raw <- input_contract$features   # T x N matrix
  features_by_voxel <- t(feature_mat_raw)      # N x T, materialised once

  n_timepoints <- nrow(feature_mat_raw)
  use_cosine <- n_timepoints > 1

  if (verbose) {
    message("G3S: Starting with ", n_voxels, " voxels, ",
            n_timepoints, " timepoints, target K=", K)
  }

  # =============================================================================
  # Phase 1: Hyper-Compression
  # =============================================================================

  if (use_cosine) {
    if (verbose) {
      message("Phase 1: SVD compression (", n_timepoints, " -> ", n_components, " dims)")
    }

    compressed <- compress_features_svd(
      feature_mat = features_by_voxel,  # Needs N x T
      n_components = n_components,
      variance_threshold = variance_threshold,
      use_irlba = use_irlba,
      use_rsvd = use_rsvd
    )

    feature_mat_compressed <- compressed$features  # N x M matrix (normalized)
    actual_components <- compressed$n_components
    variance_explained <- compressed$variance_explained

    if (verbose) {
      message("  Compressed to ", actual_components, " components (",
              round(variance_explained * 100, 1), "% variance)")
    }
  } else {
    if (verbose) {
      message("Phase 1: single timepoint detected; using Euclidean features without SVD")
    }
    feature_mat_compressed <- base::scale(features_by_voxel, center = TRUE, scale = TRUE)
    feature_mat_compressed[is.na(feature_mat_compressed)] <- 0
    actual_components <- 1
    variance_explained <- 1
  }

  # =============================================================================
  # Phase 2: Gradient-Based Seeding
  # =============================================================================

  if (verbose) {
    message("Phase 2: Finding gradient-based seeds")
  }

  # Build the exact masked-grid graph once and reuse it for seeding,
  # propagation, refinement, exact-K repair, and topology validation. `neib`
  # carries the .exact_k_graph() fields, so force_exact_k() below never has to
  # rebuild the same adjacency.
  neib <- build_grid_neighbors_g3s(mask, mask.idx, connectivity)
  coords <- neib$coords
  graph_components <- neib$components
  minimum_clusters <- max(graph_components)
  if (K < minimum_clusters) {
    .exact_k_abort(
      "disconnected_mask", K, minimum_clusters, n_voxels,
      paste0(
        "the requested ", connectivity,
        "-neighbor mask graph has ", minimum_clusters, " components"
      )
    )
  }
  k_neighbors <- ncol(neib$nn.index)

  # The masked physical scale drives both the seed inhibition radius and the
  # default compactness, so compute it once and use it for both. Deriving the
  # radius from the coordinate bounding box instead badly overestimates it on
  # non-convex masks.
  scale_info <- .g3s_spatial_scale(mask, mask.idx, K)

  seed_indices <- find_gradient_seeds_g3s(
    feature_mat = feature_mat_compressed,
    coords = coords,
    K = K,
    k_neighbors = k_neighbors,
    distance = if (use_cosine) "cosine" else "euclidean",
    knn = neib,
    spatial_scale = scale_info$scale
  )
  seed_indices <- .g3s_cover_components(
    seed_indices, graph_components, feature_mat_compressed,
    neib$nn.index, use_cosine
  )

  actual_K <- length(seed_indices)

  if (actual_K < K && verbose) {
    message("  Found ", actual_K, " spatially separated seeds (requested ", K, ")")
  }

  # =============================================================================
  # Phase 3: Geodesic Propagation
  # =============================================================================

  if (verbose) {
    message("Phase 3: Geodesic cluster propagation")
  }

  # Auto-compute compactness if not provided
  if (is.null(compactness)) {
    compactness <- scale_info$scale
    if (verbose) {
      message("  Auto-computed compactness: ", round(compactness, 2))
    }
  }
  compactness <- .cluster4d_scalar_number(
    compactness, "compactness", "cluster4d_g3s",
    lower = .Machine$double.eps
  )

  compressed_by_component <- t(feature_mat_compressed)  # M x N, materialised once

  # Call optimized C++ propagation
  labels <- g3s_propagate_cpp(
    feature_mat = compressed_by_component,
    seed_indices = as.integer(seed_indices),
    neighbor_indices = neib$nn.index,
    neighbor_dists = neib$nn.dist,
    alpha = alpha,
    compactness = compactness
  )

  if (any(labels <= 0L)) {
    stop("cluster4d_g3s: propagation left included voxels unlabeled", call. = FALSE)
  }

  # =============================================================================
  # Phase 4: Boundary Refinement (optional)
  # =============================================================================

  if (max_refinement_iter > 0) {
    if (verbose) {
      message("Phase 4: Boundary refinement (", max_refinement_iter, " iterations)")
    }

    # Refinement keeps every label to a single connected component, so the
    # exact-K engine below has almost nothing left to merge.
    labels <- refine_boundaries_g3s_cpp(
      labels = as.integer(labels),
      feature_mat = compressed_by_component,
      coords = coords,
      neighbor_indices = neib$nn.index,
      alpha = alpha,
      compactness = compactness,
      max_iter = as.integer(max_refinement_iter),
      enforce_connectivity = TRUE
    )
    if (any(labels <= 0L)) {
      stop("cluster4d_g3s: refinement left included voxels unlabeled",
           call. = FALSE)
    }
  }

  labels <- .exact_k_connected_labels(labels, neib$graph, neib$edges)
  labels <- force_exact_k(
    labels, features_by_voxel, K,
    mask = mask, connectivity = connectivity, graph_info = neib
  )

  # =============================================================================
  # Create Result Object
  # =============================================================================

  # Prepare data for standardized result
  # Use original features for center computation (same dimensionality as input)
  data_prep <- list(
    features = features_by_voxel,  # N x T for compute_cluster_centers
    coords = coords,
    mask_idx = mask.idx,
    n_voxels = n_voxels,
    dims = dim(mask),
    spacing = spacing(mask),
    geometry = input_contract$geometry
  )

  # Create standardized result
  result <- create_cluster4d_result(
    labels = labels,
    mask = mask,
    data_prep = data_prep,
    method = "g3s",
    parameters = list(
      K_requested = K,
      # force_exact_k() guarantees exactly K connected labels; the seed count
      # is reported separately because it is not the final cluster count.
      K_actual = length(unique(labels)),
      n_seeds = actual_K,
      n_components = actual_components,
      variance_threshold = variance_threshold,
      variance_explained = variance_explained,
      alpha = alpha,
      compactness = compactness,
      max_refinement_iter = max_refinement_iter,
      feature_metric = if (use_cosine) "cosine" else "euclidean",
      connectivity = connectivity
    ),
    metadata = list(
      # Masked-voxel indices of the propagation seeds. Labels are re-derived by
      # the exact-K repair afterwards, so seed i does not index cluster i.
      seed_indices = seed_indices,
      graph = list(
        contract = "exact_masked_grid",
        connectivity = connectivity,
        n_edges = nrow(neib$edges),
        n_components = minimum_clusters
      ),
      spatial = list(
        units = "physical",
        effective_dimension = scale_info$dimension,
        auto_scale = scale_info$scale
      ),
      compression_ratio = if (use_cosine) n_timepoints / actual_components else 1,
      svd_rotation = if (use_cosine) compressed$rotation else NULL,
      svd_singular_values = if (use_cosine) compressed$singular_values else NULL,
      svd_center = if (use_cosine) compressed$center else NULL,
      svd_scale = if (use_cosine) compressed$scale else NULL
    ),
    compute_centers = TRUE,
    center_method = "mean"
  )

  # Add G3S-specific class
  class(result) <- c("g3s_result", "cluster4d_result", "cluster_result", "list")

  if (verbose) {
    message("G3S complete: ", result$n_clusters, " clusters formed")
  }

  finalize_cluster4d_result(
    result, vec, mask, "g3s", result$parameters, data = data_prep
  )
}

#' Build exact masked-grid neighbors for G3S.
#'
#' @keywords internal
#' @noRd
build_grid_neighbors_g3s <- function(mask, mask_idx, connectivity,
                                     graph_info = NULL) {
  graph_info <- .exact_k_resolve_graph(graph_info, mask, connectivity)
  if (!identical(as.integer(mask_idx), as.integer(graph_info$mask_idx))) {
    stop("build_grid_neighbors_g3s: mask index contract mismatch", call. = FALSE)
  }
  n <- length(mask_idx)
  ptr <- graph_info$neighbor_ptr
  idx <- graph_info$neighbor_idx
  degree <- diff(ptr)
  max_degree <- if (n) max(degree) else 0L
  neighbor_indices <- matrix(0L, nrow = n, ncol = max_degree)
  neighbor_dists <- matrix(Inf, nrow = n, ncol = max_degree)
  coords <- .cluster4d_index_to_coord(mask, mask_idx)
  if (max_degree > 0L) {
    # Fill both padded matrices from the CSR adjacency in one vectorised pass;
    # the per-voxel R loop this replaces dominated graph construction.
    from <- rep.int(seq_len(n), degree)
    slot <- sequence(degree)
    position <- from + (slot - 1L) * n
    neighbor_indices[position] <- idx
    delta <- coords[idx, , drop = FALSE] - coords[from, , drop = FALSE]
    neighbor_dists[position] <- sqrt(.rowSums(delta * delta, length(idx), 3L))
  }
  c(
    graph_info,
    list(nn.index = neighbor_indices, nn.dist = neighbor_dists, coords = coords)
  )
}

.g3s_spatial_scale <- function(mask, mask_idx, K) {
  grid <- neuroim2::index_to_grid(mask, mask_idx)
  grid <- matrix(as.integer(grid), ncol = 3L)
  active <- vapply(
    seq_len(3L), function(axis) length(unique(grid[, axis])) > 1L,
    logical(1L)
  )
  dimension <- max(1L, sum(active))
  voxel_spacing <- as.numeric(neuroim2::spacing(mask))[seq_len(3L)]
  voxel_measure <- if (any(active)) {
    prod(voxel_spacing[active])
  } else {
    min(voxel_spacing)
  }
  scale <- (length(mask_idx) * voxel_measure / K)^(1 / dimension)
  list(scale = as.numeric(scale), dimension = as.integer(dimension))
}

.g3s_cover_components <- function(seeds, membership, feature_mat,
                                  neighbor_indices, use_cosine) {
  component_ids <- sort(unique(membership))
  seed_components <- membership[seeds]
  uncovered <- setdiff(component_ids, seed_components)
  if (!length(uncovered)) return(sort(as.integer(seeds)))

  gradient <- if (use_cosine) {
    calculate_local_gradient(t(feature_mat), neighbor_indices)
  } else {
    vapply(seq_len(nrow(feature_mat)), function(i) {
      adjacent <- neighbor_indices[i, ]
      adjacent <- adjacent[adjacent > 0L]
      if (!length(adjacent)) return(0)
      mean(rowSums(
        (feature_mat[adjacent, , drop = FALSE] -
           matrix(feature_mat[i, ], nrow = length(adjacent),
                  ncol = ncol(feature_mat), byrow = TRUE))^2
      ))
    }, numeric(1L))
  }

  for (component in uncovered) {
    seed_components <- membership[seeds]
    counts <- table(factor(seed_components, levels = component_ids))
    donor_positions <- which(counts[match(seed_components, component_ids)] > 1L)
    if (!length(donor_positions)) {
      stop("cluster4d_g3s: cannot cover every mask component", call. = FALSE)
    }
    remove_position <- donor_positions[
      order(-gradient[seeds[donor_positions]], -seeds[donor_positions])[1L]
    ]
    candidates <- which(membership == component)
    replacement <- candidates[order(gradient[candidates], candidates)[1L]]
    seeds[remove_position] <- replacement
  }
  sort(as.integer(seeds))
}


#' Print Method for G3S Results
#' 
#' @param x A g3s_result object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.g3s_result <- function(x, ...) {
  cat("G3S Clustering Result\n")
  cat("=====================\n\n")

  cat("Clusters:\n")
  cat("  Requested: ", x$parameters$K_requested, "\n")
  cat("  Actual:    ", x$n_clusters, "\n\n")

  cat("Compression:\n")
  cat("  Components: ", x$parameters$n_components, "\n")
  cat("  Variance:   ", round(x$parameters$variance_explained * 100, 1), "%\n")
  cat("  Ratio:      ", round(x$metadata$compression_ratio, 1), "x\n\n")

  cat("Parameters:\n")
  cat("  Alpha:       ", x$parameters$alpha, "\n")
  cat("  Compactness: ", round(x$parameters$compactness, 2), "\n")
  cat("  Refinement:  ", x$parameters$max_refinement_iter, " iterations\n\n")

  cat("Cluster sizes: ")
  cat(paste(summary(as.integer(table(x$cluster))), collapse = " "))
  cat("\n")

  invisible(x)
}
