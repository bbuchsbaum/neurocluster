#' Common Utilities for 4D Clustering Algorithms
#'
#' Internal functions shared across all cluster4d methods for consistency
#' and code reuse.
#'
#' @keywords internal
#' @name cluster4d_common
ensure_neurovec <- function(vec) {
  # Accept already-valid inputs
  if (inherits(vec, "NeuroVec") || inherits(vec, "SparseNeuroVec")) {
    return(vec)
  }

 # Allow a single-volume NeuroVol and wrap it as a 1-frame NeuroVec
  # This enables supervoxel/clustering algorithms to work with 3D structural images
  if (inherits(vec, "NeuroVol")) {
    arr3d <- as.array(vec)
    arr4d <- array(arr3d, dim = c(dim(vec), 1L))
    # Must use add_dim to create proper 4D NeuroSpace
    space4d <- neuroim2::add_dim(neuroim2::space(vec), 1L)
    return(neuroim2::NeuroVec(arr4d, space4d))
  }

  stop("vec must be a NeuroVec, SparseNeuroVec, or NeuroVol object")
}

.cluster4d_scalar_number <- function(x, name, method,
                                     lower = -Inf, upper = Inf,
                                     integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(method, ": ", name, " must be a finite numeric scalar", call. = FALSE)
  }
  if (integer && x != floor(x)) {
    stop(method, ": ", name, " must be an integer-valued scalar", call. = FALSE)
  }
  if (x < lower || x > upper) {
    interval <- paste0("[", lower, ", ", upper, "]")
    stop(method, ": ", name, " must be in ", interval, ", got: ", x,
         call. = FALSE)
  }
  if (integer) as.integer(x) else as.numeric(x)
}

.cluster4d_scalar_logical <- function(x, name, method) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(method, ": ", name, " must be TRUE or FALSE", call. = FALSE)
  }
  isTRUE(x)
}

# A voxel is included exactly when its mask value is finite and strictly
# positive. Non-finite mask values are rejected rather than silently changing
# the declared voxel set; negative and zero values are valid exclusions.
cluster4d_mask_array <- function(mask, method = "cluster4d") {
  if (!inherits(mask, "NeuroVol")) {
    stop(method, ": mask must be a NeuroVol object", call. = FALSE)
  }
  values <- as.array(mask)
  if (any(!is.finite(values))) {
    stop(
      method,
      ": mask values must be finite; included voxels are exactly values > 0",
      call. = FALSE
    )
  }
  array(values > 0, dim = dim(values))
}

.cluster4d_geometry <- function(x) {
  sp <- neuroim2::space(x)
  list(
    dimensions = as.integer(dim(x)[seq_len(3L)]),
    spacing = as.numeric(neuroim2::spacing(sp)[seq_len(3L)]),
    affine = unname(as.matrix(neuroim2::trans(sp)))
  )
}

.cluster4d_index_to_coord <- function(mask, indices) {
  grid <- neuroim2::index_to_grid(mask, as.integer(indices))
  grid <- matrix(as.numeric(grid), ncol = 3L)
  coords <- neuroim2::grid_to_coord(mask, grid)
  matrix(as.numeric(coords), ncol = 3L)
}

.cluster4d_same_numeric <- function(x, y, tolerance = 1e-8) {
  identical(dim(x), dim(y)) &&
    length(x) == length(y) &&
    all(is.finite(x)) && all(is.finite(y)) &&
    max(abs(x - y), 0) <= tolerance
}

#' Declared capabilities of each unified cluster4d method
#'
#' Internal, fail-closed capability table used to reject common parameters that
#' a method does not implement. A NULL connectivity set means the unified
#' connectivity argument is inactive for that method.
#'
#' @param method Unified cluster4d method name.
#' @return A list describing supported common parameters.
#' @keywords internal
#' @rdname cluster4d_common
cluster4d_method_contract <- function(method) {
  contracts <- list(
    supervoxels = list(connectivity = c(6L, 18L, 26L, 27L), iterations = TRUE, parallel = TRUE, spatial_weight = TRUE),
    snic = list(connectivity = NULL, iterations = FALSE, parallel = FALSE, spatial_weight = TRUE),
    slic = list(connectivity = c(6L, 26L), iterations = TRUE, parallel = TRUE, spatial_weight = TRUE),
    corr_slic = list(connectivity = c(6L, 26L), iterations = TRUE, parallel = FALSE, spatial_weight = TRUE),
    brs_slic = list(connectivity = c(6L, 26L), iterations = TRUE, parallel = FALSE, spatial_weight = TRUE),
    slice_msf = list(connectivity = c(6L, 26L), iterations = FALSE, parallel = FALSE, spatial_weight = TRUE),
    flash3d = list(connectivity = NULL, iterations = TRUE, parallel = FALSE, spatial_weight = TRUE),
    g3s = list(connectivity = c(6L, 18L, 26L), iterations = TRUE, parallel = FALSE, spatial_weight = TRUE),
    rena = list(connectivity = c(6L, 18L, 26L), iterations = TRUE, parallel = FALSE, spatial_weight = FALSE),
    rena_plus = list(connectivity = c(6L, 18L, 26L), iterations = TRUE, parallel = FALSE, spatial_weight = FALSE),
    mcl = list(connectivity = c(6L, 18L, 26L), iterations = TRUE, parallel = FALSE, spatial_weight = TRUE),
    acsc = list(connectivity = NULL, iterations = TRUE, parallel = FALSE, spatial_weight = TRUE),
    commute = list(connectivity = NULL, iterations = FALSE, parallel = FALSE, spatial_weight = TRUE)
  )
  contract <- contracts[[method]]
  if (is.null(contract)) {
    stop("Unknown cluster4d method contract: ", method, call. = FALSE)
  }
  contract
}

validate_cluster4d_inputs <- function(vec, mask, n_clusters, method = "cluster4d",
                                      bad_data_policy = c("error")) {
  bad_data_policy <- match.arg(bad_data_policy)

  # Check vec type
  if (!inherits(vec, "NeuroVec") && !inherits(vec, "SparseNeuroVec")) {
    stop(method, ": vec must be a NeuroVec or SparseNeuroVec object", call. = FALSE)
  }
  
  # Check mask type
  if (!inherits(mask, "NeuroVol")) {
    stop(method, ": mask must be a NeuroVol object", call. = FALSE)
  }

  vec_geometry <- .cluster4d_geometry(vec)
  mask_geometry <- .cluster4d_geometry(mask)
  if (!identical(vec_geometry$dimensions, mask_geometry$dimensions)) {
    stop(method, ": vec and mask must have identical spatial dimensions. ",
         "vec: ", paste(vec_geometry$dimensions, collapse="x"),
         ", mask: ", paste(mask_geometry$dimensions, collapse="x"),
         call. = FALSE)
  }

  if (!.cluster4d_same_numeric(vec_geometry$spacing, mask_geometry$spacing)) {
    stop(method, ": vec and mask must have identical spatial spacing", call. = FALSE)
  }

  if (!.cluster4d_same_numeric(vec_geometry$affine, mask_geometry$affine)) {
    stop(method, ": vec and mask must have identical spatial affine transforms",
         call. = FALSE)
  }

  included <- cluster4d_mask_array(mask, method)
  mask_idx <- which(included)
  if (length(mask_idx) == 0) {
    stop(method, ": mask contains no finite positive voxels", call. = FALSE)
  }

  n_clusters <- .cluster4d_scalar_number(
    n_clusters, "n_clusters", method, integer = TRUE
  )
  if (n_clusters < 1L) {
    stop(method, ": n_clusters must be positive", call. = FALSE)
  }
  if (n_clusters > length(mask_idx)) {
    stop(
      method, ": Cannot create ", n_clusters, " clusters from ",
      length(mask_idx), " masked voxels",
      call. = FALSE
    )
  }

  feature_values <- neuroim2::series(vec, mask_idx)
  bad_features <- !is.finite(feature_values)
  if (bad_data_policy == "error" && any(bad_features)) {
    stop(
      method, ": feature data contain ", sum(bad_features),
      " non-finite value(s) inside the declared mask; bad_data_policy='error'",
      call. = FALSE
    )
  }

  structure(
    list(
      mask = included,
      mask_idx = mask_idx,
      n_voxels = length(mask_idx),
      n_clusters = n_clusters,
      geometry = mask_geometry,
      bad_data_policy = bad_data_policy
    ),
    class = "cluster4d_input_contract"
  )
}

prepare_cluster4d_data <- function(vec, mask, 
                                  scale_features = TRUE,
                                  scale_coords = FALSE) {
  
  # Get mask indices
  mask_idx <- which(cluster4d_mask_array(mask))
  n_voxels <- length(mask_idx)
  
  # Extract time series - series returns T x N
  features <- series(vec, mask_idx)
  if (!is.matrix(features)) {
    # A single-frame NeuroVec is simplified to a length-N vector by neuroim2.
    # Restore the documented T x N shape before transposing to voxel rows.
    features <- matrix(features, nrow = 1L)
  }
  # Transpose to N x T for consistency
  features <- t(as.matrix(features))
  
  # Scale features if requested
  if (scale_features) {
    features <- scale(features, center = TRUE, scale = TRUE)
    # Replace any NA values (from constant columns) with 0
    features[is.na(features)] <- 0
  }
  
  # Get spatial coordinates in mm
  coords <- .cluster4d_index_to_coord(mask, mask_idx)
  
  # Optionally normalize coordinates
  if (scale_coords) {
    coords <- scale(coords, center = TRUE, scale = TRUE)
    coords[is.na(coords)] <- 0
  }
  
  list(
    features = features,
    coords = coords,
    mask_idx = mask_idx,
    n_voxels = n_voxels,
    dims = dim(mask),
    spacing = spacing(mask)
  )
}

compute_cluster_centers <- function(labels, features, coords, method = "mean") {
  
  # Handle NA and invalid labels
  valid_mask <- !is.na(labels) & labels > 0
  labels <- labels[valid_mask]
  features <- features[valid_mask, , drop = FALSE]
  coords <- coords[valid_mask, , drop = FALSE]
  
  # Get unique labels
  unique_labels <- sort(unique(labels))
  n_clusters <- length(unique_labels)
  
  if (n_clusters == 0) {
    warning("No valid clusters found")
    return(list(
      centers = matrix(nrow = 0, ncol = ncol(features)),
      coord_centers = matrix(nrow = 0, ncol = 3),
      n_clusters = 0
    ))
  }
  
  # Initialize center matrices
  centers <- matrix(0, nrow = n_clusters, ncol = ncol(features))
  coord_centers <- matrix(0, nrow = n_clusters, ncol = 3)
  
  if (method == "mean") {
    # Fast path: compute group means using rowsum() (C-accelerated) instead of per-cluster loops.
    grp <- factor(labels, levels = unique_labels)
    counts <- as.numeric(tabulate(as.integer(grp), nbins = n_clusters))
    counts[counts == 0] <- NA_real_

    # Ensure row order matches factor level order (and `counts`) by reordering.
    centers_sum <- rowsum(features, grp, reorder = TRUE)
    coord_sum <- rowsum(coords, grp, reorder = TRUE)

    centers <- centers_sum / counts
    coord_centers <- coord_sum / counts
    # rowsum() attaches dimnames; keep outputs consistent with other methods.
    dimnames(centers) <- NULL
    dimnames(coord_centers) <- NULL
  } else if (method == "medoid") {
    # Compute centers for each cluster
    for (i in seq_along(unique_labels)) {
      label <- unique_labels[i]
      cluster_mask <- labels == label
      if (sum(cluster_mask) == 0) next

      # Medoid centers (point closest to mean)
      cluster_features <- features[cluster_mask, , drop = FALSE]
      cluster_coords <- coords[cluster_mask, , drop = FALSE]

      if (nrow(cluster_features) == 1) {
        centers[i, ] <- cluster_features
        coord_centers[i, ] <- cluster_coords
      } else {
        cluster_mean <- colMeans(cluster_features)
        distances <- apply(cluster_features, 1, function(x) sum((x - cluster_mean)^2))
        medoid_idx <- which.min(distances)

        centers[i, ] <- cluster_features[medoid_idx, ]
        coord_centers[i, ] <- cluster_coords[medoid_idx, ]
      }
    }
  } else {
    stop("Unknown center method: ", method)
  }
  
  list(
    centers = centers,
    coord_centers = coord_centers,
    n_clusters = n_clusters
  )
}

.cluster4d_original_data <- function(vec, mask, method = "cluster4d") {
  input <- validate_cluster4d_inputs(vec, mask, 1L, method)
  feature_values <- neuroim2::series(vec, input$mask_idx)
  features <- if (is.matrix(feature_values)) {
    t(as.matrix(feature_values))
  } else {
    matrix(as.numeric(feature_values), nrow = length(input$mask_idx))
  }
  coords <- .cluster4d_index_to_coord(mask, input$mask_idx)
  dimnames(features) <- NULL
  dimnames(coords) <- NULL
  list(
    features = features,
    coords = coords,
    mask_idx = input$mask_idx,
    geometry = input$geometry
  )
}

.cluster4d_merge_parameters <- function(existing, canonical) {
  if (!is.list(existing)) existing <- list()
  if (!is.list(canonical)) canonical <- list()
  for (name in names(canonical)) {
    if (!is.null(canonical[[name]])) existing[[name]] <- canonical[[name]]
  }
  existing
}

# Canonical public boundary for every unified cluster4d method. Method-specific
# centers are retained only in metadata; the public centers are always means of
# final labels in the original feature space.
finalize_cluster4d_result <- function(result, vec, mask, method,
                                      parameters = list()) {
  data <- .cluster4d_original_data(
    vec, mask, paste0("cluster4d:", method, ":result")
  )

  raw_labels <- if (!is.null(result$labels)) result$labels else result$cluster
  if (!is.numeric(raw_labels) || length(raw_labels) != nrow(data$features) ||
      any(!is.finite(raw_labels)) || any(raw_labels != floor(raw_labels)) ||
      any(raw_labels <= 0)) {
    stop(
      "cluster4d:", method,
      ": result labels must be finite positive integers, one per masked voxel",
      call. = FALSE
    )
  }

  original_label_ids <- sort(unique(as.integer(raw_labels)))
  labels <- as.integer(match(as.integer(raw_labels), original_label_ids))
  actual_k <- as.integer(length(original_label_ids))
  label_ids <- seq_len(actual_k)

  center_info <- compute_cluster_centers(
    labels, data$features, data$coords, method = "mean"
  )
  centers <- unname(as.matrix(center_info$centers))
  coord_centers <- unname(as.matrix(center_info$coord_centers))

  included_mask <- neuroim2::NeuroVol(
    cluster4d_mask_array(mask, paste0("cluster4d:", method, ":result")),
    space = neuroim2::space(mask)
  )
  clusvol <- suppressWarnings(
    neuroim2::ClusteredNeuroVol(included_mask, clusters = labels)
  )

  merged_parameters <- .cluster4d_merge_parameters(
    result$parameters, parameters
  )
  if (is.null(merged_parameters$n_clusters_requested)) {
    merged_parameters$n_clusters_requested <- if (!is.null(parameters$n_clusters)) {
      as.integer(parameters$n_clusters)
    } else actual_k
  }

  provenance <- list(
    schema_version = "1.0.0",
    label_space = list(
      ids = label_ids,
      row_to_label = label_ids,
      original_ids = original_label_ids,
      relabeled = !identical(original_label_ids, label_ids)
    ),
    feature_space = list(
      representation = "original",
      summary = "mean",
      n_features = as.integer(ncol(data$features)),
      source = "masked voxel series from vec"
    ),
    coordinate_space = list(
      units = "mm",
      summary = "mean",
      source = "mask affine",
      geometry = data$geometry
    ),
    mask = list(
      inclusion_rule = "finite values strictly greater than zero",
      n_voxels = as.integer(nrow(data$features))
    )
  )

  standard_names <- c(
    "labels", "cluster", "clusvol", "centers", "coord_centers",
    "actual_k", "n_clusters", "label_ids", "method", "parameters",
    "provenance", "metadata"
  )
  extras <- result[setdiff(names(result), standard_names)]
  metadata <- if (is.list(result$metadata)) result$metadata else list()
  old_classes <- setdiff(
    class(result), c("cluster4d_result", "cluster_result", "list")
  )

  canonical <- c(
    list(
      labels = labels,
      cluster = labels,
      clusvol = clusvol,
      centers = centers,
      coord_centers = coord_centers,
      actual_k = actual_k,
      n_clusters = actual_k,
      label_ids = label_ids,
      method = method,
      parameters = merged_parameters,
      provenance = provenance,
      metadata = metadata
    ),
    extras
  )
  structure(
    canonical,
    class = unique(c(old_classes, "cluster4d_result", "cluster_result", "list"))
  )
}

create_cluster4d_result <- function(labels, mask, data_prep, 
                                   method, parameters,
                                   metadata = list(),
                                   compute_centers = TRUE,
                                   center_method = "mean") {
  
  # Preserve mask geometry while enforcing the package-wide inclusion rule.
  included_mask <- neuroim2::NeuroVol(
    cluster4d_mask_array(mask, method),
    space = neuroim2::space(mask)
  )
  clusvol <- suppressWarnings(
    ClusteredNeuroVol(included_mask, clusters = labels)
  )
  
  # Compute centers if requested
  if (compute_centers) {
    center_info <- compute_cluster_centers(
      labels, 
      data_prep$features, 
      data_prep$coords,
      method = center_method
    )
    centers <- center_info$centers
    coord_centers <- center_info$coord_centers
    n_clusters <- center_info$n_clusters
  } else {
    # Centers must be provided in metadata
    centers <- metadata$centers
    coord_centers <- metadata$coord_centers
    # Use metadata n_clusters if provided, otherwise calculate from labels
    n_clusters <- if (!is.null(metadata$n_clusters)) {
      metadata$n_clusters
    } else {
      length(unique(labels[!is.na(labels) & labels > 0]))
    }
  }
  
  # Ensure centers are present
  if (is.null(centers) || is.null(coord_centers)) {
    warning("Computing centers as they were not provided")
    center_info <- compute_cluster_centers(
      labels, 
      data_prep$features, 
      data_prep$coords,
      method = center_method
    )
    centers <- center_info$centers
    coord_centers <- center_info$coord_centers
    n_clusters <- center_info$n_clusters
  }
  
  # Create result structure
  result <- structure(
    list(
      clusvol = clusvol,
      cluster = labels,
      centers = centers,
      coord_centers = coord_centers,
      n_clusters = n_clusters,
      method = method,
      parameters = parameters,
      metadata = metadata
    ),
    class = c("cluster4d_result", "cluster_result", "list")
  )
  
  result
}

map_cluster4d_params <- function(method, ...) {
  params <- list(...)
  
  # Common mappings across methods
  if ("K" %in% names(params)) {
    if (!"n_clusters" %in% names(params)) {
      params$n_clusters <- params$K
    }
    params$K <- NULL
  }
  
  # Method-specific mappings
  if (method == "supervoxels") {
    if ("alpha" %in% names(params)) {
      if (!"spatial_weight" %in% names(params)) {
        # alpha is feature weight, so spatial = 1 - alpha
        params$spatial_weight <- 1 - params$alpha
      }
    }
    if ("iterations" %in% names(params)) {
      if (!"max_iterations" %in% names(params)) {
        params$max_iterations <- params$iterations
      }
    }
  } else if (method == "snic") {
    if ("compactness" %in% names(params)) {
      if (!"spatial_weight" %in% names(params)) {
        params$spatial_weight <- params$compactness / 10
      }
    }
    if ("max_iter" %in% names(params)) {
      stop("map_cluster4d_params:snic: max_iter is not supported", call. = FALSE)
    }
  } else if (method == "slic" || method == "corr_slic" || method == "brs_slic") {
    if ("compactness" %in% names(params)) {
      if (!"spatial_weight" %in% names(params)) {
        params$spatial_weight <- min(0.9, params$compactness / 10)
      }
    }
    if ("max_iter" %in% names(params)) {
      if (!"max_iterations" %in% names(params)) {
        params$max_iterations <- params$max_iter
      }
    }
  } else if (method == "flash3d") {
    if ("lambda_s" %in% names(params)) {
      if (!"spatial_weight" %in% names(params)) {
        params$spatial_weight <- params$lambda_s
      }
    }
    if ("rounds" %in% names(params)) {
      if (!"max_iterations" %in% names(params)) {
        params$max_iterations <- params$rounds
      }
    }
  } else if (method == "slice_msf") {
    if ("target_k_global" %in% names(params) && params$target_k_global > 0) {
      if (!"n_clusters" %in% names(params)) {
        params$n_clusters <- params$target_k_global
      }
    }
    if ("compactness" %in% names(params)) {
      if (!"spatial_weight" %in% names(params)) {
        params$spatial_weight <- min(0.9, params$compactness / 10)
      }
    }
  }
  
  params
}

#' Suggest cluster4d parameters based on data characteristics
#'
#' Provides parameter recommendations based on data size and user priorities.
#'
#' @param n_voxels Number of voxels in mask
#' @param n_timepoints Number of time points
#' @param priority What to optimize for: "speed", "quality", "memory", or "balanced"
#'
#' @return A list with suggested parameters for each method
#' 
#' @examples
#' # Get parameter suggestions for a typical fMRI dataset
#' params <- suggest_cluster4d_params(50000, 200, priority = "balanced")
#' print(params$recommended_method)
#' print(params$n_clusters)
#' 
#' # Speed-optimized parameters for large dataset
#' speed_params <- suggest_cluster4d_params(100000, 300, priority = "speed")
#' print(speed_params$recommended_method)
#' 
#' # Quality-optimized parameters for smaller dataset
#' quality_params <- suggest_cluster4d_params(10000, 150, priority = "quality")
#' print(quality_params$n_clusters)
#' 
#' # Memory-efficient parameters
#' memory_params <- suggest_cluster4d_params(200000, 100, priority = "memory")
#' print(memory_params$recommended_method)
#' 
#' @export
suggest_cluster4d_params <- function(n_voxels, n_timepoints, 
                                    priority = c("balanced", "speed", "quality", "memory")) {
  priority <- match.arg(priority)
  
  # Base recommendations
  suggestions <- list()
  
  # Estimate reasonable cluster count (roughly 1 cluster per 100-500 voxels)
  base_k <- max(10, min(1000, n_voxels / 250))
  
  if (priority == "speed") {
    suggestions$recommended_method <- if (n_voxels > 50000) "slice_msf" else "flash3d"
    suggestions$n_clusters <- round(base_k * 0.7)  # Fewer clusters for speed
    suggestions$slice_msf <- list(
      num_runs = 1,
      r = 8,
      min_size = 120
    )
    suggestions$flash3d <- list(
      rounds = 1,
      bits = 64,
      dctM = 8
    )
    suggestions$snic <- list(
      compactness = 7  # Higher for faster convergence
    )
  } else if (priority == "quality") {
    suggestions$recommended_method <- if (n_voxels < 20000) "supervoxels" else "slic"
    suggestions$n_clusters <- round(base_k * 1.2)  # More clusters for quality
    suggestions$supervoxels <- list(
      iterations = 50,
      parallel = TRUE,
      converge_thresh = 0.0005
    )
    suggestions$slic <- list(
      max_iter = 15,
      preserve_k = TRUE,
      seed_relocate = "correlation"
    )
    suggestions$slice_msf <- list(
      num_runs = 5,
      consensus = TRUE,
      use_features = TRUE
    )
  } else if (priority == "memory") {
    suggestions$recommended_method <- "snic"  # Non-iterative, low memory
    suggestions$n_clusters <- round(base_k * 0.8)
    suggestions$snic <- list(
      compactness = 5
    )
    suggestions$slice_msf <- list(
      num_runs = 1,
      r = 6  # Fewer DCT components
    )
  } else {  # balanced
    suggestions$recommended_method <- if (n_voxels > 30000) "flash3d" else "slic"
    suggestions$n_clusters <- round(base_k)
    suggestions$general <- list(
      spatial_weight = 0.5,
      max_iterations = 10,
      connectivity = 26
    )
  }
  
  # Add data size info
  suggestions$data_info <- list(
    n_voxels = n_voxels,
    n_timepoints = n_timepoints,
    estimated_memory_mb = round((n_voxels * n_timepoints * 8) / 1e6)
  )
  
  suggestions
}
