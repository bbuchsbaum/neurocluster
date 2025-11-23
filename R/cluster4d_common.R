#' Common Utilities for 4D Clustering Algorithms
#'
#' Internal functions shared across all cluster4d methods for consistency
#' and code reuse.
#'
#' @keywords internal
#' @name cluster4d_common
validate_cluster4d_inputs <- function(vec, mask, n_clusters, method = "cluster4d") {
  # Check vec type
  if (!inherits(vec, "NeuroVec") && !inherits(vec, "SparseNeuroVec")) {
    stop(method, ": vec must be a NeuroVec or SparseNeuroVec object")
  }
  
  # Check mask type
  if (!inherits(mask, "NeuroVol")) {
    stop(method, ": mask must be a NeuroVol object")
  }
  
  # Check spatial dimensions match
  vec_dims <- dim(vec)[1:3]
  mask_dims <- dim(mask)
  if (!identical(vec_dims, mask_dims)) {
    stop(method, ": vec and mask must have identical spatial dimensions. ",
         "vec: ", paste(vec_dims, collapse="x"), 
         ", mask: ", paste(mask_dims, collapse="x"))
  }
  
  # Check n_clusters
  if (!is.numeric(n_clusters) || length(n_clusters) != 1) {
    stop(method, ": n_clusters must be a single numeric value")
  }
  
  if (n_clusters <= 0) {
    stop(method, ": n_clusters must be positive, got: ", n_clusters)
  }
  
  # Check mask has valid voxels
  mask_idx <- which(mask > 0)
  if (length(mask_idx) == 0) {
    stop(method, ": No nonzero voxels in mask")
  }
  
  # Check n_clusters vs voxel count
  if (n_clusters > length(mask_idx)) {
    stop(method, sprintf(": Cannot create %d clusters from %d masked voxels", 
                        n_clusters, length(mask_idx)))
  }
  
  invisible(NULL)
}

prepare_cluster4d_data <- function(vec, mask, 
                                  scale_features = TRUE,
                                  scale_coords = FALSE) {
  
  # Get mask indices
  mask_idx <- which(mask > 0)
  n_voxels <- length(mask_idx)
  
  # Extract time series - series returns T x N
  features <- series(vec, mask_idx)
  # Transpose to N x T for consistency
  features <- t(as.matrix(features))
  
  # Scale features if requested
  if (scale_features) {
    features <- scale(features, center = TRUE, scale = TRUE)
    # Replace any NA values (from constant columns) with 0
    features[is.na(features)] <- 0
  }
  
  # Get spatial coordinates in mm
  coords <- index_to_coord(mask, mask_idx)
  coords <- as.matrix(coords)
  
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
  
  # Compute centers for each cluster
  for (i in seq_along(unique_labels)) {
    label <- unique_labels[i]
    cluster_mask <- labels == label
    
    if (sum(cluster_mask) == 0) next
    
    if (method == "mean") {
      # Mean centers
      if (sum(cluster_mask) == 1) {
        centers[i, ] <- features[cluster_mask, ]
        coord_centers[i, ] <- coords[cluster_mask, ]
      } else {
        centers[i, ] <- colMeans(features[cluster_mask, , drop = FALSE])
        coord_centers[i, ] <- colMeans(coords[cluster_mask, , drop = FALSE])
      }
    } else if (method == "medoid") {
      # Medoid centers (point closest to mean)
      cluster_features <- features[cluster_mask, , drop = FALSE]
      cluster_coords <- coords[cluster_mask, , drop = FALSE]
      
      if (nrow(cluster_features) == 1) {
        centers[i, ] <- cluster_features
        coord_centers[i, ] <- cluster_coords
      } else {
        # Find medoid in feature space
        cluster_mean <- colMeans(cluster_features)
        distances <- apply(cluster_features, 1, function(x) sum((x - cluster_mean)^2))
        medoid_idx <- which.min(distances)
        
        centers[i, ] <- cluster_features[medoid_idx, ]
        coord_centers[i, ] <- cluster_coords[medoid_idx, ]
      }
    }
  }
  
  list(
    centers = centers,
    coord_centers = coord_centers,
    n_clusters = n_clusters
  )
}

create_cluster4d_result <- function(labels, mask, data_prep, 
                                   method, parameters,
                                   metadata = list(),
                                   compute_centers = TRUE,
                                   center_method = "mean") {
  
  # Create ClusteredNeuroVol
  clusvol <- ClusteredNeuroVol(mask > 0, clusters = labels)
  
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
  } else if (method == "snic" || method == "slic") {
    if ("compactness" %in% names(params)) {
      if (!"spatial_weight" %in% names(params)) {
        # Higher compactness = more spatial weight
        # Map roughly: compactness 1-10 -> spatial_weight 0.1-0.9
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