#' Initialization Methods for 4D Clustering
#'
#' Functions for initializing cluster seeds using various strategies.
#'
#' @keywords internal
#' @name cluster4d_init

initialize_clusters <- function(coords, features, n_clusters,
                              method = c("gradient", "kmeans", "poisson", "grid", "random"),
                              mask = NULL, vec = NULL, ...) {
  method <- match.arg(method)
  
  n_voxels <- nrow(coords)
  
  # Validate n_clusters
  if (n_clusters > n_voxels) {
    stop("Cannot initialize ", n_clusters, " clusters from ", n_voxels, " voxels")
  }
  
  if (n_clusters == n_voxels) {
    # Trivial case - each voxel is its own cluster
    return(list(
      seeds = 1:n_voxels,
      initial_labels = 1:n_voxels,
      seed_coords = coords
    ))
  }
  
  # Dispatch to appropriate method
  result <- switch(method,
    gradient = init_gradient_seeds(coords, n_clusters, mask, vec, ...),
    kmeans = init_kmeans_seeds(coords, features, n_clusters, ...),
    poisson = init_poisson_seeds(coords, n_clusters, ...),
    grid = init_grid_seeds(coords, n_clusters, ...),
    random = init_random_seeds(coords, n_clusters, ...)
  )
  
  result
}

#' Gradient-based seed initialization
#'
#' Finds initial cluster centers using spatial gradient information.
#' Seeds are placed at locations with high gradient and good spatial separation.
#'
#' @param coords Coordinate matrix (n_voxels x 3)
#' @param n_clusters Number of seeds to find
#' @param mask NeuroVol mask (required for gradient)
#' @param vec NeuroVec data (required for gradient)
#' @param gradient_type Type of gradient: "spatial", "correlation", "intensity"
#'
#' @return List with seeds, initial_labels, and seed_coords
#' @keywords internal
init_gradient_seeds <- function(coords, n_clusters, mask = NULL, vec = NULL,
                               gradient_type = "spatial") {
  
  if (is.null(mask) || is.null(vec)) {
    warning("Gradient initialization requires mask and vec, falling back to kmeans")
    return(init_kmeans_seeds(coords, NULL, n_clusters))
  }
  
  # Compute gradient
  if (gradient_type == "spatial") {
    refvol <- vec[[1]]
    grad <- spatial_gradient(refvol, mask)
  } else if (gradient_type == "correlation") {
    img4d <- as.array(vec)
    mask_numeric <- array(as.numeric(as.logical(mask)), dim = dim(mask))
    grad3d <- correlation_gradient_cpp(img4d, mask_numeric)
    dim(grad3d) <- dim(mask)
    grad <- neuroim2::vol(grad3d, space = space(mask))
  } else {
    # intensity gradient
    refvol <- vec[[1]]
    grad <- spatial_gradient(refvol, mask)
  }
  
  # Get gradient values at masked voxels
  mask_idx <- which(mask > 0)
  grad_vals <- grad[mask_idx]
  
  # Find high-gradient points with spatial separation
  seeds <- find_gradient_seeds(coords, grad_vals, n_clusters)
  
  # Run kmeans from these seeds for initial labels
  seed_coords <- coords[seeds, , drop = FALSE]
  kres <- stats::kmeans(coords, centers = seed_coords, iter.max = 500)
  
  list(
    seeds = seeds,
    initial_labels = kres$cluster,
    seed_coords = seed_coords
  )
}

#' Find gradient-based seed points
#'
#' Internal function that implements the gradient-based seed selection
#' algorithm (previously find_initial_points).
#'
#' @param coords Coordinate matrix
#' @param grad_vals Gradient values at each voxel
#' @param K Number of seeds to find
#'
#' @return Vector of selected seed indices
#' @keywords internal
find_gradient_seeds <- function(coords, grad_vals, K) {
  n_voxels <- nrow(coords)
  
  if (K >= n_voxels) {
    return(1:n_voxels)
  }
  
  # Normalize gradient values
  grad_vals <- (grad_vals - min(grad_vals)) / diff(range(grad_vals))
  
  nfound <- 0
  batchsize <- min(10, K)
  sample_size <- min(1000, n_voxels)
  sel <- c()
  iter <- 1
  
  while (nfound < K) {
    # Sample candidate points
    cand <- sort(sample(1:n_voxels, min(sample_size, n_voxels)))
    gvals <- grad_vals[cand]
    
    if (iter == 1) {
      # First iteration - find points with high gradient and spatial separation
      d <- RANN::nn2(coords[cand, ])
      ord <- order(d$nn.dists[, 2] * (1 - gvals), decreasing = TRUE)
      selected <- cand[ord[1:min(batchsize, length(ord))]]
    } else {
      # Later iterations - also consider distance to already selected points
      d <- RANN::nn2(rbind(coords[cand, ], coords[sel, ]), coords[cand, ])
      ord <- order(d$nn.dists[, 2] * (1 - gvals), decreasing = TRUE)
      selected <- cand[ord[1:min(batchsize, length(ord))]]
    }
    
    sel <- c(sel, selected)
    nfound <- length(sel)
    
    if (nfound >= K) {
      sel <- sel[1:K]
      break
    }
    
    iter <- iter + 1
  }
  
  sel
}

#' K-means based seed initialization
#'
#' Uses k-means clustering to find initial seeds.
#'
#' @param coords Coordinate matrix
#' @param features Feature matrix (optional)
#' @param n_clusters Number of clusters
#' @param use_features Whether to use features in k-means
#'
#' @return List with seeds, initial_labels, and seed_coords
#' @keywords internal
init_kmeans_seeds <- function(coords, features = NULL, n_clusters, 
                             use_features = FALSE) {
  
  n_voxels <- nrow(coords)
  
  # Decide what to cluster on
  if (use_features && !is.null(features)) {
    # Combine coordinates and features
    combined <- cbind(scale(coords), scale(features))
    cluster_data <- combined
  } else {
    # Use only coordinates
    cluster_data <- coords
  }
  
  # Choose initial centers evenly along the data
  init_centers <- cluster_data[as.integer(seq(1, n_voxels, length.out = n_clusters)), , drop = FALSE]
  
  # Run k-means
  kres <- stats::kmeans(cluster_data, centers = init_centers, iter.max = 500)
  
  # Find seeds as points closest to cluster centers
  seeds <- integer(n_clusters)
  for (i in 1:n_clusters) {
    cluster_mask <- kres$cluster == i
    if (sum(cluster_mask) > 0) {
      cluster_coords <- coords[cluster_mask, , drop = FALSE]
      center_coord <- kres$centers[i, 1:3]  # First 3 dims are spatial
      distances <- rowSums((cluster_coords - matrix(center_coord, nrow = sum(cluster_mask), 
                                                   ncol = 3, byrow = TRUE))^2)
      cluster_indices <- which(cluster_mask)
      seeds[i] <- cluster_indices[which.min(distances)]
    }
  }
  
  # Remove any zeros (unfilled seeds)
  seeds <- seeds[seeds > 0]
  
  list(
    seeds = seeds,
    initial_labels = kres$cluster,
    seed_coords = coords[seeds, , drop = FALSE]
  )
}

#' Poisson disk sampling for seed initialization
#'
#' Places seeds with guaranteed minimum separation distance.
#'
#' @param coords Coordinate matrix
#' @param n_clusters Approximate number of clusters
#' @param min_distance Minimum distance between seeds (auto-computed if NULL)
#'
#' @return List with seeds, initial_labels, and seed_coords
#' @keywords internal
init_poisson_seeds <- function(coords, n_clusters, min_distance = NULL) {
  
  n_voxels <- nrow(coords)
  
  # Estimate minimum distance if not provided
  if (is.null(min_distance)) {
    # Estimate based on volume and desired cluster count
    bbox_min <- apply(coords, 2, min)
    bbox_max <- apply(coords, 2, max)
    volume <- prod(bbox_max - bbox_min)
    # Rough estimate: each cluster occupies volume/n_clusters
    cluster_volume <- volume / n_clusters
    min_distance <- (cluster_volume^(1/3)) * 0.7  # 70% of ideal spacing
  }
  
  # Poisson disk sampling
  seeds <- c()
  candidates <- 1:n_voxels
  
  while (length(seeds) < n_clusters && length(candidates) > 0) {
    # Pick a random candidate
    idx <- sample(length(candidates), 1)
    new_seed <- candidates[idx]
    
    if (length(seeds) == 0) {
      # First seed
      seeds <- c(seeds, new_seed)
    } else {
      # Check distance to existing seeds
      seed_coords <- coords[seeds, , drop = FALSE]
      new_coord <- coords[new_seed, ]
      distances <- sqrt(rowSums((seed_coords - matrix(new_coord, nrow = length(seeds), 
                                                     ncol = 3, byrow = TRUE))^2))
      
      if (min(distances) >= min_distance) {
        seeds <- c(seeds, new_seed)
      }
    }
    
    # Remove from candidates
    candidates <- candidates[-idx]
    
    # Relax constraint if we're not finding enough seeds
    if (length(candidates) < 100 && length(seeds) < n_clusters / 2) {
      min_distance <- min_distance * 0.9
    }
  }
  
  # If we didn't get enough seeds, add random ones
  if (length(seeds) < n_clusters) {
    remaining <- setdiff(1:n_voxels, seeds)
    if (length(remaining) > 0) {
      additional <- sample(remaining, min(n_clusters - length(seeds), length(remaining)))
      seeds <- c(seeds, additional)
    }
  }
  
  # Run k-means from these seeds
  seed_coords <- coords[seeds, , drop = FALSE]
  if (length(seeds) < n_clusters) {
    # Pad with duplicates if needed
    seed_coords <- rbind(seed_coords, 
                         seed_coords[sample(nrow(seed_coords), n_clusters - nrow(seed_coords), 
                                          replace = TRUE), ])
  }
  
  kres <- stats::kmeans(coords, centers = seed_coords[1:n_clusters, ], iter.max = 500)
  
  list(
    seeds = seeds,
    initial_labels = kres$cluster,
    seed_coords = coords[seeds, , drop = FALSE]
  )
}

#' Grid-based seed initialization
#'
#' Places seeds on a regular grid within the mask.
#'
#' @param coords Coordinate matrix
#' @param n_clusters Number of clusters
#'
#' @return List with seeds, initial_labels, and seed_coords
#' @keywords internal
init_grid_seeds <- function(coords, n_clusters) {
  
  n_voxels <- nrow(coords)
  
  # Create a regular grid in the bounding box
  bbox_min <- apply(coords, 2, min)
  bbox_max <- apply(coords, 2, max)
  
  # Calculate grid dimensions to get approximately n_clusters points
  grid_points_per_dim <- ceiling(n_clusters^(1/3))
  
  # Create grid
  x_seq <- seq(bbox_min[1], bbox_max[1], length.out = grid_points_per_dim)
  y_seq <- seq(bbox_min[2], bbox_max[2], length.out = grid_points_per_dim)
  z_seq <- seq(bbox_min[3], bbox_max[3], length.out = grid_points_per_dim)
  
  grid <- expand.grid(x = x_seq, y = y_seq, z = z_seq)
  
  # Find nearest voxel to each grid point
  nn_result <- RANN::nn2(coords, grid[, c("x", "y", "z")], k = 1)
  seeds <- unique(nn_result$nn.idx[, 1])
  
  # Limit to n_clusters
  if (length(seeds) > n_clusters) {
    seeds <- seeds[1:n_clusters]
  }
  
  # Run k-means from these seeds
  seed_coords <- coords[seeds, , drop = FALSE]
  if (length(seeds) < n_clusters) {
    # Pad with random points if needed
    remaining <- setdiff(1:n_voxels, seeds)
    if (length(remaining) > 0) {
      additional <- sample(remaining, min(n_clusters - length(seeds), length(remaining)))
      seeds <- c(seeds, additional)
      seed_coords <- coords[seeds, , drop = FALSE]
    }
  }
  
  kres <- stats::kmeans(coords, centers = seed_coords[1:min(n_clusters, length(seeds)), ], 
                        iter.max = 500)
  
  list(
    seeds = seeds,
    initial_labels = kres$cluster,
    seed_coords = coords[seeds, , drop = FALSE]
  )
}

#' Random seed initialization
#'
#' Randomly selects seed points.
#'
#' @param coords Coordinate matrix
#' @param n_clusters Number of clusters
#'
#' @return List with seeds, initial_labels, and seed_coords
#' @keywords internal
init_random_seeds <- function(coords, n_clusters) {
  
  n_voxels <- nrow(coords)
  
  # Random selection
  seeds <- sample(1:n_voxels, min(n_clusters, n_voxels))
  
  # Run k-means from these seeds
  seed_coords <- coords[seeds, , drop = FALSE]
  kres <- stats::kmeans(coords, centers = seed_coords, iter.max = 500)
  
  list(
    seeds = seeds,
    initial_labels = kres$cluster,
    seed_coords = seed_coords
  )
}