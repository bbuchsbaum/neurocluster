#' Tesselate a Mask Volume into K Clusters using K-means
#'
#' This function tesselates a given mask volume into K clusters using k-means
#' clustering applied to spatial coordinates. It returns a clustered mask volume object.
#'
#' @param mask A \code{NeuroVol} object representing the mask volume.
#' @param K An integer value specifying the number of clusters (default: 100).
#'
#' If \code{K} exceeds the number of nonzero voxels, a warning is issued and \code{K}
#' is set to the number of nonzero voxels.
#'
#' @return An instance of \code{ClusteredNeuroVol} representing the clustered mask volume.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a NeuroVol object 'mask' and you want to create 150 clusters
#' clustered_volume <- tesselate(mask, K = 150)
#' }
#'
#' @export
tesselate <- function(mask, K = 100) {
  mask.idx <- which(mask > 0)
  if (length(mask.idx) == 0) {
    stop("No nonzero voxels in mask.")
  }
  coords <- index_to_coord(mask, mask.idx)
  nvox <- nrow(coords)

  if (K >= nvox) {
    if (nvox == 1) {
      # Special case: only one voxel, return trivial clustering
      cluster_assignments <- rep(1, length(mask.idx))
      clusvol <- ClusteredNeuroVol(mask, cluster_assignments)
      return(clusvol)
    }
    warning("K is greater than or equal to the number of valid voxels. Setting K = ", nvox - 1)
    K <- nvox - 1
  }

  # Choose initial centers evenly along the coordinate matrix
  init_centers <- coords[as.integer(seq(1, nvox, length.out = K)), , drop = FALSE]

  kres <- stats::kmeans(coords, centers = init_centers, iter.max = 500)
  clusvol <- ClusteredNeuroVol(mask, kres$cluster)
  clusvol
}


#' Initialize clusters using gradient-based or random coordinates
#'
#' Internal function to initialize cluster assignments.
#'
#' @keywords internal
#' @noRd
init_cluster <- function(bvec, mask, coords, K, use_gradient = TRUE) {
  mask.idx <- which(mask > 0)
  nvox <- nrow(coords)

  if (K >= nvox) {
    if (nvox == 1) {
      # Special case: only one voxel, return trivial clustering
      return(rep(1, nvox))
    }
    warning("K is greater than or equal to the number of valid voxels. Setting K = ", nvox - 1)
    K <- nvox - 1
  }

  if (use_gradient && nvox >= 100) {  # Only use gradient for larger volumes
    tryCatch({
      refvol <- bvec[[1]]
      grad <- spatial_gradient(refvol, mask)
      grad_vals <- grad[mask.idx]
      valid_coords <- index_to_grid(mask, mask.idx)

      # find_initial_points function for gradient-based initialization
      init <- find_initial_points(valid_coords, grad_vals, K)

      # run kmeans with chosen seeds
      if (length(init$selected) >= K) {
        kres <- stats::kmeans(coords, centers = coords[init$selected[1:K], , drop = FALSE], iter.max = 500)
        return(kres$cluster)
      }
    }, error = function(e) {
      warning("Gradient-based initialization failed, falling back to uniform sampling: ", e$message)
    })
  }
  
  # Fallback: uniform sampling (also used when use_gradient=FALSE)
  init_centers <- coords[as.integer(seq(1, nvox, length.out = K)), , drop = FALSE]
  kres <- stats::kmeans(coords, centers = init_centers, iter.max = 500)
  kres$cluster
}


# Maintain the exact-K invariant used by the iterative supervoxel fit. An
# assignment step can empty a label even when K <= N. Empty labels are repaired
# deterministically by moving a distinct voxel whose current assignment has the
# largest combined feature/spatial loss. Donors are restricted to clusters with
# at least two voxels, so repairing one label cannot empty another.
.supervoxel_reseed_empty <- function(labels, feature_mat, coords, K,
                                     alpha, sigma1, sigma2) {
  labels <- as.integer(labels)
  counts <- tabulate(labels, nbins = K)
  empty_labels <- which(counts == 0L)
  reseeded_voxels <- integer()

  if (length(empty_labels) == 0L) {
    return(list(
      labels = labels,
      counts = counts,
      empty_labels = integer(),
      reseeded_voxels = integer()
    ))
  }

  for (empty_label in empty_labels) {
    donor_voxels <- which(counts[labels] > 1L)
    if (length(donor_voxels) == 0L) {
      stop(
        "supervoxels: exact-K reseeding is infeasible because no cluster has a spare voxel",
        call. = FALSE
      )
    }

    active_labels <- which(counts > 0L)
    feature_centers <- matrix(0, nrow = nrow(feature_mat), ncol = K)
    coord_centers <- matrix(0, nrow = ncol(coords), ncol = K)
    for (cluster in active_labels) {
      members <- which(labels == cluster)
      feature_centers[, cluster] <- rowMeans(
        feature_mat[, members, drop = FALSE]
      )
      coord_centers[, cluster] <- colMeans(
        coords[members, , drop = FALSE]
      )
    }

    donor_labels <- labels[donor_voxels]
    feature_delta <- feature_mat[, donor_voxels, drop = FALSE] -
      feature_centers[, donor_labels, drop = FALSE]
    coord_delta <- t(coords[donor_voxels, , drop = FALSE]) -
      coord_centers[, donor_labels, drop = FALSE]
    feature_loss <- 1 - exp(
      -colSums(feature_delta * feature_delta) / (2 * sigma1^2)
    )
    spatial_loss <- 1 - exp(
      -colSums(coord_delta * coord_delta) / (2 * sigma2^2)
    )
    loss <- alpha * feature_loss + (1 - alpha) * spatial_loss

    # which.max is deterministic and donor_voxels is sorted, so exact ties go
    # to the lowest masked-voxel index.
    chosen <- donor_voxels[which.max(loss)]
    old_label <- labels[chosen]
    labels[chosen] <- empty_label
    counts[old_label] <- counts[old_label] - 1L
    counts[empty_label] <- 1L
    reseeded_voxels <- c(reseeded_voxels, chosen)
  }

  if (any(counts <= 0L) || length(unique(reseeded_voxels)) !=
      length(reseeded_voxels)) {
    stop("supervoxels: internal exact-K reseed invariant failed", call. = FALSE)
  }

  list(
    labels = labels,
    counts = counts,
    empty_labels = as.integer(empty_labels),
    reseeded_voxels = as.integer(reseeded_voxels)
  )
}


.supervoxel_center_state <- function(labels, feature_mat, coords, K,
                                     alpha, sigma1, sigma2,
                                     use_medoid = FALSE,
                                     parallel = FALSE) {
  repair <- .supervoxel_reseed_empty(
    labels, feature_mat, coords, K, alpha, sigma1, sigma2
  )
  labels <- repair$labels

  if (parallel && !use_medoid) {
    center_result <- compute_centroids_parallel_fast(
      labels - 1L, feature_mat, t(coords), as.integer(K)
    )
    feature_centers <- center_result$centers
    coord_centers <- center_result$coord_centers
    counts <- as.integer(center_result$counts)
  } else {
    centers <- compute_centroids(
      feature_mat, coords, labels, medoid = use_medoid
    )
    center_ids <- suppressWarnings(as.integer(names(centers$center)))
    if (length(center_ids) != K || anyNA(center_ids) ||
        !setequal(center_ids, seq_len(K))) {
      stop("supervoxels: center labels do not cover 1:K", call. = FALSE)
    }
    order_idx <- match(seq_len(K), center_ids)
    feature_centers <- t(do.call(rbind, centers$center[order_idx]))
    coord_centers <- t(do.call(rbind, centers$centroid[order_idx]))
    counts <- tabulate(labels, nbins = K)
  }

  if (any(counts <= 0L) || length(counts) != K ||
      ncol(feature_centers) != K || ncol(coord_centers) != K ||
      any(!is.finite(feature_centers)) || any(!is.finite(coord_centers))) {
    stop(
      "supervoxels: center update violated the non-empty exact-K contract",
      call. = FALSE
    )
  }

  repair$counts <- counts
  repair$centers <- unname(feature_centers)
  repair$coord_centers <- unname(coord_centers)
  repair
}


#' Fit Supervoxel Clusters
#'
#' Internal function that performs an iterative, spatially-constrained clustering
#' of voxel features and coordinates.
#'
#' @keywords internal
#' @importFrom assertthat assert_that
#' @import neuroim2
#' @noRd
supervoxel_cluster_fit <- function(feature_mat,
                                   coords,
                                   K = min(500, nrow(coords)),
                                   sigma1 = 1,
                                   sigma2 = 5,
                                   alpha = 0.5,
                                   iterations = 25,
                                   connectivity = 26,
                                   use_medoid = FALSE,
                                   initclus = NULL,
                                   use_gradient = TRUE,
                                   parallel = TRUE,
                                   grain_size = 100,
                                   verbose = FALSE,
                                   converge_thresh = 0.001,
                                   trace = FALSE) {

  # Early parameter validation
  nvox <- nrow(coords)
  if (nvox == 0) {
    stop("No voxels to cluster (coords has 0 rows)")
  }
  if (!is.matrix(feature_mat) || ncol(feature_mat) != nvox ||
      any(!is.finite(feature_mat))) {
    stop("feature_mat must be a finite feature-by-voxel matrix")
  }
  if (!is.matrix(coords) || ncol(coords) < 1L || ncol(coords) > 3L ||
      any(!is.finite(coords))) {
    stop("coords must be a finite voxel-by-dimension matrix with 1 to 3 columns")
  }
  coord_dimensions <- ncol(coords)
  if (coord_dimensions < 3L) {
    coords <- cbind(
      coords,
      matrix(0, nrow = nvox, ncol = 3L - coord_dimensions)
    )
  }
  
  if (K <= 0) {
    stop("K must be positive, got: ", K)
  }
  
  if (K > nvox) {
    stop(sprintf("Cannot create %d clusters from %d voxels. K must be <= number of voxels.", K, nvox))
  }
  
  if (K == nvox) {
    warning(sprintf("K equals number of voxels (%d). Each voxel will be its own cluster.", nvox))
    return(list(
      clusters = seq_len(nvox),
      centers = t(feature_mat),
      coord_centers = coords[, seq_len(coord_dimensions), drop = FALSE],
      counts = rep.int(1L, nvox),
      parallel_used = FALSE,
      reseed_events = list(),
      iteration_trace = list()
    ))
  }
  
  # Adjust connectivity based on number of voxels
  max_connectivity <- min(connectivity, nvox - 1)
  if (max_connectivity != connectivity) {
    warning(sprintf("Connectivity reduced from %d to %d due to limited voxels", connectivity, max_connectivity))
    connectivity <- max_connectivity
  }
  
  # Ensure minimum connectivity
  if (connectivity < 2) {
    warning("Connectivity must be at least 2, setting to 2")
    connectivity <- 2
  }
  
  assert_that(connectivity <= 27)
  assert_that(alpha >= 0 && alpha <= 1)
  assert_that(iterations > 0)
  assert_that(converge_thresh > 0)
  assert_that(is.flag(parallel))
  assert_that(is.flag(trace))
  assert_that(length(grain_size) == 1L, is.finite(grain_size), grain_size > 0)

  # Center and scale the feature matrix
  feature_mat <- base::scale(feature_mat, center = TRUE, scale = TRUE)

  # Check for NA values after scaling
  if (any(is.na(feature_mat))) {
    warning("NA values in feature matrix after scaling, replacing with 0")
    feature_mat[is.na(feature_mat)] <- 0
  }

  # FIX 2: Adaptive sigma1 based on feature dimensionality
  # For z-scored data, ||v1 - v2||^2 = 2*T for orthogonal signals (where T = n_features)
  # We want sigma1 such that exp(-d^2/(2*sigma1^2)) gives meaningful discrimination
  # Setting sigma1 = sqrt(T) gives exp(-1) ≈ 0.37 for orthogonal signals
  n_features <- nrow(feature_mat)
  if (is.null(sigma1)) {
    sigma1 <- sqrt(n_features)
    if (verbose) message(sprintf("supervoxels: sigma1 = %.2f (adaptive, sqrt of %d features)", sigma1, n_features))
  }

  # If no initial clusters, use gradient-based seeding for better initialization
  if (is.null(initclus)) {
    grad_distance <- if (n_features == 1) "euclidean" else "cosine"
    # FIX 5: Gradient-based seeding - place initial seeds in low-gradient regions
    # Low gradient = center of homogeneous feature regions, not boundaries
    # This improves initialization quality and helps the algorithm converge faster

    seeds <- tryCatch({
      if (verbose) {
        message("supervoxels: Using functional gradient-based seeding (", grad_distance, ")...")
      }

      # Compute functional gradient - need features as N x T (rows=voxels, cols=features)
      # feature_mat is T x N, so transpose it
      feature_mat_transposed <- t(feature_mat)

      # Normalize features for cosine-based gradient computation; leave Euclidean
      # path unnormalized so magnitude differences are preserved for 3D inputs.
      feature_mat_for_grad <- if (grad_distance == "cosine") {
        row_norms <- sqrt(rowSums(feature_mat_transposed^2))
        row_norms[row_norms == 0] <- 1  # Avoid division by zero
        feature_mat_transposed / row_norms
      } else {
        feature_mat_transposed
      }

      # Find K seeds in low-gradient (stable, homogeneous) regions
      gradient_seeds <- find_gradient_seeds_g3s(
        feature_mat_for_grad,
        coords,
        K = K,
        k_neighbors = min(26, nvox - 1),  # 26-connectivity or less
        oversample_ratio = 3,
        min_separation_factor = 0.5,  # Allow seeds reasonably close together
        distance = grad_distance
      )

      if (length(gradient_seeds) < K) {
        warning(sprintf("Gradient seeding found only %d of %d seeds, falling back to uniform",
                       length(gradient_seeds), K))
        NULL
      } else {
        gradient_seeds[1:K]
      }
    }, error = function(e) {
      if (verbose) warning("Gradient-based seeding failed: ", e$message, ". Falling back to uniform.")
      NULL
    })

    # If gradient seeding worked, use those seeds; otherwise fall back to uniform grid
    if (!is.null(seeds) && length(seeds) == K) {
      # Use gradient-based seeds as k-means centers
      init_centers <- coords[seeds, , drop = FALSE]
      kres <- stats::kmeans(coords, centers = init_centers, iter.max = 500)

      if (verbose) message("supervoxels: Gradient seeding successful, ", K, " seeds placed in low-gradient regions")
    } else {
      # Fallback: uniform grid sampling
      init_centers <- coords[as.integer(seq(1, nvox, length.out = K)), , drop = FALSE]
      kres <- stats::kmeans(coords, centers = init_centers, iter.max = 500)
    }

    clusid <- sort(unique(kres$cluster))

    # CRITICAL FIX: Remap cluster IDs to be contiguous (1, 2, 3, ..., K)
    # This ensures centroids can be indexed by cluster_id - 1 in C++
    cluster_mapping <- setNames(1:length(clusid), clusid)
    curclus <- cluster_mapping[as.character(kres$cluster)]

  } else {
    # Use user-supplied initclus
    assert_that(length(initclus) == nvox)
    clusid <- sort(unique(initclus))
    assert_that(length(clusid) == K)
    
    # CRITICAL FIX: Ensure user-supplied clusters are also contiguous
    # This maintains consistency with the k-means initialization path
    cluster_mapping <- setNames(1:length(clusid), clusid)
    curclus <- cluster_mapping[as.character(initclus)]

  }

  use_parallel_assignment <- parallel && nvox > 1000L
  use_parallel_centers <- parallel && nvox > 1000L && K > 50L
  center_state <- .supervoxel_center_state(
    curclus, feature_mat, coords, K, alpha, sigma1, sigma2,
    use_medoid = use_medoid,
    parallel = use_parallel_centers
  )
  curclus <- center_state$labels
  cluster_counts <- center_state$counts
  num_centroids <- center_state$centers
  sp_centroids <- center_state$coord_centers

  # Find connectivity-based neighbors
  neib <- FNN::get.knn(coords, k = connectivity)
  # dthresh is the median distance among the 'connectivity'-th neighbor
  dthresh <- stats::median(neib$nn.dist[, connectivity, drop = FALSE])
  if (verbose) message("dthresh: ", dthresh)

  coords_t <- t(coords)

  # Precompute and sanitize neighbor indices once (constant across iterations).
  nn_indices <- neib$nn.index - 1L  # Convert to 0-based
  if (any(nn_indices < 0 | nn_indices >= nvox, na.rm = TRUE)) {
    warning("Invalid neighbor indices found. Coercing to valid range.")
    nn_indices[is.na(nn_indices) | nn_indices < 0] <- 0L
    nn_indices[nn_indices >= nvox] <- nvox - 1L
  }
  nn_dist <- neib$nn.dist

  iter <- 1
  switches <- 1
  iter.max <- iterations
  newclus <- curclus
  prev_switches <- Inf
  no_improvement_count <- 0
  reseed_events <- list()
  iteration_trace <- list()
  
  # FIX 1: Removed cheap_iters - use full alpha from iteration 1

  # The old "cheap-then-exact" strategy (alpha=0 for first iteration) destroyed
  # feature information when k-means initialization didn't align with feature clusters.
  # Now we use full alpha from the start to properly weight features vs spatial.
  original_alpha <- alpha

  # Adaptive spatial binning (small K needs wider search)
  window_factor <- if (K < 20) {
    max(4.0, 3.0 + 10.0 / K)
  } else {
    3.0
  }
  bin_expand <- if (K < 10 || nvox < 5000) 2L else 1L

  grain_size_eff <- as.integer(grain_size)

  while (iter <= iter.max && switches > 0) {
    # FIX 1: Always use full alpha (no more cheap_iters strategy)
    current_alpha <- original_alpha
    
    # FIX 7: Convert curclus to 0-based for C++ (C++ uses 0-based centroid indexing)
    # The C++ code iterates centroids as k=0 to K-1 and stores these 0-based indices
    # in the spatial bins. When it returns assignments, they are also 0-based.
    curclus_0based <- curclus - 1L

    # Use spatially-binned assignment for efficient candidate selection
    if (use_parallel_assignment) {
      newclus_0based <- fused_assignment_parallel_binned(
        nn_indices, nn_dist, curclus_0based,
        coords_t, num_centroids, sp_centroids, cluster_counts, feature_mat,
        dthresh, sigma1, sigma2, current_alpha,
        grain_size = grain_size_eff,
        window_factor = window_factor,
        bin_expand = bin_expand
      )
    } else {
      # Even sequential mode benefits from spatial binning
      newclus_0based <- fused_assignment_binned(
        nn_indices, nn_dist, curclus_0based,
        coords_t, num_centroids, sp_centroids, cluster_counts, feature_mat,
        dthresh, sigma1, sigma2, current_alpha,
        window_factor = window_factor,
        bin_expand = bin_expand
      )
    }

    assigned <- newclus_0based + 1L
    center_state <- .supervoxel_center_state(
      assigned, feature_mat, coords, K, alpha, sigma1, sigma2,
      use_medoid = use_medoid,
      parallel = use_parallel_centers
    )
    newclus <- center_state$labels
    cluster_counts <- center_state$counts
    num_centroids <- center_state$centers
    sp_centroids <- center_state$coord_centers

    if (length(center_state$empty_labels) > 0L) {
      reseed_events[[length(reseed_events) + 1L]] <- list(
        iteration = as.integer(iter),
        labels = center_state$empty_labels,
        voxels = center_state$reseeded_voxels
      )
    }

    # Reseeding is part of the assignment transition, so convergence cannot
    # accept a state with absent labels.
    switches <- sum(newclus != curclus)
    curclus <- newclus
    if (trace) {
      iteration_trace[[length(iteration_trace) + 1L]] <- list(
        iteration = as.integer(iter),
        counts = cluster_counts,
        n_feature_centers = as.integer(ncol(num_centroids)),
        n_coord_centers = as.integer(ncol(sp_centroids)),
        empty_labels = center_state$empty_labels,
        reseeded_voxels = center_state$reseeded_voxels
      )
    }

    # Check for convergence
    switch_ratio <- switches / nvox
    if (verbose) {
      message(
        "supervoxels_fit: iter ",
        iter,
        " -- num switches = ",
        switches,
        " (",
        round(switch_ratio * 100, 2),
        "% of voxels)"
      )
    }
    
    # Early convergence check
    if (switch_ratio < converge_thresh) {
      if (verbose) {
        message("supervoxels_fit: converged at iteration ", iter, " (switch ratio < ", converge_thresh, ")")
      }
      break
    }
    
    # Check if we're stuck (no improvement)
    if (switches >= prev_switches) {
      no_improvement_count <- no_improvement_count + 1
      if (no_improvement_count >= 3) {
        if (verbose) {
          message("supervoxels_fit: stopping - no improvement for 3 iterations")
        }
        break
      }
    } else {
      no_improvement_count <- 0
    }
    
    prev_switches <- switches
    iter <- iter + 1
  }

  # Return the final clusters and the final centroids
  list(
    clusters = curclus,
    centers = t(num_centroids),
    coord_centers = t(sp_centroids)[, seq_len(coord_dimensions), drop = FALSE],
    counts = cluster_counts,
    parallel_used = use_parallel_assignment,
    reseed_events = reseed_events,
    iteration_trace = iteration_trace
  )
}


#' K-nearest-neighbor shrink
#'
#' Replace each voxel by the mean of its k nearest neighbors in its local spatial neighborhood.
#'
#' @param bvec A \code{\link[neuroim2:NeuroVec-class]{NeuroVec}} instance (the data).
#' @param mask A \code{\link[neuroim2:NeuroVol-class]{NeuroVol}} mask defining the voxels to include. If numeric, nonzero = included.
#' @param k The number of nearest neighbors to average over.
#' @param connectivity The number of spatial neighbors to include in the search around each voxel.
#'
#' @return A \code{SparseNeuroVec} or similar object with the smoothed data.
#'
#' @examples
#' \dontrun{
#' mask <- neuroim2::NeuroVol(array(1, c(20,20,20)), neuroim2::NeuroSpace(c(20,20,20)))
#' bvec <- replicate(10,
#'                   neuroim2::NeuroVol(array(runif(20*20*20), c(20,20,20)),
#'                            neuroim2::NeuroSpace(c(20,20,20))),
#'                   simplify=FALSE)
#' bvec <- do.call(neuroim2::concat, bvec)
#' sbvec <- knn_shrink(bvec, mask, k=3)
#' }
#'
#' @export
knn_shrink <- function(bvec, mask, k = 5, connectivity = 27) {
  assert_that(inherits(bvec, "NeuroVec"))
  assert_that(inherits(mask, "NeuroVol"))
  assert_that(k >= 1)
  assert_that(connectivity >= k)

  mask <- as.logical(mask)
  mask.idx <- which(mask)
  coords <- index_to_coord(mask, mask.idx)
  feature_mat <- series(bvec, mask.idx)

  # find k-nearest neighbors within 'connectivity' radius
  neib <- FNN::get.knn(coords, k = connectivity)

  # for each voxel i, average the voxel and (k-1) neighbors
  sfeature_mat <- t(do.call(rbind, purrr::map(seq_len(nrow(neib$nn.index)), function(i) {
    rowMeans(feature_mat[, c(i, neib$nn.index[i, 1:(k-1)])])
  })))

  SparseNeuroVec(sfeature_mat, space(bvec), mask = mask)
}


#' Supervoxel Clustering (3D volumes)
#'
#' Cluster a \code{NeuroVec} instance into a set of spatially constrained clusters.
#'
#' @note Consider using \code{\link{cluster4d}} with \code{method = "supervoxels"} for a
#' standardized interface across all clustering methods.
#'
#' @param bvec A \code{\link[neuroim2:NeuroVec-class]{NeuroVec}} instance supplying the data to cluster.
#'   Can also be a 3D \code{\link[neuroim2:NeuroVol-class]{NeuroVol}} for structural image segmentation,
#'   which will be automatically converted to a single-timepoint NeuroVec internally.
#' @param mask A \code{\link[neuroim2:NeuroVol-class]{NeuroVol}} mask defining the voxels to include. If numeric, nonzero = included.
#' @param K The number of clusters to find (default 500).
#' @param sigma1 The bandwidth of the heat kernel for the data vectors.
#' @param sigma2 The bandwidth of the heat kernel for the coordinate vectors.
#' @param iterations The maximum number of cluster iterations.
#' @param connectivity The number of nearest neighbors defining the neighborhood.
#' @param use_medoid Logical; whether to use medoids rather than means for cluster centers.
#' @param use_gradient Logical; use the image gradient to initialize clusters if possible.
#' @param alpha The relative weighting of data similarity vs spatial similarity;
#'   \code{alpha=1} = all data weighting, \code{alpha=0} = purely spatial weighting.
#' @param parallel Logical; whether to use parallel processing for cluster assignment updates.
#'   Default is TRUE. Parallel processing is automatically disabled for small datasets (<1000 voxels).
#' @param grain_size Integer; the minimum number of voxels to process per parallel task.
#'   Default is 100. Smaller values provide better load balancing but increase overhead.
#' @param verbose Logical; whether to print detailed progress messages including convergence
#'   metrics. Default is FALSE.
#' @param converge_thresh Numeric; convergence threshold as proportion of voxels switching
#'   clusters. Algorithm stops when switch ratio falls below this value. Default is 0.001
#'   (0.1% of voxels).
#'
#' @return A \code{list} (of class \code{cluster_result}) with elements:
#'   \item{clusvol}{\code{ClusteredNeuroVol} containing the final clustering.}
#'   \item{cluster}{Integer vector of cluster assignments for each voxel.}
#'   \item{centers}{K-by-T matrix of final cluster means in the original feature space.}
#'   \item{coord_centers}{K-by-3 matrix of final spatial means in physical units.}
#'   \item{metadata}{Algorithm diagnostics, including whether parallel assignment
#'     was requested and used and any deterministic empty-label reseeds.}
#'
#' @importFrom neuroim2 series
#' @export
#'
#' @examples
#' \dontrun{
#' mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
#' bvec <- replicate(10,
#'                   NeuroVol(array(runif(20*20*20), c(20,20,20)),
#'                            NeuroSpace(c(20,20,20))),
#'                   simplify=FALSE)
#' bvec <- do.call(concat, bvec)
#' cres1 <- supervoxels(bvec, mask, K=100, sigma1=1, sigma2=3)
#' }
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Scale input data (\code{bvec}) so each feature dimension is centered and scaled.
#'   \item If \code{use_gradient = TRUE}, initialize cluster seeds using gradient-based heuristics.
#'   \item Run an iterative, spatially-constrained clustering that updates voxel assignments based on
#'         both feature similarity (bandwidth \code{sigma1}) and spatial proximity (bandwidth \code{sigma2}),
#'         weighted by \code{alpha}.
#'   \item Return the final clusters, plus the feature-space and coordinate-space centers.
#' }
#'
#' @details
#' ## Exact-K and parallel contracts
#'
#' Every assignment center must have a positive membership count. If an update
#' empties one or more labels, distinct voxels are deterministically moved from
#' non-singleton clusters in descending combined feature/spatial loss order.
#' Empty all-zero centers are therefore never passed to an assignment kernel.
#' Since `K <= N` is required, this reseed policy is feasible and the returned
#' labels, center rows, and requested K agree on every convergence path.
#'
#' With `parallel = TRUE`, voxel assignment uses RcppParallel when there are
#' more than 1000 included voxels; smaller inputs use the same sequential
#' scoring kernel. Parallel and sequential kernels implement the same scoring
#' and deterministic tie rule. Centroid reduction may also run in parallel for
#' sufficiently large K. `metadata$algorithm$parallel_used` reports whether the
#' parallel assignment path actually ran. No architecture-specific fallback is
#' applied.
#'
#' @param num_threads Optional integer to override the number of threads used by
#'   RcppParallel (defaults to package/global setting). Ignored if `parallel = FALSE`.
supervoxels <- function(bvec, mask,
                        K = 500,
                        sigma1 = NULL,
                        sigma2 = 2.5,
                        iterations = 50,
                        connectivity = 27,
                        use_medoid = FALSE,
                        use_gradient = TRUE,
                        alpha = 0.5,
                        parallel = TRUE,
                        grain_size = 100,
                        num_threads = NULL,
                        verbose = FALSE,
                        converge_thresh = 0.001) {

  # Support 3D NeuroVol input by converting to single-timepoint NeuroVec

  # This allows supervoxel segmentation of structural images (e.g., T1-weighted MRI)
  if (inherits(bvec, "NeuroVol") && length(dim(bvec)) == 3) {
    if (verbose) message("supervoxels: converting 3D NeuroVol to single-timepoint NeuroVec")
    arr4d <- array(as.array(bvec), dim = c(dim(bvec), 1))
    bvec <- neuroim2::NeuroVec(arr4d, neuroim2::add_dim(neuroim2::space(bvec), 1))
  }

  mask.idx <- which(mask > 0)
  if (length(mask.idx) == 0) {
    stop("No nonzero voxels in mask.")
  }

  # Thread control for RcppParallel
  if (parallel && !is.null(num_threads)) {
    old_opts <- RcppParallel::setThreadOptions(numThreads = num_threads)
    on.exit({
      tryCatch({
        if (is.list(old_opts) && !is.null(old_opts$numThreads)) {
          if (!is.null(old_opts$stackSize)) {
            RcppParallel::setThreadOptions(numThreads = old_opts$numThreads, stackSize = old_opts$stackSize)
          } else {
            RcppParallel::setThreadOptions(numThreads = old_opts$numThreads)
          }
        } else if (is.null(old_opts)) {
          tryCatch(RcppParallel::setThreadOptions(numThreads = "auto"), error = function(e) NULL)
        } else if (is.character(old_opts)) {
          tryCatch(RcppParallel::setThreadOptions(numThreads = old_opts), error = function(e) NULL)
        } else {
          RcppParallel::setThreadOptions(numThreads = as.integer(old_opts)[1])
        }
      }, error = function(e) NULL)
    }, add = TRUE)
    if (verbose) message("supervoxels: using ", num_threads, " threads (RcppParallel)")
  }
  
  # Early check for K vs number of voxels
  if (K > length(mask.idx)) {
    stop(sprintf("Cannot create %d clusters from %d masked voxels. K must be <= number of masked voxels.", 
                 K, length(mask.idx)))
  }

  # coordinate grid in mm units
  coords <- index_to_coord(mask, mask.idx)

  # gather scaled features (time x voxels); keep both orientations
  # NOTE: neuroim2::series returns a vector (not matrix) for single-timepoint data,
  # so we ensure it's always a T x N matrix
  feature_mat <- neuroim2::series(bvec, mask.idx)        # T x N (or vector if T=1)
  if (!is.matrix(feature_mat)) {
    # Single timepoint: convert vector to 1 x N matrix
    feature_mat <- matrix(feature_mat, nrow = 1)
  }
  feature_mat_vox <- t(as.matrix(feature_mat))           # N x T for downstream summaries

  # sigma1 will be set adaptively in supervoxel_cluster_fit if NULL

  # FIX 6: Let supervoxel_cluster_fit handle initialization with functional gradient seeding
  # Previously, init_cluster was called here which used spatial gradient from a single volume.
  # Now we pass initclus=NULL to trigger the new functional gradient seeding in supervoxel_cluster_fit,
  # which considers the full feature space for better seed placement.

  # run the iterative fit
  ret <- supervoxel_cluster_fit(feature_mat, coords, K = K, sigma1 = sigma1, sigma2 = sigma2,
                                iterations = iterations, connectivity = connectivity,
                                use_medoid = use_medoid, alpha = alpha,
                                initclus = NULL, use_gradient = use_gradient,
                                parallel = parallel, grain_size = grain_size,
                                verbose = verbose, converge_thresh = converge_thresh)

  # build the final ClusteredNeuroVol with consistent logical mask (only positive values are TRUE)
  logical_mask <- mask > 0
  kvol <- ClusteredNeuroVol(logical_mask, clusters = ret$clusters)

  # Prepare data for standardized result
  data_prep <- list(
    features = feature_mat_vox,
    coords = coords,
    mask_idx = mask.idx,
    n_voxels = length(mask.idx),
    dims = dim(mask),
    spacing = spacing(mask)
  )
  
  # Create standardized result
  result <- create_cluster4d_result(
    labels = ret$clusters,
    mask = mask,
    data_prep = data_prep,
    method = "supervoxels",
    parameters = list(
      K = K,
      sigma1 = sigma1,
      sigma2 = sigma2,
      iterations = iterations,
      connectivity = connectivity,
      use_medoid = use_medoid,
      use_gradient = use_gradient,
      alpha = alpha,
      parallel = parallel,
      grain_size = grain_size,
      converge_thresh = converge_thresh
    ),
    metadata = list(
      algorithm = list(
        exact_k_policy = "deterministic_distinct_farthest_loss_reseed",
        parallel_requested = parallel,
        parallel_used = ret$parallel_used,
        reseed_events = ret$reseed_events
      )
    ),
    # Public centers are means in the original input space, recomputed from the
    # final repaired labels rather than exposing scaled working centers.
    compute_centers = TRUE
  )
  
  # Ensure backward compatibility with old class
  class(result) <- c("cluster_result", "list")
  result
}


#' Supervoxel Clustering in Time
#'
#' Cluster feature matrix (rows = time points) in a "supervoxel" style but over temporal dimension.
#'
#' @param feature_mat A matrix (nrows = time points, ncols = features) or vice versa.
#' @param K Number of clusters.
#' @param sigma1 Heat kernel bandwidth for feature similarity (data vectors).
#' @param sigma2 Heat kernel bandwidth for spatial similarity (coordinate vectors).
#' @param iterations Maximum number of cluster iterations.
#' @param TR Repetition time (seconds).
#' @param filter List specifying optional frequency filters, e.g., \code{list(lp=0.1, hp=0)}.
#' @param use_medoid Whether to use medoids for cluster centers.
#' @param nreps Number of repeated initializations.
#'
#' @return A list of cluster results (one per repetition), each of which
#'   has the same structure as \code{supervoxel_cluster_fit()}.
#'
#' @export
#' @examples
#' \dontrun{
#' feature_mat <- matrix(rnorm(100 * 10), 100, 10)
#' library(future)
#' plan(multicore)
#' cres <- supervoxel_cluster_time(t(feature_mat), K=5)
#' }
supervoxel_cluster_time <- function(feature_mat,
                                    K = min(nrow(feature_mat), 100),
                                    sigma1 = 1,
                                    sigma2 = 3,
                                    iterations = 50,
                                    TR = 2,
                                    filter = list(lp = 0, hp = 0),
                                    use_medoid = FALSE,
                                    nreps = 5) {

  # optional filtering (filter_mat is in experimental/; warn if unavailable)
  if (filter$lp > 0 || filter$hp > 0) {
    if (exists("filter_mat", mode = "function")) {
      message("supervoxel_cluster_time: filtering time series")
      feature_mat <- get("filter_mat", mode = "function")(feature_mat, filter$lp, filter$hp)
    } else {
      warning("filter_mat not available; skipping time-series filtering")
    }
  }

  nels <- nrow(feature_mat)
  if (nels == 0) stop("No rows in feature_mat.")

  # create a 'coordinate' vector for time (0, TR, 2*TR, ...)
  coords <- matrix(seq(0, by = TR, length.out = nels), ncol = 1)

  # run multiple starts in parallel (requires 'future' and 'furrr')
  fits <- furrr::future_map(seq_len(nreps), function(i) {
    initsamps <- sort(sample(seq_len(nels), K))
    initcenters <- feature_mat[initsamps, , drop = FALSE]

    # quick assignment by nearest "time" coordinate
    curclus <- FNN::get.knnx(coords[initsamps, , drop = FALSE], coords, k = 1)$nn.index[,1]

    ret <- supervoxel_cluster_fit(
      t(feature_mat), coords,
      sigma1 = sigma1,
      sigma2 = sigma2,
      K = K,
      initclus = curclus,
      iterations = iterations,
      connectivity = 3,
      use_medoid = use_medoid
    )

    class(ret) <- c("cluster_result_time", "cluster_result", "list")
    ret
  })

  fits
}


#' Supervoxel Clustering on a Surface
#'
#' Cluster feature data on a cortical surface or mesh using a supervoxel-like approach.
#'
#' @param bsurf A \code{NeuroSurface} or similar object with geometry, coords, and data.
#' @param K Number of clusters.
#' @param sigma1 Heat kernel bandwidth for feature similarity (data vectors).
#' @param sigma2 Heat kernel bandwidth for spatial similarity (coordinate vectors).
#' @param iterations Max iterations.
#' @param connectivity Neighborhood size on the surface (e.g., # of nearest mesh neighbors).
#' @param use_medoid Whether to use medoids for cluster centers.
#'
#' @return A \code{list} with:
#' \describe{
#'   \item{clusvol}{A \code{NeuroSurface} storing the final clustering result.}
#'   \item{clusters}{Integer vector of cluster assignments (one per vertex).}
#'   \item{centers}{Matrix of cluster centers.}
#'   \item{coord_centers}{Matrix of spatial centroid coordinates.}
#'   \item{index_sets}{List of vertex indices for each cluster.}
#' }
#'
#' @export
supervoxel_cluster_surface <- function(bsurf,
                                       K = 500,
                                       sigma1 = 1,
                                       sigma2 = 5,
                                       iterations = 50,
                                       connectivity = 6,
                                       use_medoid = FALSE) {

  if (!requireNamespace("neurosurf", quietly = TRUE)) {
    stop("supervoxel_cluster_surface() requires the 'neurosurf' package. Please install it to use surface clustering.")
  }

  mask.idx <- neurosurf::indices(bsurf)
  coords <- neurosurf::coords(bsurf)[mask.idx, , drop = FALSE]
  feature_mat <- neuroim2::series(bsurf, mask.idx)

  ret <- supervoxel_cluster_fit(feature_mat, coords, K = K,
                                sigma1 = sigma1, sigma2 = sigma2,
                                iterations = iterations,
                                connectivity = connectivity,
                                use_medoid = use_medoid)

  # create a NeuroSurface storing clusters
  kvol <- neurosurf::NeuroSurface(
    geometry = neurosurf::geometry(bsurf),
    indices = mask.idx,
    data = ret$clusters
  )

  # also provide the sets of vertex indices for each cluster
  index_sets <- split(mask.idx, ret$clusters)

  out <- list(
    clusvol = kvol,
    clusters = ret$clusters,
    centers = ret$centers,
    coord_centers = ret$coord_centers,
    index_sets = index_sets
  )
  class(out) <- c("cluster_result_surface", "cluster_result", "list")
  out
}
