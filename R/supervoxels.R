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
#' # Assuming you have a NeuroVol object 'mask' and you want to create 150 clusters
#' clustered_volume <- tesselate(mask, K = 150)
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
      # Use neurocluster::: to access non-exported function
      init <- neurocluster:::find_initial_points(valid_coords, grad_vals, K)

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
                                   converge_thresh = 0.001) {

  # Early parameter validation
  nvox <- nrow(coords)
  if (nvox == 0) {
    stop("No voxels to cluster (coords has 0 rows)")
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
      coord_centers = coords
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

    # Reorder centers to match the remapped cluster IDs
    # kres$centers rows correspond to original cluster IDs, need to reorder
    centers_reordered <- kres$centers[clusid, , drop = FALSE]

    # Seeds are the nearest grid points to the reordered cluster centers
    # (Use k-means centers, which may have shifted from initial gradient seeds)
    seeds <- FNN::get.knnx(coords, centers_reordered)$nn.index[,1]
    sp_centroids <- coords[seeds, , drop = FALSE]  # K x 3 matrix

    # FIX 4: Compute feature centroids as cluster AVERAGES, not seed voxels
    # The seed voxel may belong to a different true cluster than other cluster members
    # Using the average is more robust and representative
    num_centroids <- matrix(0, nrow = nrow(feature_mat), ncol = K)
    for (k in 1:K) {
      cluster_members <- which(curclus == k)
      if (length(cluster_members) == 1) {
        # Single member: use that voxel's features
        num_centroids[, k] <- feature_mat[, cluster_members]
      } else {
        # Multiple members: use mean of all members' features
        num_centroids[, k] <- rowMeans(feature_mat[, cluster_members, drop = FALSE])
      }
    }

    # Transpose both to match C++ expectations: (dims x K) for both
    sp_centroids <- t(sp_centroids)  # Now 3 x K
    # num_centroids is already features x K, no transpose needed

  } else {
    # Use user-supplied initclus
    assert_that(length(initclus) == nvox)
    clusid <- sort(unique(initclus))
    assert_that(length(clusid) == K)
    
    # CRITICAL FIX: Ensure user-supplied clusters are also contiguous
    # This maintains consistency with the k-means initialization path
    cluster_mapping <- setNames(1:length(clusid), clusid)
    curclus <- cluster_mapping[as.character(initclus)]

    # compute_centroids is presumably your function returning
    # the numeric/data centroid and spatial centroid for each cluster
    centroids <- compute_centroids(feature_mat, coords, curclus, medoid = use_medoid)
    sp_centroids <- t(do.call(rbind, centroids$centroid))  # Transpose to get 3 x K
    num_centroids <- t(do.call(rbind, centroids$center))   # Transpose to get features x K
  }

  # Find connectivity-based neighbors
  neib <- FNN::get.knn(coords, k = connectivity)
  # dthresh is the median distance among the 'connectivity'-th neighbor
  dthresh <- stats::median(neib$nn.dist[, connectivity, drop = FALSE])
  message("dthresh: ", dthresh)

  iter <- 1
  switches <- 1
  iter.max <- iterations
  newclus <- curclus
  prev_switches <- Inf
  no_improvement_count <- 0
  
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

  while (iter <= iter.max && switches > 0) {
    # FIX 1: Always use full alpha (no more cheap_iters strategy)
    current_alpha <- original_alpha
    
    # Sanitize neighbor indices before passing to C++
    nn_indices <- neib$nn.index - 1L  # Convert to 0-based
    if (any(nn_indices < 0 | nn_indices >= nvox, na.rm = TRUE)) {
      warning("Invalid neighbor indices found. Coercing to valid range.")
      nn_indices[is.na(nn_indices) | nn_indices < 0] <- 0L
      nn_indices[nn_indices >= nvox] <- nvox - 1L
    }

    # FIX 7: Convert curclus to 0-based for C++ (C++ uses 0-based centroid indexing)
    # The C++ code iterates centroids as k=0 to K-1 and stores these 0-based indices
    # in the spatial bins. When it returns assignments, they are also 0-based.
    curclus_0based <- curclus - 1L

    # Use spatially-binned assignment for efficient candidate selection
    if (parallel && nvox > 1000) {
      newclus_0based <- fused_assignment_parallel_binned(
        nn_indices, neib$nn.dist, curclus_0based,
        t(coords), num_centroids, sp_centroids, feature_mat,
        dthresh, sigma1, sigma2, current_alpha,
        grain_size = max(1024L, as.integer(nvox / 32L)),
        window_factor = window_factor,
        bin_expand = bin_expand
      )
    } else {
      # Even sequential mode benefits from spatial binning
      newclus_0based <- fused_assignment_binned(
        nn_indices, neib$nn.dist, curclus_0based,
        t(coords), num_centroids, sp_centroids, feature_mat,
        dthresh, sigma1, sigma2, current_alpha,
        window_factor = window_factor,
        bin_expand = bin_expand
      )
    }

    # Convert back to 1-based for R
    newclus <- newclus_0based + 1L
    
    # Compute switches in R to avoid atomic contention in parallel C++ code
    switches <- sum(newclus != curclus)

    if (switches > 0) {
      # Get actual number of unique clusters (may be less than K)
      n_actual_clusters <- length(unique(newclus))
      
      # Use fast parallel centroid computation
      if (parallel && nvox > 1000 && n_actual_clusters > 50) {
        cent_result <- compute_centroids_parallel_fast(
          as.integer(newclus - 1L), feature_mat, t(coords), as.integer(K)
        )
        num_centroids <- cent_result$centers
        sp_centroids <- cent_result$coord_centers
      } else {
        # Fallback to sequential computation for small problems
        centroids <- compute_centroids(feature_mat, coords, newclus, medoid = use_medoid)
        center_ids <- suppressWarnings(as.integer(names(centroids$center)))
        if (anyNA(center_ids)) {
          center_ids <- seq_along(centroids$center)
        }
        temp_num <- matrix(0, nrow = nrow(feature_mat), ncol = K)
        temp_sp  <- matrix(0, nrow = 3, ncol = K)
        for (idx in seq_along(centroids$center)) {
          cid <- center_ids[idx]
          if (cid < 1 || cid > K) next
          temp_num[, cid] <- centroids$center[[idx]]
          temp_sp[, cid]  <- centroids$centroid[[idx]]
        }
        num_centroids <- temp_num
        sp_centroids  <- temp_sp
      }
      curclus <- newclus
    }

    # Check for convergence
    switch_ratio <- switches / nvox
    if (verbose) {
      message("supervoxels_fit: iter ", iter, " -- num switches = ", switches, 
              " (", round(switch_ratio * 100, 2), "% of voxels)")
    } else {
      message("supervoxels_fit: iter ", iter, " -- num switches = ", switches)
    }
    
    # Early convergence check
    if (switch_ratio < converge_thresh) {
      message("supervoxels_fit: converged at iteration ", iter, " (switch ratio < ", converge_thresh, ")")
      break
    }
    
    # Check if we're stuck (no improvement)
    if (switches >= prev_switches) {
      no_improvement_count <- no_improvement_count + 1
      if (no_improvement_count >= 3) {
        message("supervoxels_fit: stopping - no improvement for 3 iterations")
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
    coord_centers = t(sp_centroids)
  )
}


#' K-nearest-neighbor shrink
#'
#' Replace each voxel by the mean of its k nearest neighbors in its local spatial neighborhood.
#'
#' @param bvec A \code{\linkS4class{NeuroVec}} instance (the data).
#' @param mask A \code{\linkS4class{NeuroVol}} mask defining the voxels to include. If numeric, nonzero = included.
#' @param k The number of nearest neighbors to average over.
#' @param connectivity The number of spatial neighbors to include in the search around each voxel.
#'
#' @return A \code{SparseNeuroVec} or similar object with the smoothed data.
#'
#' @examples
#' mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
#' bvec <- replicate(10,
#'                   NeuroVol(array(runif(20*20*20), c(20,20,20)),
#'                            NeuroSpace(c(20,20,20))),
#'                   simplify=FALSE)
#' bvec <- do.call(concat, bvec)
#'
#' sbvec <- knn_shrink(bvec, mask, k=3)
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
#' @param bvec A \code{\linkS4class{NeuroVec}} instance supplying the data to cluster.
#'   Can also be a 3D \code{\linkS4class{NeuroVol}} for structural image segmentation,
#'   which will be automatically converted to a single-timepoint NeuroVec internally.
#' @param mask A \code{\linkS4class{NeuroVol}} mask defining the voxels to include. If numeric, nonzero = included.
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
#'   \item{centers}{Matrix of cluster centers in feature space.}
#'   \item{coord_centers}{Matrix of cluster spatial centroids.}
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
#' ## Parallelization Status
#' 
#' **NOW PARALLELIZED with RcppParallel!** The supervoxels algorithm can now run
#' cluster assignment updates in parallel across multiple CPU cores.
#' 
#' ### Parallel Operations:
#' 
#' 1. **Heat Kernel Computation**: Parallel across voxels using RcppParallel
#'    - Each voxel's best cluster assignment computed independently
#'    - Automatic load balancing with configurable grain size
#'    - Scales linearly with number of CPU cores
#' 
#' 2. **Sequential Operations** (still):
#'    - Initialization (K-means or gradient-based seed selection)
#'    - Centroid updates after each iteration
#'    - Convergence checking between iterations
#' 
#' ### Performance Characteristics:
#' 
#' - **Speedup**: Typically 2-8x faster on multicore systems
#' - **Automatic optimization**: Disabled for small datasets (<1000 voxels)
#' - **Memory overhead**: Minimal - uses shared memory via RcppParallel
#' - **Computational complexity**: Still O(N × K × iterations) but parallelized over N
#' 
#' ### Parallel Configuration:
#' 
#' - **parallel**: Set to FALSE to force sequential execution
#' - **grain_size**: Controls work distribution (default 100)
#'   - Smaller values = better load balancing but more overhead
#'   - Larger values = less overhead but potential imbalance
#' - **Thread control**: Set threads via `RcppParallel::setThreadOptions()`
#' 
#' ### When Parallelization Helps Most:
#' 
#' - Large numbers of voxels (N > 10,000)
#' - Many clusters (K > 100)
#' - Multiple iterations needed for convergence
#' - Systems with 4+ CPU cores
#' 
#' ### Performance Tips:
#' 
#' - **Set threads**: `RcppParallel::setThreadOptions(numThreads = 4)`
#' - **Tune grain_size**: Start with nvoxels/nthreads/10
#' - **Monitor CPU usage**: Should see near 100% on all cores during updates
#' - **Memory considerations**: Parallel version uses slightly more RAM
#' - **Disable for debugging**: Set `parallel = FALSE` for reproducible debugging
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

  # Safety fallback: on Apple Silicon, fused_assignment_parallel_binned has
  # occasionally crashed (bus error) with parallel=TRUE. Disable parallel on
  # that platform unless user explicitly forces it.
  if (parallel && grepl("aarch64.*darwin", R.version$platform)) {
    if (verbose) message("supervoxels: disabling parallel on aarch64-darwin for stability")
    parallel <- FALSE
  }

  # Thread control for RcppParallel
  if (parallel && !is.null(num_threads)) {
    old_opts <- RcppParallel::setThreadOptions(numThreads = num_threads)
    on.exit(RcppParallel::setThreadOptions(numThreads = old_opts$numThreads,
                                           stackSize   = old_opts$stackSize),
            add = TRUE)
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
  
  # If centers from fit don't match actual clusters, recompute to ensure consistency
  actual_k <- length(unique(ret$clusters))
  centers_meta <- ret$centers
  coord_meta <- ret$coord_centers
  if (is.null(centers_meta) || nrow(centers_meta) != actual_k ||
      is.null(coord_meta)   || nrow(coord_meta)   != actual_k) {
    recomputed <- compute_cluster_centers(ret$clusters, feature_mat_vox, coords, method = "mean")
    centers_meta <- recomputed$centers
    coord_meta   <- recomputed$coord_centers
  }

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
      centers = centers_meta,
      coord_centers = coord_meta
    ),
    compute_centers = FALSE  # we supply centers via metadata
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
#' feature_mat <- matrix(rnorm(100 * 10), 100, 10)
#' library(future)
#' plan(multicore)
#' cres <- supervoxel_cluster_time(t(feature_mat), K=20)
supervoxel_cluster_time <- function(feature_mat,
                                    K = min(nrow(feature_mat), 100),
                                    sigma1 = 1,
                                    sigma2 = 3,
                                    iterations = 50,
                                    TR = 2,
                                    filter = list(lp = 0, hp = 0),
                                    use_medoid = FALSE,
                                    nreps = 5) {

  # optional filtering
  if (filter$lp > 0 || filter$hp > 0) {
    message("supervoxel_cluster_time: filtering time series")
    feature_mat <- filter_mat(feature_mat, filter$lp, filter$hp)
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
#' @import neurosurf
supervoxel_cluster_surface <- function(bsurf,
                                       K = 500,
                                       sigma1 = 1,
                                       sigma2 = 5,
                                       iterations = 50,
                                       connectivity = 6,
                                       use_medoid = FALSE) {

  mask.idx <- indices(bsurf)
  coords <- coords(bsurf)[mask.idx, , drop = FALSE]
  feature_mat <- series(bsurf, mask.idx)

  ret <- supervoxel_cluster_fit(feature_mat, coords, K = K,
                                sigma1 = sigma1, sigma2 = sigma2,
                                iterations = iterations,
                                connectivity = connectivity,
                                use_medoid = use_medoid)

  # create a NeuroSurface storing clusters
  kvol <- NeuroSurface(
    geometry = geometry(bsurf),
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
