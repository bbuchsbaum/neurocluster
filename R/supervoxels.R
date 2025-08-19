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
      valid_coords <- index_to_grid(mask, mask.idx)

      # find_initial_points function for gradient-based initialization
      # Use neurocluster::: to access non-exported function
      init <- neurocluster:::find_initial_points(valid_coords, grad, K)

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
                                   use_gradient = TRUE) {

  message("supervoxel_fit: sigma1 = ", sigma1, " sigma2 = ", sigma2)

  # Adjust connectivity based on number of voxels
  max_connectivity <- min(connectivity, nrow(coords) - 1)
  if (max_connectivity != connectivity) {
    warning(sprintf("Connectivity reduced from %d to %d due to limited voxels", connectivity, max_connectivity))
    connectivity <- max_connectivity
  }
  
  assert_that(connectivity > 1 & connectivity <= 27)
  assert_that(alpha >= 0 && alpha <= 1)

  # Center and scale the feature matrix
  feature_mat <- base::scale(feature_mat, center = TRUE, scale = TRUE)

  nvox <- nrow(coords)
  if (K >= nvox) {
    if (nvox == 1) {
      # Special case: only one voxel, return trivial result
      return(list(
        clusters = rep(1, nvox),
        centers = matrix(feature_mat[, 1], ncol = 1),
        coord_centers = matrix(coords[1, ], nrow = 1)
      ))
    }
    warning("K is greater than or equal to the number of data points. Setting K = ", nvox - 1)
    K <- nvox - 1
  }

  # If no initial clusters, do naive coordinate-based kmeans
  if (is.null(initclus)) {
    init_centers <- coords[as.integer(seq(1, nvox, length.out = K)), , drop = FALSE]
    kres <- stats::kmeans(coords, centers = init_centers, iter.max = 500)

    clusid <- sort(unique(kres$cluster))
    curclus <- kres$cluster

    # seeds are the nearest grid points to the cluster centers
    seeds <- FNN::get.knnx(coords, kres$centers)$nn.index[,1]
    sp_centroids <- coords[seeds, , drop = FALSE]

    # feature_mat is (#observations x #features) if transposed, so check usage carefully
    num_centroids <- feature_mat[, seeds, drop = FALSE]

  } else {
    # Use user-supplied initclus
    assert_that(length(initclus) == nvox)
    clusid <- sort(unique(initclus))
    assert_that(length(clusid) == K)
    curclus <- initclus

    # compute_centroids is presumably your function returning
    # the numeric/data centroid and spatial centroid for each cluster
    centroids <- compute_centroids(feature_mat, coords, curclus, medoid = use_medoid)
    sp_centroids <- do.call(rbind, centroids$centroid)
    num_centroids <- do.call(rbind, centroids$center)
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

  while (iter <= iter.max && switches > 0) {
    # find candidate clusters for each voxel using neighborhood distance threshold
    candlist <- find_candidates(neib$nn.index - 1, neib$nn.dist, curclus, dthresh)

    # choose best among the candidate clusters for each voxel
    newclus <- best_candidate(candlist, curclus, t(coords), t(num_centroids),
                              t(sp_centroids), feature_mat, sigma1, sigma2, alpha)
    switches <- attr(newclus, "nswitches")

    if (switches > 0) {
      # recompute centroids for updated assignment
      centroids <- compute_centroids(feature_mat, coords, newclus, medoid = use_medoid)
      sp_centroids <- do.call(rbind, centroids$centroid)
      num_centroids <- do.call(rbind, centroids$center)
      curclus <- newclus
    }

    message("supervoxels_fit: iter ", iter, " -- num switches = ", switches)
    iter <- iter + 1
  }

  # Return the final clusters and the final centroids
  list(
    clusters = newclus,
    centers = do.call(rbind, centroids$center),
    coord_centers = do.call(rbind, centroids$centroid)
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
#' @param bvec A \code{\linkS4class{NeuroVec}} instance supplying the data to cluster.
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
#' mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
#' bvec <- replicate(10,
#'                   NeuroVol(array(runif(20*20*20), c(20,20,20)),
#'                            NeuroSpace(c(20,20,20))),
#'                   simplify=FALSE)
#' bvec <- do.call(concat, bvec)
#' cres1 <- supervoxels(bvec, mask, K=100, sigma1=1, sigma2=3)
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
supervoxels <- function(bvec, mask,
                        K = 500,
                        sigma1 = 1,
                        sigma2 = 2.5,
                        iterations = 50,
                        connectivity = 27,
                        use_medoid = FALSE,
                        use_gradient = TRUE,
                        alpha = 0.5) {

  mask.idx <- which(mask > 0)
  if (length(mask.idx) == 0) {
    stop("No nonzero voxels in mask.")
  }

  # coordinate grid in mm units
  coords <- index_to_coord(mask, mask.idx)

  # initialization of cluster assignments
  clusinit <- init_cluster(bvec, mask, coords, K, use_gradient)

  # gather scaled features
  feature_mat <- neuroim2::series(bvec, mask.idx)

  # run the iterative fit
  ret <- supervoxel_cluster_fit(feature_mat, coords, K = K, sigma1 = sigma1, sigma2 = sigma2,
                                iterations = iterations, connectivity = connectivity,
                                use_medoid = use_medoid, alpha = alpha,
                                initclus = clusinit, use_gradient = use_gradient)

  # build the final ClusteredNeuroVol with consistent logical mask (only positive values are TRUE)
  logical_mask <- mask > 0
  kvol <- ClusteredNeuroVol(logical_mask, clusters = ret$clusters)

  # Return a structured list
  structure(
    list(
      clusvol = kvol,
      cluster = ret$clusters,
      centers = ret$centers,        # <-- corrected (was ret$center in original)
      coord_centers = ret$coord_centers
    ),
    class = c("cluster_result", "list")
  )
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


