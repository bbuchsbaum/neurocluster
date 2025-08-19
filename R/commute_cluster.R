
#' @keywords internal
#' @noRd
commute_cluster_fit <- function(X, cds, K, ncomp=ceiling(sqrt(K*2)),
                                   alpha=.2, sigma1=.73, sigma2=5,
                                   connectivity=27,
                                   weight_mode=c("binary", "heat")) {

  weight_mode <- match.arg(weight_mode)

  csds <- matrixStats::colSds(X)
  bad <- which(is.na(csds) | csds==0)

  if (length(bad) > 0) {
    warning(paste(length(bad), "features have NA or 0 stdev. Random vectors"))
    X[, bad] <- matrix(rnorm(length(bad)*nrow(X)), nrow(X), length(bad))
  }

  X <- base::scale(X)
  W <- neighborweights::weighted_spatial_adjacency(cds, t(X),
                                                   dthresh=sigma2*2,
                                                   nnk=connectivity,
                                                   wsigma=sigma1,
                                                   sigma=sigma2,
                                                   alpha=alpha,
                                                   weight_mode=weight_mode,
                                                   include_diagonal=FALSE,
                                                   stochastic=TRUE)


  message("commute_cluster: computing commute-time embedding.")
  
  # Handle numerical stability issues for commute_cluster_fit
  ct <- tryCatch({
    neighborweights::commute_time_distance(W, ncomp=ncomp)
  }, error = function(e) {
    if (grepl("TridiagEigen", e$message)) {
      stop("Eigenvalue decomposition failed in commute clustering. ",
           "This often occurs with perfectly correlated signals or singular weight matrices. ",
           "Try: 1) Using fewer components (ncomp), 2) Adjusting alpha parameter, ",
           "3) Checking for duplicate time series, or 4) Adding small noise to your data. ",
           "Original error: ", e$message)
    } else {
      stop(e)
    }
  })

  message("commute_cluster: performing k-means with ", K, " clusters.")
  kres <- kmeans(ct$cds, center=K, iter.max=100)
}



#' Commute Time Clustering
#'
#' The commute_cluster function performs spatially constrained clustering on a \code{NeuroVec} instance
#' using the commute time distance and K-means clustering.
#'
#' @param bvec A \code{NeuroVec} instance supplying the data to cluster.
#' @param mask A \code{NeuroVol} mask defining the voxels to include in the clustering result.
#' If the mask contains \code{numeric} data, nonzero values will define the included voxels.
#' If the mask is a \code{\linkS4class{LogicalNeuroVol}}, then \code{TRUE} will define the set
#' of included voxels.
#' @param K The number of clusters to find. Default is 100.
#' @param ncomp The number of components to use for the commute time embedding. Default is the ceiling of \code{sqrt(K2)}.
#' @param alpha A numeric value controlling the balance between spatial and feature similarity. Default is 0.5.
#' @param sigma1 A numeric value controlling the spatial weighting function. Default is 0.73.
#' @param sigma2 A numeric value controlling the feature weighting function. Default is 5.
#' @param connectivity An integer representing the number of nearest neighbors to consider when constructing the similarity graph.
#' Default is 27.
#' @param weight_mode A character string indicating the type of weight function for the similarity graph. Options are "binary" and "heat".
#' Default is "heat".
#'
#' @return A \code{list} of class \code{commute_time_cluster_result} with the following elements:
#' \describe{
#' \item{clusvol}{An instance of type \linkS4class{ClusteredNeuroVol}.}
#' \item{cluster}{A vector of cluster indices equal to the number of voxels in the mask.}
#' \item{centers}{A matrix of cluster centers with each column representing the feature vector for a cluster.}
#' \item{coord_centers}{A matrix of spatial coordinates with each row corresponding to a cluster.}
#' }
#'
#' @examples
#' mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
#' vec <- replicate(10, NeuroVol(array(runif(202020), c(20,20,20)),
#' NeuroSpace(c(20,20,20))), simplify=FALSE)
#' vec <- do.call(concat, vec)
#'
#' commute_res <- commute_cluster(vec, mask, K=100)
#'
#' @seealso \code{\link{snic}}
#' @importFrom neuroim2 NeuroVec NeuroVol
#' @importFrom matrixStats colSds
#' @importFrom neighborweights weighted_spatial_adjacency commute_time_distance
#' @importFrom stats kmeans
#' @import assertthat
#'
#' @export
commute_cluster <- function(bvec, mask,
                            K=100,
                            ncomp=ceiling(sqrt(K*2)),
                            alpha=.5,
                            sigma1=.73,
                            sigma2=5,
                            connectivity=27,
                            weight_mode=c("binary", "heat")) {

  weight_mode <- match.arg(weight_mode)

  mask.idx <- which(mask > 0)
  grid <- index_to_coord(mask, mask.idx)

  feature_mat <- neuroim2::series(bvec, mask.idx)
  csds <- matrixStats::colSds(feature_mat)
  bad <- which(is.na(csds) | csds==0)

  if (length(bad) > 0) {
    warning(paste(length(bad), "have time-courses have NA or 0 stdev. Random vectors"))
    feature_mat[, bad] <- matrix(rnorm(length(bad)*nrow(feature_mat)), nrow(feature_mat), length(bad))
  }

  feature_mat <- base::scale(feature_mat)

  message("commute_cluster: computing similarity matrix.")

  W <- neighborweights::weighted_spatial_adjacency(grid, t(feature_mat),
                                                   dthresh=sigma2*6,
                                                   nnk=connectivity,
                                                   wsigma=sigma1,
                                                   sigma=sigma2,
                                                   alpha=alpha,
                                                   weight_mode=weight_mode,
                                                   include_diagonal=FALSE,
                                                   stochastic=TRUE)

  # W might be a sparse matrix from Matrix package
  if (inherits(W, "Matrix")) {
    W <- (W + Matrix::t(W))/2
  } else {
    W <- (W + t(W))/2
  }
  message("commute_cluster: computing commute-time embedding.")
  
  # Handle numerical stability issues for main commute_cluster function
  ct <- tryCatch({
    neighborweights::commute_time_distance(W, ncomp=ncomp)
  }, error = function(e) {
    if (grepl("TridiagEigen", e$message)) {
      stop("Eigenvalue decomposition failed in commute clustering. ",
           "This often occurs with perfectly correlated signals or singular weight matrices. ",
           "Try: 1) Using fewer components (ncomp), 2) Adjusting alpha parameter, ",
           "3) Checking for duplicate time series, or 4) Adding small noise to your data. ",
           "Original error: ", e$message)
    } else {
      stop(e)
    }
  })

  message("commute_cluster: performing k-means with ", K, " clusters.")
  kres <- kmeans(ct$cds, center=K, iter.max=500)
  # Create ClusteredNeuroVol with consistent logical mask (only positive values are TRUE)
  logical_mask <- mask > 0
  kvol <- ClusteredNeuroVol(logical_mask,kres$cluster)

  message("commute_cluster: computing final centroids")
  centroids <- compute_centroids(feature_mat, grid, kres$cluster, medoid=FALSE)

  ret <- structure(
    list(clusvol=kvol,
         cluster=kres$cluster,
         centers=centroids$center,
         coord_centers=centroids$centroid),
    class=c("commute_time_cluster_result", "cluster_result", "list"))

  ret
}


