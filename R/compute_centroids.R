#' Compute Centroids of Clusters
#'
#' The compute_centroids function calculates the center and centroid of each cluster
#' given a feature matrix, grid, and cluster assignments.
#'
#' @param feature_mat A matrix of features, where columns represent data points
#'                    and rows represent features.
#' @param grid A matrix representing the spatial grid of data points.
#' @param assignment A vector containing the cluster assignment for each data point.
#' @param medoid A logical value indicating whether to calculate medoids instead
#'               of means for cluster centers and centroids. Default is FALSE.
#'
#' @return A list containing two elements:
#'         \item{center}{A matrix containing the centers of each cluster.}
#'         \item{centroid}{A matrix containing the centroids of each cluster.}
#'
#' @examples
#' # Simple synthetic example
#' feature_mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' grid <- cbind(x = runif(10), y = runif(10), z = runif(10))
#' assignment <- rep(1:2, each = 5)
#' centroids <- compute_centroids(feature_mat, grid, assignment)
#' print(names(centroids))
#' 
#' \dontrun{
#'   # Larger example with real neuroimaging data
#'   # Assuming `feature_mat`, `grid`, and `assignment` are available
#'   centroids <- compute_centroids(feature_mat, grid, assignment)
#'   # To compute medoids instead of means
#'   medoids <- compute_centroids(feature_mat, grid, assignment, medoid=TRUE)
#' }
#'
#' @export
compute_centroids <- function(feature_mat, grid, assignment, medoid=FALSE) {

  csplit <- split(1:length(assignment), assignment)

  if (!medoid) {
    purrr::transpose(purrr::map(csplit, function(id) {
      id <- id[id >= 1 & id <= ncol(feature_mat)]
      if (length(id) == 0) {
        return(list(center = rep(NA_real_, nrow(feature_mat)),
                    centroid = rep(NA_real_, ncol(grid))))
      }
      mat <- feature_mat[, id, drop=FALSE]
      coords <- grid[id, , drop=FALSE]
      list(center=rowMeans(mat), centroid=colMeans(coords))
    }))
  } else {
    purrr::transpose(purrr::map(csplit, function(id) {
      id <- id[id >= 1 & id <= ncol(feature_mat)]
      if (length(id) == 0) {
        return(list(center = rep(NA_real_, nrow(feature_mat)),
                    centroid = rep(NA_real_, ncol(grid))))
      }
      mat <- feature_mat[, id, drop=FALSE]
      coords <- grid[id,,drop=FALSE]
      coords_dist <- as.matrix(dist(coords))
      coords_medoid_ind <- which.min(rowSums(coords_dist))
      Dmat <- 1-cor(mat)
      mat_medoid_ind <- which.min(rowSums(Dmat))
      list(center=mat[,mat_medoid_ind], centroid=coords[coords_medoid_ind,])
    }))

  }
}
