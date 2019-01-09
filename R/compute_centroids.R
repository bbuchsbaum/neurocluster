

#' compute_centroids
#'
#' Get centroids for a matrix and set of assignments
#'
#' @param feature_mat a matrix for features where each column is a feature and each row is an observation
#' @param grid a matrix of coordinates, where each row is a coordinates associated with the nth feature.
#' @param assignment a vector of integer indices indicating the cluster assignment
#' @param medoid compute the medoid rather than the mean of the cluster
#' @keywords internal
#' @importFrom purrr map
#' @export
compute_centroids <- function(feature_mat, grid, assignment, medoid=FALSE) {

  csplit <- split(1:length(assignment), assignment)

  if (!medoid) {
    purrr::transpose(purrr::map(csplit, function(id) {
      mat <- feature_mat[, id]
      coords <- grid[id,,drop=FALSE]
      list(center=rowMeans(mat), centroid=colMeans(coords))
    }))
  } else {
    purrr::transpose(purrr::map(csplit, function(id) {
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
