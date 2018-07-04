#' meta_clust
#'
#' hierarchical clustering of an existing cluster solution.
#'
#' @export
meta_hclust <- function(x, cuts) {
  UseMethod("mean", x)
}



#' merge a set of clustering results using a consensus clustering algorithm
#'
#' @export
merge_clus <- function(x, ...) {
  UseMethod("merge_clus", x)
}
