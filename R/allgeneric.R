#' Meta Clustering for Cluster Results (S4 Generic)
#'
#' The meta_clust generic function is used to define methods for different
#' classes of input objects.
#'
#' @param x A clustering result, typically an object of class \code{"cluster_result"}.
#' @param cuts The number of cluster cuts to consider.
#' @param ... Additional arguments passed to methods.
#'
#' @return Depends on the method called.
#'
#' @export
meta_clust <- function(x, cuts, ...) {
  UseMethod("meta_clust", x)
}



#' Merge Clustering Results Using a Consensus Clustering Algorithm
#'
#' The merge_clus function combines a set of clustering results using a consensus clustering algorithm.
#'
#' @param x A clustering result, typically a list or an object of class \code{"cluster_result"}.
#' @param method A character string indicating the consensus clustering algorithm to use. Default is "SE".
#' See \code{\link[clue]{cl_consensus}} for available methods.
#' @param ... Additional clustering results to be merged.
#'
#' @return A \code{\linkS4class{ClusteredNeuroVol}} instance.
#'
#' @examples
#' # Assuming clustering1, clustering2, and clustering3 are objects of class "cluster_result"
#' merged_clustering <- merge_clus(clustering1, clustering2, clustering3, method="SE")
#'
#' @seealso \code{\link[clue]{cl_consensus}}, \code{\link[clue]{as.cl_hard_partition}}, \code{\link[clue]{cl_ensemble}}
#' @importFrom clue cl_consensus as.cl_hard_partition cl_ensemble
#' @export
merge_clus <- function(x, method, ...) {
  UseMethod("merge_clus", x, method)
}
