#' Meta Clustering for Cluster Results
#'
#' The meta_clust function performs meta clustering on a given clustering result
#' by applying hierarchical clustering or other clustering algorithms.
#'
#' @param x A clustering result, typically an object of class \code{"cluster_result"}.
#' @param cuts The number of cluster cuts to consider. Default is the minimum
#'             of half the number of centers and 2.
#' @param algo A character string indicating the clustering algorithm to use.
#'             Default is "hclust" (hierarchical clustering).
#' @param hclust_method A character string specifying the agglomeration method
#'                      to use for hierarchical clustering. Default is "ward.D".
#'
#' @return A list containing:
#'         \item{cvols}{A list of \code{\linkS4class{ClusteredNeuroVol}} instances.}
#'         \item{cuts}{The number of cluster cuts.}
#'         \item{cutmat}{A matrix representing the cluster assignments for each cut.}
#'         \item{hclus}{The hierarchical clustering result.}
#'
#' @seealso \code{\link[hclust]{hclust}}, \code{\link[stats]{cutree}}
#' @export
meta_clust.cluster_result <- function(x, cuts=min(as.integer(length(x$centers)/2),2),
                                      algo="hclust", hclust_method="ward.D") {

  orig <- x$cluster

  cvols <- if (algo == "hclust") {
    cen <- do.call(rbind, x$centers)
    #D <- 1 - cor(t(x$centers))

    D <- 1 - cor(t(cen))
    hres <- hclust(as.dist(D), method=hclust_method)
    cmat <- do.call(cbind, lapply(cuts, function(i) {
      cind <- cutree(hres, i)
      cind[orig]
    }))

    cvols <- map(1:ncol(cmat), function(i) {
      ClusteredNeuroVol(x$clusvol@mask, cluster=cmat[,i])
    })

    c(cvols, list(x$clusvol))
  } else {
    stop()
  }
  # } else if (algo == "apcluster") {
  #   S <- cor(t(x$centers))
  #   apres <- apcluster(S, details=TRUE)
  #   out <- integer(nrow(x$centers))
  #   clmap <- do.call(rbind, imap(apres@clusters, function(cl,i) {
  #     cbind(cl,i)
  #   }))
  #
  #   cltab <- as.list(clmap[,2])
  #   names(cltab) <- as.character(clmap[,1])
  #   clusters <- unlist(cltab[as.character(orig)])
  # }

  list(cvols=cvols, cuts=cuts, cutmat=cmat, hclus=hres)

}

#' Merge Clustering Results for ClusteredNeuroVol Objects
#'
#' This method of merge_clus is specifically designed to merge clustering results for \code{ClusteredNeuroVol} objects.
#'
#' @param x A \code{ClusteredNeuroVol} object or an object of class \code{"cluster_result"}.
#' @param method A character string indicating the consensus clustering algorithm to use. Default is "SE".
#' See \code{\link[clue]{cl_consensus}} for available methods.
#' @param ... Additional clustering results to be merged.
#'
#' @return A \code{\linkS4class{ClusteredNeuroVol}} instance.
#'
#' @seealso \code{\link[clue]{cl_consensus}}, \code{\link[clue]{as.cl_hard_partition}}, \code{\link[clue]{cl_ensemble}}
#' @importFrom clue cl_consensus as.cl_hard_partition cl_ensemble
#' @importFrom assertthat assert_that map_lgl
#' @export
merge_clus.cluster_result <- function(x, method="SE", ...) {
  args <- c(list(x), list(...))
  assert_that(all(map_lgl(args, ~ inherits(., "cluster_result"))))

  ens <- do.call(clue::cl_ensemble, args)
  cons <- clue::cl_consensus(ens, method="SE")
  hpart <- clue::as.cl_hard_partition(cons)
  as.integer(hpart$.Data)
  #ClusteredNeuroVol(x$clusvol@mask, clusters=hpart$.Data)
}


#' @export
merge_clus.cluster_result_time <- function(x, ...) {
  args <- c(list(x), list(...))
  assert_that(all(map_lgl(args, ~ inherits(., "cluster_result"))))

  ens <- do.call(clue::cl_ensemble, args)
  cons <- clue::cl_consensus(ens)
  hpart <- clue::as.cl_hard_partition(cons)
  as.integer(hpart$.Data)
}

#' @export
cl_class_ids.cluster_result <- function(x) {
  x$cluster
}


#' @export
is.cl_partition.cluster_result <- function(x) {
  TRUE
}

