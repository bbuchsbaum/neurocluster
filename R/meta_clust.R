#' Meta Clustering for Cluster Results
#'
#' The meta_clust function performs meta clustering on a given clustering result
#' by applying hierarchical clustering or other clustering algorithms.
#'
#' @param x A clustering result, typically an object of class \code{"cluster_result"}.
#' @param cuts The number of cluster cuts to consider. Default is the minimum
#'             of half the number of centers and 2.
#' @param ... Additional arguments:
#'   \describe{
#'     \item{algo}{A character string indicating the clustering algorithm to use.
#'                 Default is "hclust" (hierarchical clustering).}
#'     \item{hclust_method}{A character string specifying the agglomeration method
#'                          to use for hierarchical clustering. Default is "ward.D".}
#'   }
#'
#' @return A list containing:
#'         \item{cvols}{A list of \code{\linkS4class{ClusteredNeuroVol}} instances.}
#'         \item{cuts}{The number of cluster cuts.}
#'         \item{cutmat}{A matrix representing the cluster assignments for each cut.}
#'         \item{hclus}{The hierarchical clustering result.}
#'
#' @details
#' ## Parallelization Status
#' 
#' **Currently NOT parallelized.** Meta-clustering performs hierarchical clustering
#' on cluster centers using sequential R functions.
#' 
#' ### Sequential Operations:
#' 
#' 1. **Distance Matrix Computation**: 
#'    - Computes correlation distance (1 - cor) between all cluster centers
#'    - O(K²) pairwise correlations where K = number of input clusters
#'    - Sequential double loop in base R
#' 
#' 2. **Hierarchical Clustering**: `hclust()`
#'    - Standard hierarchical agglomeration (Ward.D by default)
#'    - Sequential merging of closest clusters
#'    - O(K² log K) complexity
#' 
#' 3. **Tree Cutting**: `cutree()`
#'    - Cuts dendrogram at multiple heights
#'    - Creates nested cluster assignments
#'    - Linear in number of clusters
#' 
#' ### Why Not Parallelized:
#' 
#' - **Small scale**: Usually operates on hundreds of clusters, not thousands of voxels
#' - **R built-ins**: Uses standard `hclust()` which is sequential
#' - **Fast enough**: Typically completes in seconds even for large K
#' - **Memory efficient**: Only stores K×K distance matrix
#' 
#' ### Potential for Parallelization:
#' 
#' - Distance matrix computation could use parallel correlation
#' - Some hierarchical clustering variants support parallelization
#' - Multiple cuts could be computed in parallel
#' 
#' ### Performance Characteristics:
#' 
#' - **Input size**: K clusters from initial clustering (typically 100-1000)
#' - **Complexity**: O(K²) for distances, O(K² log K) for clustering
#' - **Memory**: O(K²) for distance matrix
#' - **Speed**: Usually < 1 second for K < 500
#' 
#' ### Performance Tips:
#' 
#' - **Pre-reduce clusters**: Start with fewer initial clusters if meta-clustering is slow
#' - **Use appropriate linkage**: Ward.D is slower but often gives better results
#' - **Limit cuts**: Fewer cut levels = faster processing
#' - **Alternative**: Use consensus clustering (`merge_clus()`) for different approach
#' 
#' ### Use Cases:
#' 
#' Meta-clustering is useful for:
#' - Creating multi-resolution parcellations
#' - Hierarchical organization of functional regions
#' - Reducing large numbers of clusters to interpretable groups
#' - Finding stable cluster boundaries across scales
#'
#' @seealso \code{\link[stats]{hclust}}, \code{\link[stats]{cutree}}
#' @export
meta_clust.cluster_result <- function(x, cuts=min(as.integer(length(x$centers)/2),2), ...) {
  # Extract additional arguments
  dots <- list(...)
  algo <- if (!is.null(dots$algo)) dots$algo else "hclust"
  hclust_method <- if (!is.null(dots$hclust_method)) dots$hclust_method else "ward.D"

  orig <- x$cluster

  # Check if centers exist (some methods like SNIC don't provide centers)
  if (is.null(x$centers)) {
    stop("Cannot perform meta clustering: the clustering result does not contain cluster centers. ",
         "Methods like SNIC do not compute explicit cluster centers.")
  }

  cvols <- if (algo == "hclust") {
    # x$centers is already a matrix, no need to rbind
    cen <- x$centers
    #D <- 1 - cor(t(x$centers))

    D <- 1 - cor(t(cen))
    hres <- hclust(as.dist(D), method=hclust_method)
    cmat <- do.call(cbind, lapply(cuts, function(i) {
      cind <- cutree(hres, i)
      cind[orig]
    }))

    cvols <- lapply(1:ncol(cmat), function(i) {
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
#' @importFrom assertthat assert_that
#' @importFrom purrr map_lgl
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

#' Extract Class IDs from Cluster Result
#'
#' This function extracts the cluster class identifiers from a cluster result object.
#' It is a method for the \code{\link[clue]{cl_class_ids}} generic from the \code{clue} package.
#'
#' @param x A \code{cluster_result} object containing clustering information.
#'
#' @return An integer vector of cluster assignments, one for each data point.
#'
#' @seealso \code{\link[clue]{cl_class_ids}} for the generic function.
#' @method cl_class_ids cluster_result
#' @export
cl_class_ids.cluster_result <- function(x) {
  x$cluster
}


#' Test if Object is a Partition
#'
#' This function tests whether a cluster result object represents a partition.
#' It is a method for the \code{\link[clue]{is.cl_partition}} generic from the \code{clue} package.
#' For \code{cluster_result} objects, this always returns \code{TRUE} since cluster results
#' represent valid partitions where each data point belongs to exactly one cluster.
#'
#' @param x A \code{cluster_result} object.
#'
#' @return \code{TRUE}, indicating that cluster results are always valid partitions.
#'
#' @seealso \code{\link[clue]{is.cl_partition}} for the generic function.
#' @method is.cl_partition cluster_result
#' @export
is.cl_partition.cluster_result <- function(x) {
  TRUE
}

