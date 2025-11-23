#' Meta Clustering for Cluster Results
#'
#' The meta_clust function performs meta clustering on a given clustering result
#' by applying hierarchical clustering on the cluster centers.
#'
#' @param x A clustering result, typically an object of class \code{"cluster_result"}.
#' @param cuts Integer vector specifying the number of cluster cuts to consider.
#'             Default is \code{NULL}, which generates cuts at 2, 5, and 10 clusters
#'             (or fewer depending on the number of input clusters).
#' @param ... Additional arguments:
#'   \describe{
#'     \item{algo}{A character string indicating the clustering algorithm to use.
#'                 Default is "hclust" (hierarchical clustering).}
#'     \item{dist_method}{Character string: "correlation" (default) or "euclidean"
#'                        for distance calculation between cluster centers.}
#'     \item{hclust_method}{A character string specifying the agglomeration method
#'                          to use for hierarchical clustering. Default is "ward.D".}
#'   }
#'
#' @return A list containing:
#'         \item{cvols}{A list of \code{\linkS4class{ClusteredNeuroVol}} instances.}
#'         \item{cuts}{The number of cluster cuts.}
#'         \item{cutmat}{A matrix representing the cluster assignments for each cut.}
#'         \item{hclus}{The hierarchical clustering result.}
#'         \item{original_result}{The original clustering result (optional reference).}
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
meta_clust.cluster_result <- function(x, cuts = NULL, ...) {
  # Extract additional arguments
  dots <- list(...)
  algo <- if (!is.null(dots$algo)) dots$algo else "hclust"
  hclust_method <- if (!is.null(dots$hclust_method)) dots$hclust_method else "ward.D"
  dist_method <- if (!is.null(dots$dist_method)) dots$dist_method else "correlation"

  # Check if centers exist (some methods like SNIC don't provide centers)
  if (is.null(x$centers)) {
    stop("Cannot perform meta clustering: input 'x' does not contain cluster centers. ",
         "Methods like SNIC do not compute explicit cluster centers.")
  }

  # Fix default cuts logic based on actual number of clusters (rows)
  n_clusters <- nrow(x$centers)
  if (is.null(cuts)) {
    # Default: a few representative cuts up to half the number of clusters
    max_cut <- max(2, floor(n_clusters / 2))
    cuts <- sort(unique(c(2, min(5, max_cut), min(10, max_cut))))
    cuts <- cuts[cuts < n_clusters]
  }

  if (any(cuts >= n_clusters) || any(cuts < 2)) {
    stop("Cuts must be integers between 2 and ", n_clusters - 1, ".")
  }

  orig_labels <- x$cluster

  if (algo == "hclust") {
    # Compute Distance Matrix
    if (dist_method == "correlation") {
      # 1 - Pearson Correlation
      # Transpose needed because cor() works on columns (variables)
      D <- as.dist(1 - stats::cor(t(x$centers)))
    } else if (dist_method == "euclidean") {
      D <- stats::dist(x$centers)
    } else {
      stop("Unknown dist_method '", dist_method, "'. Use 'correlation' or 'euclidean'.")
    }

    # Hierarchical Clustering
    hres <- stats::hclust(D, method = hclust_method)

    # Cut tree at specified levels
    # cutree returns a matrix if 'k' is a vector
    cut_assignments <- stats::cutree(hres, k = cuts)

    # Ensure matrix format even if single cut
    if (is.vector(cut_assignments)) {
      cut_assignments <- matrix(cut_assignments, ncol = 1)
    }
    colnames(cut_assignments) <- paste0("k_", cuts)

    # Map meta-clusters back to voxel space
    # cut_assignments is (K_clusters x N_cuts)
    # orig_labels is (N_voxels) containing values 1..K
    # Fast mapping via matrix indexing
    cmat <- cut_assignments[orig_labels, , drop = FALSE]

    # Generate ClusteredNeuroVol objects
    cvols <- lapply(seq_along(cuts), function(i) {
      new_labels <- cmat[, i]
      neuroim2::ClusteredNeuroVol(x$clusvol@mask, cluster = new_labels)
    })

    # Return structure
    list(
      cvols = cvols,
      cuts = cuts,
      cutmat = cmat,
      hclus = hres,
      original_result = x
    )

  } else {
    stop(sprintf("Algorithm '%s' is not supported. Currently only 'hclust' is implemented.", algo))
  }
}

#' Merge Clustering Results for ClusteredNeuroVol Objects
#'
#' This method of merge_clus is specifically designed to merge clustering results for \code{ClusteredNeuroVol} objects
#' by performing consensus clustering across multiple clustering results.
#'
#' @param x A \code{ClusteredNeuroVol} object or an object of class \code{"cluster_result"}.
#' @param method A character string indicating the consensus clustering algorithm to use. Default is "SE".
#' See \code{\link[clue]{cl_consensus}} for available methods.
#' @param ... Additional clustering results to be merged.
#'
#' @return A \code{\linkS4class{ClusteredNeuroVol}} instance representing the consensus partition.
#'
#' @seealso \code{\link[clue]{cl_consensus}}, \code{\link[clue]{as.cl_hard_partition}}, \code{\link[clue]{cl_ensemble}}
#' @importFrom clue cl_consensus as.cl_hard_partition cl_ensemble
#' @importFrom assertthat assert_that
#' @export
merge_clus.cluster_result <- function(x, method="SE", ...) {
  args <- c(list(x), list(...))

  # Use vapply for lighter dependency checking than purrr
  is_cluster_result <- vapply(args, function(obj) inherits(obj, "cluster_result"), logical(1))
  assertthat::assert_that(all(is_cluster_result),
                          msg = "All arguments must be of class 'cluster_result'")

  # Create Ensemble
  ens <- do.call(clue::cl_ensemble, args)

  # Compute Consensus
  cons <- clue::cl_consensus(ens, method = method)

  # Extract Hard Partition
  hpart <- clue::as.cl_hard_partition(cons)
  consensus_labels <- as.integer(hpart$.Data)

  # Return ClusteredNeuroVol as promised by documentation
  neuroim2::ClusteredNeuroVol(x$clusvol@mask, cluster = consensus_labels)
}


#' @export
merge_clus.cluster_result_time <- function(x, ...) {
  # This variant returns integer vector (for internal/time-series use)
  args <- c(list(x), list(...))

  # Use vapply for lighter dependency checking
  is_cluster_result <- vapply(args, function(obj) inherits(obj, "cluster_result"), logical(1))
  assertthat::assert_that(all(is_cluster_result))

  ens <- do.call(clue::cl_ensemble, args)
  # Default consensus method typically used if not specified
  cons <- clue::cl_consensus(ens)
  hpart <- clue::as.cl_hard_partition(cons)

  # Return integer vector for this variant
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
#' @importFrom clue cl_class_ids
cl_class_ids.cluster_result <- function(x) {
  as.integer(x$cluster)
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
#' @importFrom clue is.cl_partition
is.cl_partition.cluster_result <- function(x) {
  TRUE
}
