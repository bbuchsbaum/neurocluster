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

  orig_labels <- x$cluster

  # Fix default cuts logic based on actual number of clusters (rows)
  n_clusters <- nrow(x$centers)
  labels_used <- sort(unique(orig_labels))
  # ensure centers align with labels_used; if not, assume centers rows correspond to labels_used order
  if (length(labels_used) != n_clusters) {
    labels_used <- seq_len(n_clusters)
  }
  if (is.null(cuts)) {
    # Default: a few representative cuts up to half the number of clusters
    max_cut <- max(2, floor(n_clusters / 2))
    cuts <- sort(unique(c(2, min(5, max_cut), min(10, max_cut))))
    cuts <- cuts[cuts < n_clusters]
  }

  if (any(cuts >= n_clusters) || any(cuts < 2)) {
    stop("Cuts must be integers between 2 and ", n_clusters - 1, ".")
  }

  if (algo == "hclust") {
    # Compute Distance Matrix
    if (dist_method == "correlation") {
      # 1 - Pearson Correlation
      # Transpose needed because cor() works on columns (variables)
      C <- stats::cor(t(x$centers))
      C[!is.finite(C)] <- 0  # handle zero-variance clusters
      diag(C) <- 1
      D <- as.dist(1 - C)
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
    # orig_labels is (N_voxels) containing values 1..K (contiguous)
    label_index <- match(orig_labels, labels_used)
    cmat <- cut_assignments[label_index, , drop = FALSE]

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

  # All inputs must have the same number of elements
  n_elems <- vapply(args, function(obj) length(obj$cluster), integer(1))
  assertthat::assert_that(length(unique(n_elems)) == 1,
                          msg = "All cluster_result inputs must have the same number of labeled elements.")

  # Create Ensemble
  ens <- do.call(clue::cl_ensemble, args)

  # Compute Consensus
  cons <- clue::cl_consensus(ens, method = method)

  # Extract Hard Partition
  hpart <- clue::as.cl_hard_partition(cons)
  consensus_labels <- as.integer(hpart$.Data)

  # Return ClusteredNeuroVol as promised by documentation
  mask_obj <- if (!is.null(x$clusvol)) x$clusvol@mask else NULL
  if (!is.null(mask_obj)) {
    neuroim2::ClusteredNeuroVol(mask_obj, cluster = consensus_labels)
  } else {
    # Fallback: return a cluster_result-like list when mask missing
    list(cluster = consensus_labels, method = "consensus", n_clusters = length(unique(consensus_labels)))
  }
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

#' Lightweight cluster-quality metrics
#'
#' Computes quick diagnostic metrics for a clustering result, supporting both
#' data-driven and spatial assessments.
#'
#' @param result A \code{cluster_result} (or compatible list with \code{cluster} and optionally \code{centers}, \code{coord_centers}).
#' @param feature_mat Optional matrix of features (voxels x time/features). If NULL, tries \code{result$metadata$features} or \code{result$data_prep$features}.
#' @param coords Optional matrix of coordinates (voxels x 3). If NULL, tries \code{result$metadata$coords} or \code{result$data_prep$coords}.
#' @param truth Optional integer vector of ground-truth labels for ARI.
#'
#' @return Named list with fields (when available):
#' \describe{
#'   \item{n_clusters}{Number of clusters.}
#'   \item{size_summary}{Min/median/max cluster size.}
#'   \item{ari_truth}{Adjusted Rand Index vs. \code{truth} (if provided).}
#'   \item{within_data_corr}{Mean Fisher-z of voxel-to-centroid correlations (data space).}
#'   \item{between_data_corr}{Mean Fisher-z of centroid-to-centroid correlations (data space).}
#'   \item{within_spatial_dist}{Mean Euclidean distance voxel-to-centroid (spatial).}
#' }
#' @export
cluster_metrics <- function(result,
                            feature_mat = NULL,
                            coords = NULL,
                            truth = NULL) {

  labs <- result$cluster
  stopifnot(!is.null(labs))
  n <- length(labs)
  n_clusters <- length(unique(labs))

  # Pull features/coords if not provided
  if (is.null(feature_mat)) {
    feature_mat <- result$metadata$features %||% result$data_prep$features
  }
  if (is.null(coords)) {
    coords <- result$metadata$coords %||% result$data_prep$coords
  }

  fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))

  # Size summary
  sz <- tabulate(labs, nbins = n_clusters)
  size_summary <- c(min = min(sz), median = stats::median(sz), max = max(sz))

  out <- list(
    n_clusters = n_clusters,
    size_summary = size_summary
  )

  # ARI if truth supplied
  if (!is.null(truth) && length(truth) == n) {
    out$ari_truth <- adj_rand_index_internal(labs, truth)
  }

  # Data-based metrics
  if (!is.null(feature_mat)) {
    # Ensure matrix is voxels x features
    if (nrow(feature_mat) != n && ncol(feature_mat) == n) {
      feature_mat <- t(feature_mat)
    }
    stopifnot(nrow(feature_mat) == n)

    # Use provided centers if present, else compute means
    centers <- result$centers
    if (is.null(centers)) {
      centers <- vapply(seq_len(n_clusters), function(k) colMeans(feature_mat[labs == k, , drop = FALSE]), numeric(ncol(feature_mat)))
      centers <- t(centers)
    }
    if (nrow(centers) != n_clusters) {
      # attempt to reorder/trim to contiguous ids
      centers <- centers[seq_len(n_clusters), , drop = FALSE]
    }

    # voxel-to-centroid correlation
    # normalize
    fm <- scale(feature_mat, center = TRUE, scale = TRUE)
    cm <- scale(centers, center = TRUE, scale = TRUE)
    # efficient dot product correlation
    denom_v <- sqrt(rowSums(fm^2))
    denom_c <- sqrt(rowSums(cm^2))
    denom_v[denom_v == 0] <- 1
    denom_c[denom_c == 0] <- 1
    # compute per voxel corr with its centroid
    within_corr <- numeric(n)
    for (k in seq_len(n_clusters)) {
      idx <- which(labs == k)
      if (length(idx) == 0) next
      within_corr[idx] <- as.numeric(fm[idx, , drop = FALSE] %*% cm[k, ]) /
        (denom_v[idx] * denom_c[k])
    }
    within_corr[!is.finite(within_corr)] <- 0
    out$within_data_corr <- mean(fisher_z(within_corr))

    # centroid-to-centroid correlation
    cc <- stats::cor(t(cm))
    cc[!is.finite(cc)] <- 0
    out$between_data_corr <- mean(fisher_z(cc[upper.tri(cc)]))
  }

  # Spatial metrics
  if (!is.null(coords)) {
    if (nrow(coords) != n && ncol(coords) == n) coords <- t(coords)
    stopifnot(nrow(coords) == n)
    ccent <- result$coord_centers
    if (is.null(ccent)) {
      ccent <- vapply(seq_len(n_clusters), function(k) colMeans(coords[labs == k, , drop = FALSE]), numeric(ncol(coords)))
      ccent <- t(ccent)
    }
    if (nrow(ccent) != n_clusters) {
      ccent <- ccent[seq_len(n_clusters), , drop = FALSE]
    }
    within_dist <- numeric(n)
    for (k in seq_len(n_clusters)) {
      idx <- which(labs == k)
      if (length(idx) == 0) next
      diff <- coords[idx, , drop = FALSE] - matrix(rep(ccent[k, ], each = length(idx)), ncol = ncol(coords), byrow = TRUE)
      within_dist[idx] <- sqrt(rowSums(diff^2))
    }
    out$within_spatial_dist <- mean(within_dist)
  }

  out
}

# Internal ARI (duplicated here to avoid extra deps)
adj_rand_index_internal <- function(labels1, labels2) {
  if (length(labels1) != length(labels2)) return(NA_real_)
  tab <- table(labels1, labels2)
  sum_comb <- sum(choose(tab, 2))
  sum_rows <- sum(choose(rowSums(tab), 2))
  sum_cols <- sum(choose(colSums(tab), 2))
  n <- length(labels1)
  expected <- sum_rows * sum_cols / choose(n, 2)
  max_idx <- 0.5 * (sum_rows + sum_cols)
  if (max_idx == expected) return(0)
  (sum_comb - expected) / (max_idx - expected)
}

# Null-coalescing helper
`%||%` <- function(a, b) if (!is.null(a)) a else b

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
