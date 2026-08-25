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
#'         \item{cvols}{A list of \code{\link[neuroim2:ClusteredNeuroVol-class]{ClusteredNeuroVol}} instances.}
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
#'    - O(K^2) pairwise correlations where K = number of input clusters
#'    - Sequential double loop in base R
#' 
#' 2. **Hierarchical Clustering**: `hclust()`
#'    - Standard hierarchical agglomeration (Ward.D by default)
#'    - Sequential merging of closest clusters
#'    - O(K^2 log K) complexity
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
#' - **Memory efficient**: Only stores K x K distance matrix
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
#' - **Complexity**: O(K^2) for distances, O(K^2 log K) for clustering
#' - **Memory**: O(K^2) for distance matrix
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
#' @return A \code{\link[neuroim2:ClusteredNeuroVol-class]{ClusteredNeuroVol}} instance representing the consensus partition.
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

#' Explicit cluster-quality estimands
#'
#' Computes partition, within-cluster temporal, and physical spatial estimands
#' from final voxel labels. Temporal coherence is the mean Pearson correlation
#' over all unordered within-cluster voxel pairs. Spatial dispersion is the
#' root mean squared Euclidean distance from each voxel coordinate to its final
#' cluster centroid; lower values are more compact.
#'
#' @param result A `cluster_result` with one final label per voxel.
#' @param feature_mat Optional finite numeric matrix with voxels in rows and at
#'   least two time/features in columns. A transposed matrix is accepted when
#'   only its columns match the number of labels. Every row must be
#'   non-constant for Pearson correlation.
#' @param coords Optional finite numeric matrix with voxels in rows and spatial
#'   dimensions in columns. When omitted for a canonical `cluster4d_result`,
#'   physical millimetre coordinates are reconstructed from its mask and
#'   affine provenance.
#' @param truth Optional partition labels over exactly the same voxels.
#'
#' @return A named list. Always contains `n_clusters` and `size_summary`.
#'   With `truth`, it adds adjusted Rand, variation of information in bits, and
#'   pairwise Dice. With `feature_mat`, it adds
#'   `temporal_pairwise_correlation`. With coordinates, it adds
#'   `spatial_rms_distance_mm` (or `spatial_rms_distance` when explicit
#'   coordinates are not declared to be millimetres by a canonical result).
#' @export
cluster_metrics <- function(result,
                            feature_mat = NULL,
                            coords = NULL,
                            truth = NULL) {
  if (is.null(result$cluster)) {
    stop("result must contain final voxel labels in $cluster", call. = FALSE)
  }
  labs <- .partition_label_codes(result$cluster, "result$cluster")
  n <- length(labs)
  n_clusters <- max(labs)
  sizes <- tabulate(labs, nbins = n_clusters)
  size_summary <- c(
    min = min(sizes), median = stats::median(sizes), max = max(sizes)
  )
  out <- list(
    n_clusters = n_clusters,
    size_summary = size_summary
  )

  if (!is.null(truth)) {
    agreement <- .partition_metrics(labs, truth)
    out$ari_truth <- agreement$ari
    out$variation_of_information_truth_bits <-
      agreement$variation_of_information_bits
    out$pairwise_dice_truth <- agreement$pairwise_dice
  }

  if (!is.null(feature_mat)) {
    out$temporal_pairwise_correlation <-
      .temporal_pairwise_correlation(labs, feature_mat)
  }

  coords_are_mm <- FALSE
  if (is.null(coords) && inherits(result, "cluster4d_result")) {
    support <- .cluster4d_result_support(result)
    coords <- support$coords
    coords_are_mm <- TRUE
  }
  if (!is.null(coords)) {
    spatial_name <- if (coords_are_mm ||
                        (inherits(result, "cluster4d_result") &&
                         identical(
                           result$provenance$coordinate_space$units, "mm"
                         ))) {
      "spatial_rms_distance_mm"
    } else {
      "spatial_rms_distance"
    }
    out[[spatial_name]] <- .spatial_rms_dispersion(labs, coords)
  }
  out
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
