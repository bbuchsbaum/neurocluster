#' Adjacency-preserving exact-K repair
#'
#' Internal machinery for converting a masked-voxel partition to exactly K
#' connected labels. Over-target partitions use only adjacent Ward merges;
#' under-target partitions split deterministic non-articulation leaves from a
#' spanning tree. Disconnected pieces of an input label are normalized before
#' either operation.
#'
#' @keywords internal
#' @name exact_k
NULL

.exact_k_abort <- function(reason, target_k, minimum_k, maximum_k, detail) {
  condition <- structure(
    list(
      message = paste0("Exact K is infeasible: ", detail),
      call = NULL,
      reason = reason,
      target_k = as.integer(target_k),
      minimum_k = as.integer(minimum_k),
      maximum_k = as.integer(maximum_k)
    ),
    class = c("cluster4d_exact_k_infeasible", "error", "condition")
  )
  stop(condition)
}

.exact_k_graph <- function(mask, connectivity) {
  method <- "force_exact_k"
  connectivity <- .cluster4d_scalar_number(
    connectivity, "connectivity", method, integer = TRUE
  )
  if (!connectivity %in% c(6L, 18L, 26L)) {
    stop("force_exact_k: connectivity must be 6, 18, or 26", call. = FALSE)
  }
  mask_idx <- which(cluster4d_mask_array(mask, method))
  adjacency <- build_grid_adjacency(mask, mask_idx, connectivity)
  summary <- Matrix::summary(adjacency)
  keep <- summary$i < summary$j
  edges <- cbind(
    as.integer(summary$i[keep]),
    as.integer(summary$j[keep])
  )
  if (!length(edges)) edges <- matrix(integer(), ncol = 2L)

  graph <- igraph::make_empty_graph(length(mask_idx), directed = FALSE)
  if (nrow(edges)) {
    graph <- igraph::add_edges(graph, as.vector(t(edges)))
  }
  from <- c(edges[, 1L], edges[, 2L])
  to <- c(edges[, 2L], edges[, 1L])
  neighbors <- split(
    to,
    factor(from, levels = seq_len(length(mask_idx)))
  )
  neighbors <- lapply(neighbors, function(x) sort(as.integer(x)))

  list(
    mask_idx = mask_idx,
    edges = edges,
    graph = graph,
    neighbors = neighbors,
    connectivity = connectivity
  )
}

.exact_k_relabel <- function(labels) {
  ids <- sort(unique(labels))
  as.integer(match(labels, ids))
}

.exact_k_connected_labels <- function(labels, graph, edges) {
  same <- if (nrow(edges)) labels[edges[, 1L]] == labels[edges[, 2L]] else logical()
  label_graph <- igraph::make_empty_graph(length(labels), directed = FALSE)
  if (any(same)) {
    label_graph <- igraph::add_edges(
      label_graph, as.vector(t(edges[same, , drop = FALSE]))
    )
  }
  .exact_k_relabel(as.integer(igraph::components(label_graph)$membership))
}

.exact_k_cluster_summaries <- function(labels, feature_mat) {
  k <- max(labels)
  counts <- tabulate(labels, nbins = k)
  sums <- rowsum(feature_mat, labels, reorder = TRUE)
  centers <- sums / counts
  dimnames(centers) <- NULL
  list(counts = counts, centers = centers)
}

.exact_k_boundary_pairs <- function(labels, edges) {
  if (!nrow(edges)) return(matrix(integer(), ncol = 2L))
  left <- labels[edges[, 1L]]
  right <- labels[edges[, 2L]]
  boundary <- left != right
  if (!any(boundary)) return(matrix(integer(), ncol = 2L))
  pairs <- cbind(
    pmin(left[boundary], right[boundary]),
    pmax(left[boundary], right[boundary])
  )
  unique(pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE])
}

.exact_k_select_merge <- function(labels, feature_mat, edges) {
  pairs <- .exact_k_boundary_pairs(labels, edges)
  if (!nrow(pairs)) return(NULL)
  summaries <- .exact_k_cluster_summaries(labels, feature_mat)
  left <- pairs[, 1L]
  right <- pairs[, 2L]
  differences <- summaries$centers[left, , drop = FALSE] -
    summaries$centers[right, , drop = FALSE]
  costs <- summaries$counts[left] * summaries$counts[right] /
    (summaries$counts[left] + summaries$counts[right]) *
    rowSums(differences^2)
  selected <- order(costs, left, right)[1L]
  list(
    left = as.integer(left[selected]),
    right = as.integer(right[selected]),
    cost = as.numeric(costs[selected])
  )
}

.exact_k_tree_leaves <- function(nodes, neighbors, labels, cluster_id) {
  root <- min(nodes)
  parent <- integer(length(labels))
  visited <- logical(length(labels))
  queue <- integer(length(nodes))
  head <- 1L
  tail <- 1L
  queue[tail] <- root
  visited[root] <- TRUE

  while (head <= tail) {
    node <- queue[head]
    head <- head + 1L
    adjacent <- neighbors[[node]]
    if (!length(adjacent)) next
    adjacent <- adjacent[labels[adjacent] == cluster_id & !visited[adjacent]]
    for (next_node in adjacent) {
      visited[next_node] <- TRUE
      parent[next_node] <- node
      tail <- tail + 1L
      queue[tail] <- next_node
    }
  }
  if (!all(visited[nodes])) {
    .exact_k_abort(
      "disconnected_input_label", max(labels), 1L, length(labels),
      paste0("label ", cluster_id, " is not connected after normalization")
    )
  }
  child_count <- tabulate(parent[parent > 0L], nbins = length(labels))
  leaves <- nodes[child_count[nodes] == 0L & nodes != root]
  if (!length(leaves) && length(nodes) == 2L) leaves <- max(nodes)
  leaves
}

.exact_k_select_split <- function(labels, feature_mat, neighbors) {
  candidates <- vector("list", max(labels))
  for (cluster_id in seq_len(max(labels))) {
    nodes <- which(labels == cluster_id)
    if (length(nodes) <= 1L) next
    leaves <- .exact_k_tree_leaves(nodes, neighbors, labels, cluster_id)
    if (!length(leaves)) next
    center <- colMeans(feature_mat[nodes, , drop = FALSE])
    differences <- feature_mat[leaves, , drop = FALSE] -
      matrix(center, nrow = length(leaves), ncol = ncol(feature_mat), byrow = TRUE)
    gains <- length(nodes) / (length(nodes) - 1) * rowSums(differences^2)
    best <- order(-gains, leaves)[1L]
    candidates[[cluster_id]] <- list(
      cluster = as.integer(cluster_id),
      voxel = as.integer(leaves[best]),
      gain = as.numeric(gains[best])
    )
  }
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) return(NULL)
  gains <- vapply(candidates, `[[`, numeric(1), "gain")
  clusters <- vapply(candidates, `[[`, integer(1), "cluster")
  voxels <- vapply(candidates, `[[`, integer(1), "voxel")
  candidates[[order(-gains, clusters, voxels)[1L]]]
}

#' @rdname exact_k
#' @param labels Positive integer labels in masked-voxel order.
#' @param feature_mat Numeric matrix with voxels in rows.
#' @param K_target Desired number of connected labels.
#' @param mask NeuroVol defining voxel geometry and inclusion.
#' @param connectivity Grid connectivity, one of 6, 18, or 26.
#' @return Contiguous positive integer labels with exactly `K_target`
#'   connected components, or a `cluster4d_exact_k_infeasible` condition.
#' @keywords internal
force_exact_k <- function(labels, feature_mat, K_target, mask,
                          connectivity = 6L) {
  if (!inherits(mask, "NeuroVol")) {
    stop("force_exact_k: mask must be a NeuroVol", call. = FALSE)
  }
  if (!is.numeric(feature_mat)) {
    stop("force_exact_k: feature_mat must be numeric", call. = FALSE)
  }
  feature_mat <- as.matrix(feature_mat)
  if (any(!is.finite(feature_mat))) {
    stop("force_exact_k: feature_mat must be finite", call. = FALSE)
  }
  if (!is.numeric(labels) || any(!is.finite(labels)) ||
      any(labels != floor(labels)) || any(labels <= 0)) {
    stop("force_exact_k: labels must be finite positive integers", call. = FALSE)
  }
  labels <- as.integer(labels)
  graph_info <- .exact_k_graph(mask, connectivity)
  n_voxels <- length(graph_info$mask_idx)
  if (length(labels) != n_voxels || nrow(feature_mat) != n_voxels) {
    stop(
      "force_exact_k: labels and feature rows must equal included mask voxels",
      call. = FALSE
    )
  }
  K_target <- .cluster4d_scalar_number(
    K_target, "K_target", "force_exact_k",
    lower = 1, upper = n_voxels, integer = TRUE
  )

  mask_components <- as.integer(igraph::components(graph_info$graph)$membership)
  minimum_k <- length(unique(mask_components))
  if (K_target < minimum_k) {
    .exact_k_abort(
      "disconnected_mask_components", K_target, minimum_k, n_voxels,
      paste0(
        "target ", K_target, " is below the ", minimum_k,
        " disconnected mask components required by the topology"
      )
    )
  }

  labels <- .exact_k_connected_labels(
    labels, graph_info$graph, graph_info$edges
  )
  current_k <- max(labels)

  while (current_k > K_target) {
    merge <- .exact_k_select_merge(labels, feature_mat, graph_info$edges)
    if (is.null(merge)) {
      .exact_k_abort(
        "no_adjacent_merge", K_target, minimum_k, n_voxels,
        "no adjacent cluster pair remains"
      )
    }
    labels[labels == merge$right] <- merge$left
    labels <- .exact_k_relabel(labels)
    current_k <- current_k - 1L
  }

  while (current_k < K_target) {
    split <- .exact_k_select_split(
      labels, feature_mat, graph_info$neighbors
    )
    if (is.null(split)) {
      .exact_k_abort(
        "no_connected_split", K_target, minimum_k, n_voxels,
        "no connected cluster with more than one voxel remains"
      )
    }
    current_k <- current_k + 1L
    labels[split$voxel] <- current_k
  }

  .exact_k_relabel(labels)
}
