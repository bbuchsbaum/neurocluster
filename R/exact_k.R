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

#' Build the masked-grid graph used by every exact-K caller.
#'
#' Returns the upper-triangular edge list, a CSR adjacency (`neighbor_ptr` /
#' `neighbor_idx`), the weak-component membership, and an igraph handle for
#' callers that still want one. The edge list and CSR come straight from C++;
#' no N-by-N sparse matrix is materialised, so there is no 65535-voxel ceiling.
#'
#' @keywords internal
#' @noRd
.exact_k_graph <- function(mask, connectivity) {
  method <- "force_exact_k"
  connectivity <- .cluster4d_scalar_number(
    connectivity, "connectivity", method, integer = TRUE
  )
  if (!connectivity %in% c(6L, 18L, 26L)) {
    stop("force_exact_k: connectivity must be 6, 18, or 26", call. = FALSE)
  }
  mask_idx <- which(cluster4d_mask_array(mask, method))
  native <- build_grid_edges_cpp(
    mask_idx = as.integer(mask_idx),
    dims = as.integer(dim(mask)),
    connectivity = as.integer(connectivity)
  )
  edges <- native$edges
  if (!nrow(edges)) edges <- matrix(integer(), ncol = 2L)

  graph <- igraph::make_empty_graph(length(mask_idx), directed = FALSE)
  if (nrow(edges)) {
    graph <- igraph::add_edges(graph, as.vector(t(edges)))
  }

  list(
    dims = as.integer(dim(mask)),
    mask_idx = mask_idx,
    edges = edges,
    graph = graph,
    neighbor_ptr = native$neighbor_ptr,
    neighbor_idx = native$neighbor_idx,
    components = native$components,
    n_voxels = length(mask_idx),
    connectivity = connectivity
  )
}

#' Validate that a caller-supplied graph matches the mask and connectivity.
#' @keywords internal
#' @noRd
.exact_k_resolve_graph <- function(graph_info, mask, connectivity) {
  if (is.null(graph_info)) return(.exact_k_graph(mask, connectivity))
  required <- c("dims", "mask_idx", "edges", "neighbor_ptr", "neighbor_idx",
                "components", "n_voxels", "connectivity")
  if (!is.list(graph_info) || !all(required %in% names(graph_info))) {
    stop("force_exact_k: graph_info is not an .exact_k_graph() result",
         call. = FALSE)
  }
  if (!identical(as.integer(graph_info$connectivity), as.integer(connectivity))) {
    stop("force_exact_k: graph_info connectivity does not match", call. = FALSE)
  }
  if (!identical(as.integer(graph_info$dims), as.integer(dim(mask)))) {
    stop("force_exact_k: graph_info dimensions do not match mask", call. = FALSE)
  }
  mask_idx <- which(cluster4d_mask_array(mask, "force_exact_k"))
  if (!identical(as.integer(graph_info$mask_idx), as.integer(mask_idx))) {
    stop("force_exact_k: graph_info mask does not match", call. = FALSE)
  }
  if (!identical(as.integer(graph_info$n_voxels), as.integer(length(mask_idx)))) {
    stop("force_exact_k: graph_info voxel count does not match mask", call. = FALSE)
  }
  graph_info
}

.exact_k_relabel <- function(labels) {
  ids <- sort(unique(labels))
  as.integer(match(labels, ids))
}

.exact_k_connected_labels <- function(labels, graph, edges) {
  # `graph` is accepted for backwards compatibility; components come from a
  # single union-find pass over the edge list, which needs no igraph object.
  edge_label_components_cpp(as.integer(labels), edges)
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

#' Ward cost of merging each (lo, hi) cluster pair.
#' @keywords internal
#' @noRd
.exact_k_pair_cost <- function(sums, counts, lo, hi) {
  if (!length(lo)) return(numeric())
  n_left <- counts[lo]
  n_right <- counts[hi]
  differences <- sums[lo, , drop = FALSE] / n_left -
    sums[hi, , drop = FALSE] / n_right
  n_left * n_right / (n_left + n_right) *
    .rowSums(differences * differences, length(lo), ncol(sums))
}

#' Merge adjacent clusters down to `K_target`.
#'
#' Same greedy adjacent-Ward rule and the same (cost, left, right) tie-break as
#' `.exact_k_select_merge()`, but the cluster sums, counts, and the boundary
#' pair set are all carried forward across merges instead of being rebuilt from
#' every edge on every iteration. That turns an O(merges * E) loop into
#' O(E + merges * boundary_degree).
#'
#' @keywords internal
#' @noRd
.exact_k_merge_to_target <- function(labels, feature_mat, K_target, edges,
                                     minimum_k, n_voxels) {
  k <- max(labels)
  if (k <= K_target) return(labels)

  counts <- tabulate(labels, nbins = k)
  sums <- rowsum(feature_mat, labels, reorder = TRUE)
  dimnames(sums) <- NULL

  left <- labels[edges[, 1L]]
  right <- labels[edges[, 2L]]
  boundary <- left != right
  lo <- pmin.int(left[boundary], right[boundary])
  hi <- pmax.int(left[boundary], right[boundary])
  keep <- !duplicated(as.numeric(lo) * (k + 1) + hi)
  lo <- lo[keep]
  hi <- hi[keep]
  cost <- .exact_k_pair_cost(sums, counts, lo, hi)

  while (k > K_target) {
    if (!length(lo)) {
      .exact_k_abort(
        "no_adjacent_merge", K_target, minimum_k, n_voxels,
        "no adjacent cluster pair remains"
      )
    }
    best <- which(cost == min(cost))
    if (length(best) > 1L) best <- best[order(lo[best], hi[best])[1L]]
    target <- lo[[best]]
    donor <- hi[[best]]

    counts[target] <- counts[target] + counts[donor]
    sums[target, ] <- sums[target, ] + sums[donor, ]
    labels[labels == donor] <- target

    # Redirect the donor's pairs onto the target, then drop the pair that just
    # became a self-loop plus any duplicates the redirection created.
    lo[lo == donor] <- target
    hi[hi == donor] <- target
    flipped <- lo > hi
    if (any(flipped)) {
      swap <- lo[flipped]
      lo[flipped] <- hi[flipped]
      hi[flipped] <- swap
    }
    alive <- lo != hi
    lo <- lo[alive]
    hi <- hi[alive]
    cost <- cost[alive]
    touched <- which(lo == target | hi == target)
    if (length(touched) > 1L) {
      duplicated_pairs <- duplicated(
        as.numeric(lo[touched]) * (k + 1) + hi[touched]
      )
      if (any(duplicated_pairs)) {
        drop <- touched[duplicated_pairs]
        lo <- lo[-drop]
        hi <- hi[-drop]
        cost <- cost[-drop]
      }
    }

    # Compact the label space exactly as .exact_k_relabel() would: the donor id
    # disappears and every larger id shifts down by one.
    labels <- labels - (labels > donor)
    lo <- lo - (lo > donor)
    hi <- hi - (hi > donor)
    target <- target - (target > donor)
    counts <- counts[-donor]
    sums <- sums[-donor, , drop = FALSE]
    k <- k - 1L

    touched <- which(lo == target | hi == target)
    if (length(touched)) {
      cost[touched] <- .exact_k_pair_cost(
        sums, counts, lo[touched], hi[touched]
      )
    }
  }
  as.integer(labels)
}

.exact_k_select_split <- function(labels, feature_mat, ptr, idx) {
  n <- length(labels)
  # Scratch is allocated once and modified in place; passing it into a helper
  # would make R copy all three vectors on every cluster, which is what made
  # this O(K * N) rather than O(N + E).
  visited <- logical(n)
  parent <- integer(n)
  child_count <- integer(n)
  queue <- integer(n)

  best_gain <- -Inf
  best <- NULL
  nodes_by_cluster <- split(seq_len(n), labels)

  for (cluster_id in seq_len(max(labels))) {
    nodes <- nodes_by_cluster[[as.character(cluster_id)]]
    if (is.null(nodes) || length(nodes) <= 1L) next

    # Breadth-first spanning tree of the cluster, rooted at its lowest voxel.
    root <- nodes[[1L]]
    head <- 1L
    tail <- 1L
    queue[tail] <- root
    visited[root] <- TRUE
    while (head <= tail) {
      node <- queue[head]
      head <- head + 1L
      from <- ptr[[node]]
      to <- ptr[[node + 1L]]
      if (to <= from) next
      adjacent <- idx[(from + 1L):to]
      adjacent <- adjacent[labels[adjacent] == cluster_id & !visited[adjacent]]
      for (next_node in adjacent) {
        visited[next_node] <- TRUE
        parent[next_node] <- node
        child_count[node] <- child_count[node] + 1L
        tail <- tail + 1L
        queue[tail] <- next_node
      }
    }

    if (!all(visited[nodes])) {
      .exact_k_abort(
        "disconnected_input_label", max(labels), 1L, n,
        paste0("label ", cluster_id, " is not connected after normalization")
      )
    }
    leaves <- nodes[child_count[nodes] == 0L & nodes != root]
    if (!length(leaves) && length(nodes) == 2L) leaves <- max(nodes)

    visited[nodes] <- FALSE
    parent[nodes] <- 0L
    child_count[nodes] <- 0L
    if (!length(leaves)) next

    center <- colMeans(feature_mat[nodes, , drop = FALSE])
    differences <- feature_mat[leaves, , drop = FALSE] -
      matrix(center, nrow = length(leaves), ncol = ncol(feature_mat),
             byrow = TRUE)
    gains <- length(nodes) / (length(nodes) - 1) *
      .rowSums(differences * differences, length(leaves), ncol(feature_mat))
    pick <- order(-gains, leaves)[1L]
    gain <- gains[[pick]]
    # Clusters are visited in ascending id and only a strictly larger gain
    # wins, which reproduces order(-gains, cluster, voxel)[1].
    if (gain > best_gain) {
      best_gain <- gain
      best <- list(
        cluster = as.integer(cluster_id),
        voxel = as.integer(leaves[[pick]]),
        gain = as.numeric(gain)
      )
    }
  }
  best
}

#' @rdname exact_k
#' @param labels Positive integer labels in masked-voxel order.
#' @param feature_mat Numeric matrix with voxels in rows.
#' @param K_target Desired number of connected labels.
#' @param mask NeuroVol defining voxel geometry and inclusion.
#' @param connectivity Grid connectivity, one of 6, 18, or 26.
#' @param graph_info Optional prebuilt `.exact_k_graph()` result for the same
#'   mask and connectivity. Supplying it avoids rebuilding the masked-grid
#'   graph, which callers that already hold one would otherwise pay for twice.
#'   The receipt is rejected if its dimensions, included voxels, or connectivity
#'   do not match `mask`.
#' @return Contiguous positive integer labels with exactly `K_target`
#'   connected components, or a `cluster4d_exact_k_infeasible` condition.
#' @keywords internal
force_exact_k <- function(labels, feature_mat, K_target, mask,
                          connectivity = 6L, graph_info = NULL) {
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
  graph_info <- .exact_k_resolve_graph(graph_info, mask, connectivity)
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

  minimum_k <- length(unique(as.integer(graph_info$components)))
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

  if (current_k > K_target) {
    labels <- .exact_k_merge_to_target(
      labels, feature_mat, K_target, graph_info$edges, minimum_k, n_voxels
    )
    current_k <- max(labels)
  }

  while (current_k < K_target) {
    split <- .exact_k_select_split(
      labels, feature_mat, graph_info$neighbor_ptr, graph_info$neighbor_idx
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
