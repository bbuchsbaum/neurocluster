#' Adjacency-preserving exact-K repair
#'
#' Internal machinery for converting a masked-voxel partition to exactly K
#' connected labels. Over-target partitions use only adjacent Ward merges;
#' under-target partitions use deterministic feature-weighted spanning-tree
#' bisections chosen by Ward SSE reduction. Disconnected pieces of an input
#' label are normalized before either operation.
#'
#' @keywords internal
#' @name exact_k
NULL

.exact_k_abort <- function(reason, target_k, minimum_k, maximum_k, detail,
                           n_voxels = maximum_k, min_cluster_size = 1L,
                           current_k = NA_integer_) {
  condition <- structure(
    list(
      message = paste0("Exact K is infeasible: ", detail),
      call = NULL,
      reason = reason,
      target_k = as.integer(target_k),
      minimum_k = as.integer(minimum_k),
      maximum_k = as.integer(maximum_k),
      n_voxels = as.integer(n_voxels),
      min_cluster_size = as.integer(min_cluster_size),
      current_k = as.integer(current_k)
    ),
    class = c("cluster4d_exact_k_infeasible", "error", "condition")
  )
  stop(condition)
}

.exact_k_new_trace <- function(enabled) {
  if (!isTRUE(enabled)) return(NULL)
  trace <- new.env(parent = emptyenv())
  trace$operations <- list()
  trace
}

.exact_k_trace_add <- function(trace, operation) {
  if (is.null(trace)) return(invisible(NULL))
  trace$operations[[length(trace$operations) + 1L]] <- operation
  invisible(NULL)
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
                                     minimum_k, n_voxels,
                                     maximum_k = n_voxels,
                                     min_cluster_size = 1L, trace = NULL) {
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
        "no_adjacent_merge", K_target, minimum_k, maximum_k,
        "no adjacent cluster pair remains",
        n_voxels, min_cluster_size, k
      )
    }
    best <- which(cost == min(cost))
    if (length(best) > 1L) best <- best[order(lo[best], hi[best])[1L]]
    target <- lo[[best]]
    donor <- hi[[best]]

    .exact_k_trace_add(trace, list(
      type = "adjacent_ward_merge",
      left_label = as.integer(target),
      right_label = as.integer(donor),
      left_size = as.integer(counts[[target]]),
      right_size = as.integer(counts[[donor]]),
      ward_cost = as.numeric(cost[[best]])
    ))

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

.exact_k_merge_pair <- function(labels, left, right) {
  target <- min(left, right)
  donor <- max(left, right)
  labels[labels == donor] <- target
  .exact_k_relabel(labels)
}

.exact_k_prepare_minimum_size <- function(labels, feature_mat, edges, K_target,
                                          min_cluster_size, minimum_k,
                                          maximum_k, n_voxels, trace = NULL) {
  if (min_cluster_size <= 1L) return(labels)

  repeat {
    counts <- tabulate(labels, nbins = max(labels))
    undersized <- which(counts < min_cluster_size)
    if (!length(undersized)) break
    cluster <- undersized[order(counts[undersized], undersized)[1L]]
    pairs <- .exact_k_boundary_pairs(labels, edges)
    candidates <- which(pairs[, 1L] == cluster | pairs[, 2L] == cluster)
    if (!length(candidates)) {
      .exact_k_abort(
        "undersized_component", K_target, minimum_k, maximum_k,
        paste0(
          "cluster ", cluster, " has size ", counts[[cluster]],
          " below min_cluster_size=", min_cluster_size,
          " and has no adjacent cluster"
        ),
        n_voxels, min_cluster_size, max(labels)
      )
    }
    pairs <- pairs[candidates, , drop = FALSE]
    summaries <- .exact_k_cluster_summaries(labels, feature_mat)
    costs <- .exact_k_pair_cost(
      sweep(summaries$centers, 1L, summaries$counts, "*"),
      summaries$counts, pairs[, 1L], pairs[, 2L]
    )
    before_capacity <- sum(counts %/% min_cluster_size)
    merged_sizes <- counts[pairs[, 1L]] + counts[pairs[, 2L]]
    capacity_gain <- merged_sizes %/% min_cluster_size -
      counts[pairs[, 1L]] %/% min_cluster_size -
      counts[pairs[, 2L]] %/% min_cluster_size
    other <- ifelse(pairs[, 1L] == cluster, pairs[, 2L], pairs[, 1L])
    selected <- order(-capacity_gain, costs, other, pairs[, 1L], pairs[, 2L])[1L]
    .exact_k_trace_add(trace, list(
      type = "minimum_size_merge",
      left_label = as.integer(pairs[selected, 1L]),
      right_label = as.integer(pairs[selected, 2L]),
      left_size = as.integer(counts[[pairs[selected, 1L]]]),
      right_size = as.integer(counts[[pairs[selected, 2L]]]),
      capacity_gain = as.integer(capacity_gain[[selected]]),
      ward_cost = as.numeric(costs[[selected]])
    ))
    labels <- .exact_k_merge_pair(
      labels, pairs[selected, 1L], pairs[selected, 2L]
    )
    if (max(labels) < minimum_k ||
        sum(tabulate(labels) %/% min_cluster_size) < before_capacity) {
      .exact_k_abort(
        "minimum_size_repair_failed", K_target, minimum_k, maximum_k,
        "minimum-size normalization violated topology or split capacity",
        n_voxels, min_cluster_size, max(labels)
      )
    }
  }
  labels
}

.exact_k_consolidate_split_capacity <- function(
    labels, feature_mat, edges, K_target, min_cluster_size,
    minimum_k, maximum_k, n_voxels, trace = NULL) {
  if (min_cluster_size <= 1L) return(labels)
  repeat {
    counts <- tabulate(labels, nbins = max(labels))
    capacity <- sum(counts %/% min_cluster_size)
    if (capacity >= K_target) return(labels)
    pairs <- .exact_k_boundary_pairs(labels, edges)
    if (!nrow(pairs) || max(labels) <= minimum_k) {
      .exact_k_abort(
        "insufficient_split_capacity", K_target, minimum_k, maximum_k,
        paste0(
          "current clusters can support only ", capacity,
          " clusters of size at least ", min_cluster_size
        ),
        n_voxels, min_cluster_size, max(labels)
      )
    }
    summaries <- .exact_k_cluster_summaries(labels, feature_mat)
    costs <- .exact_k_pair_cost(
      sweep(summaries$centers, 1L, summaries$counts, "*"),
      summaries$counts, pairs[, 1L], pairs[, 2L]
    )
    merged_sizes <- counts[pairs[, 1L]] + counts[pairs[, 2L]]
    capacity_gain <- merged_sizes %/% min_cluster_size -
      counts[pairs[, 1L]] %/% min_cluster_size -
      counts[pairs[, 2L]] %/% min_cluster_size
    selected <- order(
      -capacity_gain, costs, pairs[, 1L], pairs[, 2L]
    )[1L]
    .exact_k_trace_add(trace, list(
      type = "split_capacity_merge",
      left_label = as.integer(pairs[selected, 1L]),
      right_label = as.integer(pairs[selected, 2L]),
      left_size = as.integer(counts[[pairs[selected, 1L]]]),
      right_size = as.integer(counts[[pairs[selected, 2L]]]),
      capacity_gain = as.integer(capacity_gain[[selected]]),
      ward_cost = as.numeric(costs[[selected]])
    ))
    labels <- .exact_k_merge_pair(
      labels, pairs[selected, 1L], pairs[selected, 2L]
    )
  }
}

.exact_k_feature_mst <- function(nodes, labels, cluster_id, feature_mat, edges) {
  inside <- labels[edges[, 1L]] == cluster_id &
    labels[edges[, 2L]] == cluster_id
  induced <- edges[inside, , drop = FALSE]
  if (nrow(induced) < length(nodes) - 1L) return(NULL)

  differences <- feature_mat[induced[, 1L], , drop = FALSE] -
    feature_mat[induced[, 2L], , drop = FALSE]
  costs <- rowSums(differences * differences)
  ordered <- order(
    costs,
    pmin(induced[, 1L], induced[, 2L]),
    pmax(induced[, 1L], induced[, 2L])
  )
  induced <- induced[ordered, , drop = FALSE]

  local_edges <- matrix(
    match(as.vector(induced), nodes), ncol = 2L, dimnames = NULL
  )
  parent <- seq_along(nodes)
  rank <- integer(length(nodes))
  find_root <- function(node) {
    while (node != parent[[node]]) {
      parent[[node]] <<- parent[[parent[[node]]]]
      node <- parent[[node]]
    }
    node
  }
  tree <- matrix(0L, nrow = length(nodes) - 1L, ncol = 2L)
  tree_count <- 0L
  for (row in seq_len(nrow(induced))) {
    left <- local_edges[row, 1L]
    right <- local_edges[row, 2L]
    root_left <- find_root(left)
    root_right <- find_root(right)
    if (root_left == root_right) next
    if (rank[[root_left]] < rank[[root_right]]) {
      swap <- root_left
      root_left <- root_right
      root_right <- swap
    }
    parent[[root_right]] <- root_left
    if (rank[[root_left]] == rank[[root_right]]) {
      rank[[root_left]] <- rank[[root_left]] + 1L
    }
    tree_count <- tree_count + 1L
    tree[tree_count, ] <- c(left, right)
    if (tree_count == length(nodes) - 1L) break
  }
  if (tree_count != length(nodes) - 1L) return(NULL)
  tree
}

.exact_k_cluster_split_candidate <- function(
    nodes, labels, cluster_id, feature_mat, edges, min_cluster_size,
    other_capacity, K_target) {
  if (length(nodes) < 2L * min_cluster_size) return(NULL)
  tree <- .exact_k_feature_mst(nodes, labels, cluster_id, feature_mat, edges)
  if (is.null(tree)) return(NULL)

  adjacency <- vector("list", length(nodes))
  for (row in seq_len(nrow(tree))) {
    left <- tree[row, 1L]
    right <- tree[row, 2L]
    adjacency[[left]] <- c(adjacency[[left]], right)
    adjacency[[right]] <- c(adjacency[[right]], left)
  }
  adjacency <- lapply(adjacency, function(neighbors) {
    neighbors[order(nodes[neighbors])]
  })

  root <- which.min(nodes)
  parent <- integer(length(nodes))
  traversal <- integer(length(nodes))
  queue <- integer(length(nodes))
  head <- 1L
  tail <- 1L
  queue[[tail]] <- root
  seen <- logical(length(nodes))
  seen[[root]] <- TRUE
  visited <- 0L
  while (head <= tail) {
    node <- queue[[head]]
    head <- head + 1L
    visited <- visited + 1L
    traversal[[visited]] <- node
    for (neighbor in adjacency[[node]]) {
      if (seen[[neighbor]]) next
      seen[[neighbor]] <- TRUE
      parent[[neighbor]] <- node
      tail <- tail + 1L
      queue[[tail]] <- neighbor
    }
  }
  if (visited != length(nodes)) return(NULL)

  subtree_count <- rep(1L, length(nodes))
  subtree_sum <- feature_mat[nodes, , drop = FALSE]
  subtree_min <- nodes
  for (node in rev(traversal[-1L])) {
    ancestor <- parent[[node]]
    subtree_count[[ancestor]] <- subtree_count[[ancestor]] +
      subtree_count[[node]]
    subtree_sum[ancestor, ] <- subtree_sum[ancestor, ] + subtree_sum[node, ]
    subtree_min[[ancestor]] <- min(subtree_min[[ancestor]], subtree_min[[node]])
  }

  candidates <- traversal[-1L]
  feasible <- subtree_count[candidates] >= min_cluster_size &
    length(nodes) - subtree_count[candidates] >= min_cluster_size
  capacity_after <- other_capacity +
    subtree_count[candidates] %/% min_cluster_size +
    (length(nodes) - subtree_count[candidates]) %/% min_cluster_size
  feasible <- feasible & capacity_after >= K_target
  candidates <- candidates[feasible]
  if (!length(candidates)) return(NULL)
  left_count <- subtree_count[candidates]
  right_count <- length(nodes) - left_count
  left_mean <- subtree_sum[candidates, , drop = FALSE] / left_count
  total_sum <- colSums(feature_mat[nodes, , drop = FALSE])
  right_mean <- (
    matrix(total_sum, nrow = length(candidates), ncol = ncol(feature_mat),
           byrow = TRUE) - subtree_sum[candidates, , drop = FALSE]
  ) / right_count
  mean_difference <- left_mean - right_mean
  gains <- left_count * right_count / length(nodes) *
    rowSums(mean_difference * mean_difference)
  selected <- order(
    -gains, subtree_min[candidates], nodes[candidates], nodes[parent[candidates]]
  )[1L]
  child <- candidates[[selected]]

  descendants <- integer(length(nodes))
  descendant_count <- 1L
  descendants[[descendant_count]] <- child
  cursor <- 1L
  while (cursor <= descendant_count) {
    node <- descendants[[cursor]]
    cursor <- cursor + 1L
    children <- adjacency[[node]][parent[adjacency[[node]]] == node]
    if (length(children)) {
      range <- descendant_count + seq_along(children)
      descendants[range] <- children
      descendant_count <- descendant_count + length(children)
    }
  }
  descendants <- descendants[seq_len(descendant_count)]
  list(
    cluster = as.integer(cluster_id),
    nodes = as.integer(nodes[descendants]),
    gain = as.numeric(gains[[selected]]),
    split_min = as.integer(subtree_min[[child]]),
    child = as.integer(nodes[[child]]),
    parent = as.integer(nodes[[parent[[child]]]])
  )
}

.exact_k_select_split <- function(labels, feature_mat, edges,
                                  min_cluster_size = 1L,
                                  K_target = max(labels) + 1L) {
  nodes_by_cluster <- split(seq_along(labels), labels)
  cluster_sizes <- tabulate(labels, nbins = max(labels))
  total_capacity <- sum(cluster_sizes %/% min_cluster_size)
  candidates <- lapply(seq_len(max(labels)), function(cluster_id) {
    nodes <- nodes_by_cluster[[as.character(cluster_id)]]
    if (is.null(nodes)) return(NULL)
    .exact_k_cluster_split_candidate(
      nodes, labels, cluster_id, feature_mat, edges, min_cluster_size,
      total_capacity - length(nodes) %/% min_cluster_size, K_target
    )
  })
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) return(NULL)
  gains <- vapply(candidates, `[[`, numeric(1L), "gain")
  clusters <- vapply(candidates, `[[`, integer(1L), "cluster")
  split_min <- vapply(candidates, `[[`, integer(1L), "split_min")
  children <- vapply(candidates, `[[`, integer(1L), "child")
  parents <- vapply(candidates, `[[`, integer(1L), "parent")
  candidates[[order(-gains, clusters, split_min, children, parents)[1L]]]
}

#' Split to a target while caching the best candidate per unchanged cluster.
#'
#' A split changes only its source cluster and the newly created cluster, so
#' rescanning every other cluster would repeat the same feature-tree work. The
#' cache is rebuilt at most once if the remaining minimum-size capacity becomes
#' exactly tight; before that transition every feasible tree cut can consume at
#' most one unit of capacity, and afterwards only capacity-preserving cuts are
#' eligible.
#'
#' @keywords internal
#' @noRd
.exact_k_split_to_target <- function(
    labels, feature_mat, K_target, edges, min_cluster_size,
    minimum_k, maximum_k, n_voxels, trace = NULL) {
  current_k <- max(labels)
  if (current_k >= K_target) return(labels)

  nodes_by_cluster <- split(seq_along(labels), labels)
  cluster_sizes <- tabulate(labels, nbins = current_k)
  total_capacity <- sum(cluster_sizes %/% min_cluster_size)

  candidate_for <- function(cluster_id) {
    nodes <- nodes_by_cluster[[cluster_id]]
    if (is.null(nodes)) return(NULL)
    .exact_k_cluster_split_candidate(
      nodes, labels, cluster_id, feature_mat, edges, min_cluster_size,
      total_capacity - length(nodes) %/% min_cluster_size, K_target
    )
  }
  rebuild_cache <- function() {
    lapply(seq_len(current_k), candidate_for)
  }
  candidate_cache <- rebuild_cache()

  while (current_k < K_target) {
    available <- which(!vapply(candidate_cache, is.null, logical(1L)))
    if (!length(available)) {
      .exact_k_abort(
        "no_connected_split", K_target, minimum_k, maximum_k,
        paste0(
          "no connected tree cut can create two clusters of size at least ",
          min_cluster_size, " while retaining capacity for requested K"
        ),
        n_voxels, min_cluster_size, current_k
      )
    }
    candidates <- candidate_cache[available]
    gains <- vapply(candidates, `[[`, numeric(1L), "gain")
    clusters <- vapply(candidates, `[[`, integer(1L), "cluster")
    split_min <- vapply(candidates, `[[`, integer(1L), "split_min")
    children <- vapply(candidates, `[[`, integer(1L), "child")
    parents <- vapply(candidates, `[[`, integer(1L), "parent")
    selected <- order(-gains, clusters, split_min, children, parents)[1L]
    split <- candidates[[selected]]

    source <- split$cluster
    source_nodes <- nodes_by_cluster[[source]]
    new_nodes <- sort(split$nodes)
    remainder_nodes <- setdiff(source_nodes, new_nodes)
    source_size <- length(source_nodes)
    current_k <- current_k + 1L
    .exact_k_trace_add(trace, list(
      type = "feature_tree_ward_split",
      source_label = as.integer(source),
      new_label = as.integer(current_k),
      source_size = as.integer(source_size),
      split_size = as.integer(length(new_nodes)),
      remainder_size = as.integer(length(remainder_nodes)),
      ward_gain = as.numeric(split$gain),
      split_min_voxel = as.integer(split$split_min),
      child_voxel = as.integer(split$child),
      parent_voxel = as.integer(split$parent)
    ))

    labels[new_nodes] <- current_k
    nodes_by_cluster[[source]] <- as.integer(remainder_nodes)
    nodes_by_cluster[[current_k]] <- as.integer(new_nodes)
    previous_capacity <- total_capacity
    total_capacity <- total_capacity - source_size %/% min_cluster_size +
      length(remainder_nodes) %/% min_cluster_size +
      length(new_nodes) %/% min_cluster_size

    if (previous_capacity > K_target && total_capacity == K_target) {
      candidate_cache <- rebuild_cache()
    } else {
      candidate_cache[source] <- list(candidate_for(source))
      candidate_cache[current_k] <- list(candidate_for(current_k))
    }
  }
  as.integer(labels)
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
#' @param min_cluster_size Positive minimum size of every output cluster.
#'   The compatibility default is one. A target that cannot satisfy this bound
#'   within each disconnected mask component fails with a structured condition.
#' @param record_operations Whether to attach an `exact_k_repair` metadata
#'   attribute containing pre/post sizes and deterministic repair operations.
#' @return Contiguous positive integer labels with exactly `K_target`
#'   connected components, or a `cluster4d_exact_k_infeasible` condition.
#' @keywords internal
force_exact_k <- function(labels, feature_mat, K_target, mask,
                          connectivity = 6L, graph_info = NULL,
                          min_cluster_size = 1L, record_operations = FALSE) {
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
  min_cluster_size <- .cluster4d_scalar_number(
    min_cluster_size, "min_cluster_size", "force_exact_k",
    lower = 1, upper = n_voxels, integer = TRUE
  )

  if (!is.logical(record_operations) || length(record_operations) != 1L ||
      is.na(record_operations)) {
    stop("force_exact_k: record_operations must be TRUE or FALSE", call. = FALSE)
  }

  trace <- .exact_k_new_trace(record_operations)
  input_labels <- .exact_k_relabel(labels)
  input_sizes <- tabulate(input_labels, nbins = max(input_labels))
  labels <- .exact_k_connected_labels(
    labels, graph_info$graph, graph_info$edges
  )
  normalized_k <- max(labels)
  normalized_sizes <- tabulate(labels, nbins = normalized_k)
  if (normalized_k != max(input_labels)) {
    .exact_k_trace_add(trace, list(
      type = "connectivity_normalization",
      input_k = as.integer(max(input_labels)),
      normalized_k = as.integer(normalized_k),
      input_sizes = as.integer(input_sizes),
      normalized_sizes = as.integer(normalized_sizes)
    ))
  }

  minimum_k <- length(unique(as.integer(graph_info$components)))
  component_sizes <- tabulate(as.integer(graph_info$components), nbins = minimum_k)
  maximum_k <- sum(component_sizes %/% min_cluster_size)
  if (K_target < minimum_k) {
    .exact_k_abort(
      "disconnected_mask_components", K_target, minimum_k, maximum_k,
      paste0(
        "target ", K_target, " is below the ", minimum_k,
        " disconnected mask components required by the topology"
      ),
      n_voxels, min_cluster_size, normalized_k
    )
  }
  if (any(component_sizes < min_cluster_size) || K_target > maximum_k) {
    .exact_k_abort(
      "minimum_cluster_size", K_target, minimum_k, maximum_k,
      paste0(
        "target ", K_target, " with min_cluster_size=", min_cluster_size,
        " is infeasible for component sizes ",
        paste(component_sizes, collapse = ",")
      ),
      n_voxels, min_cluster_size, normalized_k
    )
  }

  labels <- .exact_k_prepare_minimum_size(
    labels, feature_mat, graph_info$edges, K_target, min_cluster_size,
    minimum_k, maximum_k, n_voxels, trace
  )
  current_k <- max(labels)

  if (current_k > K_target) {
    labels <- .exact_k_merge_to_target(
      labels, feature_mat, K_target, graph_info$edges, minimum_k, n_voxels,
      maximum_k, min_cluster_size, trace
    )
    current_k <- max(labels)
  }

  if (current_k < K_target) {
    labels <- .exact_k_consolidate_split_capacity(
      labels, feature_mat, graph_info$edges, K_target, min_cluster_size,
      minimum_k, maximum_k, n_voxels, trace
    )
    current_k <- max(labels)
  }

  if (current_k < K_target) {
    labels <- .exact_k_split_to_target(
      labels, feature_mat, K_target, graph_info$edges, min_cluster_size,
      minimum_k, maximum_k, n_voxels, trace
    )
  }

  labels <- .exact_k_relabel(labels)
  sizes <- tabulate(labels, nbins = K_target)
  if (length(sizes) != K_target || any(sizes < min_cluster_size)) {
    .exact_k_abort(
      "postcondition_failed", K_target, minimum_k, maximum_k,
      paste0("output sizes were ", paste(sizes, collapse = ",")),
      n_voxels, min_cluster_size, max(labels)
    )
  }
  if (isTRUE(record_operations)) {
    operation_types <- vapply(
      trace$operations, function(operation) operation$type, character(1L)
    )
    has_merge <- any(grepl("merge$", operation_types))
    has_split <- any(grepl("split$", operation_types))
    direction <- if (has_merge && has_split) {
      "mixed"
    } else if (has_merge) {
      "merge"
    } else if (has_split) {
      "split"
    } else {
      "none"
    }
    attr(labels, "exact_k_repair") <- list(
      requested_k = as.integer(K_target),
      min_cluster_size = as.integer(min_cluster_size),
      connectivity = as.integer(connectivity),
      input_k = as.integer(max(input_labels)),
      input_sizes = as.integer(input_sizes),
      normalized_k = as.integer(normalized_k),
      normalized_sizes = as.integer(normalized_sizes),
      final_k = as.integer(K_target),
      final_sizes = as.integer(sizes),
      direction = direction,
      operations = trace$operations
    )
  }
  labels
}
