exact_k_mask <- function(values, spacing = c(1, 1, 1)) {
  neuroim2::NeuroVol(
    values,
    neuroim2::NeuroSpace(dim(values), spacing = spacing, origin = c(0, 0, 0))
  )
}

exact_k_offsets <- function(connectivity) {
  offsets <- as.matrix(expand.grid(-1:1, -1:1, -1:1))
  offsets <- offsets[rowSums(abs(offsets)) > 0, , drop = FALSE]
  if (connectivity == 6L) {
    offsets[rowSums(abs(offsets)) == 1L, , drop = FALSE]
  } else if (connectivity == 18L) {
    offsets[rowSums(abs(offsets)) <= 2L, , drop = FALSE]
  } else offsets
}

independent_label_components <- function(labels, mask, connectivity) {
  mask_idx <- which(as.array(mask) > 0)
  coords <- arrayInd(mask_idx, dim(mask))
  offsets <- exact_k_offsets(connectivity)
  key <- function(x) paste(x[, 1], x[, 2], x[, 3], sep = ":")
  lookup <- setNames(seq_len(nrow(coords)), key(coords))
  component_count <- 0L
  visited <- logical(length(labels))

  for (start in seq_along(labels)) {
    if (visited[start]) next
    component_count <- component_count + 1L
    queue <- start
    visited[start] <- TRUE
    while (length(queue)) {
      node <- queue[[1L]]
      queue <- queue[-1L]
      neighbor_coords <- sweep(offsets, 2, coords[node, ], "+")
      inside <- apply(
        neighbor_coords, 1,
        function(x) all(x >= 1L & x <= dim(mask))
      )
      candidates <- unname(lookup[key(neighbor_coords[inside, , drop = FALSE])])
      candidates <- candidates[
        !is.na(candidates) & !visited[candidates] & labels[candidates] == labels[node]
      ]
      if (length(candidates)) {
        visited[candidates] <- TRUE
        queue <- c(queue, candidates)
      }
    }
  }
  component_count
}

expect_exact_k_topology <- function(labels, mask, target, connectivity) {
  expect_identical(sort(unique(labels)), seq_len(target))
  expect_identical(
    independent_label_components(labels, mask, connectivity),
    target
  )
}

independent_grid_edges <- function(mask, connectivity) {
  coords <- arrayInd(which(as.array(mask) > 0), dim(mask))
  if (nrow(coords) < 2L) return(matrix(integer(), ncol = 2L))
  pairs <- utils::combn(seq_len(nrow(coords)), 2L)
  delta <- abs(coords[pairs[1L, ], , drop = FALSE] -
                 coords[pairs[2L, ], , drop = FALSE])
  chebyshev <- apply(delta, 1L, max)
  manhattan <- rowSums(delta)
  adjacent <- chebyshev == 1L & if (connectivity == 6L) {
    manhattan == 1L
  } else if (connectivity == 18L) {
    manhattan <= 2L
  } else {
    TRUE
  }
  t(pairs[, adjacent, drop = FALSE])
}

slow_tree_split_oracle <- function(labels, features, mask, connectivity,
                                   min_cluster_size, target_k) {
  edges <- independent_grid_edges(mask, connectivity)
  sizes <- tabulate(labels, nbins = max(labels))
  total_capacity <- sum(sizes %/% min_cluster_size)
  candidates <- list()

  for (cluster_id in seq_len(max(labels))) {
    nodes <- which(labels == cluster_id)
    induced <- edges[
      labels[edges[, 1L]] == cluster_id & labels[edges[, 2L]] == cluster_id,
      , drop = FALSE
    ]
    if (length(nodes) < 2L * min_cluster_size) next
    stopifnot(nrow(induced) == length(nodes) - 1L)
    root <- min(nodes)

    for (cut in seq_len(nrow(induced))) {
      remaining <- induced[-cut, , drop = FALSE]
      visited <- induced[cut, 1L]
      frontier <- visited
      while (length(frontier)) {
        node <- frontier[[1L]]
        frontier <- frontier[-1L]
        incident <- remaining[
          remaining[, 1L] == node | remaining[, 2L] == node, , drop = FALSE
        ]
        neighbors <- setdiff(as.integer(incident), node)
        neighbors <- setdiff(neighbors, visited)
        if (length(neighbors)) {
          visited <- c(visited, neighbors)
          frontier <- c(frontier, neighbors)
        }
      }
      split_nodes <- sort(unique(visited))
      if (root %in% split_nodes) split_nodes <- setdiff(nodes, split_nodes)
      remainder <- setdiff(nodes, split_nodes)
      if (length(split_nodes) < min_cluster_size ||
          length(remainder) < min_cluster_size) next
      capacity_after <- total_capacity - length(nodes) %/% min_cluster_size +
        length(split_nodes) %/% min_cluster_size +
        length(remainder) %/% min_cluster_size
      if (capacity_after < target_k) next

      split_mean <- colMeans(features[split_nodes, , drop = FALSE])
      remainder_mean <- colMeans(features[remainder, , drop = FALSE])
      gain <- length(split_nodes) * length(remainder) / length(nodes) *
        sum((split_mean - remainder_mean)^2)
      cut_edge <- induced[cut, ]
      child <- cut_edge[cut_edge %in% split_nodes]
      parent <- cut_edge[!cut_edge %in% split_nodes]
      candidates[[length(candidates) + 1L]] <- list(
        cluster = as.integer(cluster_id),
        nodes = as.integer(split_nodes),
        gain = as.numeric(gain),
        split_min = as.integer(min(split_nodes)),
        child = as.integer(child),
        parent = as.integer(parent)
      )
    }
  }
  if (!length(candidates)) return(NULL)
  gains <- vapply(candidates, `[[`, numeric(1L), "gain")
  clusters <- vapply(candidates, `[[`, integer(1L), "cluster")
  split_min <- vapply(candidates, `[[`, integer(1L), "split_min")
  children <- vapply(candidates, `[[`, integer(1L), "child")
  parents <- vapply(candidates, `[[`, integer(1L), "parent")
  candidates[[order(-gains, clusters, split_min, children, parents)[1L]]]
}

test_that("feature-tree splits match an independent exhaustive tree-cut oracle", {
  line <- exact_k_mask(array(1, c(6L, 1L, 1L)))
  branch_values <- array(0, c(3L, 3L, 1L))
  branch_values[cbind(c(1L, 2L, 2L, 2L, 3L),
                      c(2L, 1L, 2L, 3L, 2L), 1L)] <- 1L
  branch <- exact_k_mask(branch_values)
  fixtures <- list(
    list(mask = line, labels = rep(1L, 6L), levels = c(-1, 0, 1),
         minimum = c(1L, 2L), target = 2L),
    list(mask = branch, labels = rep(1L, 5L), levels = c(0, 1),
         minimum = c(1L, 2L), target = 2L),
    list(mask = line, labels = c(1L, 1L, 1L, 2L, 2L, 2L), levels = c(0, 1),
         minimum = 1L, target = 3L)
  )

  for (fixture in fixtures) {
    assignments <- as.matrix(expand.grid(
      rep(list(fixture$levels), length(fixture$labels))
    ))
    graph <- .exact_k_graph(fixture$mask, 6L)
    for (minimum in fixture$minimum) {
      for (row in seq_len(nrow(assignments))) {
        features <- matrix(assignments[row, ], ncol = 1L)
        expected <- slow_tree_split_oracle(
          fixture$labels, features, fixture$mask, 6L,
          minimum, fixture$target
        )
        actual <- .exact_k_select_split(
          fixture$labels, features, graph$edges, minimum, fixture$target
        )
        expect_equal(actual$gain, expected$gain, tolerance = 1e-12)
        expect_identical(actual$cluster, expected$cluster)
        expect_identical(sort(actual$nodes), sort(expected$nodes))
        expect_identical(
          actual[c("split_min", "child", "parent")],
          expected[c("split_min", "child", "parent")]
        )
      }
    }
  }
})

test_that("minimum-size feasibility is explicit and preserved through awkward repairs", {
  mask <- exact_k_mask(array(1, c(12L, 1L, 1L)))
  features <- cbind(seq_len(12L), rep(c(0, 1), 6L))

  split <- force_exact_k(
    rep(1L, 12L), features, 3L, mask, 6L, min_cluster_size = 3L
  )
  expect_exact_k_topology(split, mask, 3L, 6L)
  expect_true(all(tabulate(split) >= 3L))

  awkward <- force_exact_k(
    c(1L, 2L, 2L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L),
    features, 3L, mask, 6L, min_cluster_size = 3L
  )
  expect_exact_k_topology(awkward, mask, 3L, 6L)
  expect_true(all(tabulate(awkward) >= 3L))

  impossible <- tryCatch(
    force_exact_k(
      rep(1L, 12L), features, 5L, mask, 6L, min_cluster_size = 3L
    ),
    cluster4d_exact_k_infeasible = identity
  )
  expect_s3_class(impossible, "cluster4d_exact_k_infeasible")
  expect_identical(impossible$reason, "minimum_cluster_size")
  expect_identical(impossible$target_k, 5L)
  expect_identical(impossible$n_voxels, 12L)
  expect_identical(impossible$min_cluster_size, 3L)
  expect_identical(impossible$current_k, 1L)

  disconnected_values <- array(0, c(9L, 1L, 1L))
  disconnected_values[c(1L, 2L, 5L:9L)] <- 1L
  disconnected <- exact_k_mask(disconnected_values)
  graph_impossible <- tryCatch(
    force_exact_k(
      rep(1L, 7L), matrix(seq_len(7L), ncol = 1L), 2L,
      disconnected, 6L, min_cluster_size = 3L
    ),
    cluster4d_exact_k_infeasible = identity
  )
  expect_s3_class(graph_impossible, "cluster4d_exact_k_infeasible")
  expect_identical(graph_impossible$reason, "minimum_cluster_size")
  expect_identical(graph_impossible$target_k, 2L)
  expect_identical(graph_impossible$n_voxels, 7L)
  expect_identical(graph_impossible$min_cluster_size, 3L)
  expect_identical(graph_impossible$current_k, 2L)
})

test_that("random connected masks satisfy exact K, size, and connectivity contracts", {
  set.seed(4831)
  dims <- c(6L, 6L, 4L)
  for (trial in seq_len(40L)) {
    connectivity <- sample(c(6L, 18L, 26L), 1L)
    offsets <- exact_k_offsets(connectivity)
    wanted <- sample(15:32, 1L)
    occupied <- matrix(c(3L, 3L, 2L), nrow = 1L)
    while (nrow(occupied) < wanted) {
      origin <- occupied[sample(seq_len(nrow(occupied)), 1L), ]
      candidate <- origin + offsets[sample(seq_len(nrow(offsets)), 1L), ]
      if (all(candidate >= 1L & candidate <= dims) &&
          !any(apply(occupied, 1L, function(x) all(x == candidate)))) {
        occupied <- rbind(occupied, candidate)
      }
    }
    values <- array(0L, dims)
    values[occupied] <- 1L
    mask <- exact_k_mask(values)
    n <- sum(values)
    features <- matrix(rnorm(n * 3L), ncol = 3L)
    initial <- sample(seq_len(min(6L, n)), n, replace = TRUE)
    minimum <- sample(1:3, 1L)
    target <- sample(seq_len(min(5L, n %/% minimum)), 1L)

    result <- force_exact_k(
      initial, features, target, mask, connectivity,
      min_cluster_size = minimum
    )
    expect_exact_k_topology(result, mask, target, connectivity)
    expect_true(all(tabulate(result) >= minimum), info = paste("trial", trial))
  }
})

test_that("small-grid targets are exhaustively feasible exactly when topology permits", {
  dims <- c(2L, 2L, 1L)
  for (bits in seq_len(2^prod(dims) - 1L)) {
    included <- as.logical(intToBits(bits)[seq_len(prod(dims))])
    mask <- exact_k_mask(array(as.integer(included), dims))
    n <- sum(included)
    features <- matrix(seq_len(n), ncol = 1L)
    initial <- rep(1L, n)

    for (connectivity in c(6L, 18L, 26L)) {
      minimum_k <- independent_label_components(
        rep(1L, n), mask, connectivity
      )
      # Count components independently by asking whether singleton IDs can be
      # merged to a target; the smallest successful target is the oracle min.
      feasible <- integer()
      for (target in seq_len(n)) {
        attempt <- tryCatch(
          force_exact_k(initial, features, target, mask, connectivity),
          cluster4d_exact_k_infeasible = identity
        )
        if (inherits(attempt, "cluster4d_exact_k_infeasible")) {
          expect_identical(attempt$reason, "disconnected_mask_components")
        } else {
          feasible <- c(feasible, target)
          expect_exact_k_topology(attempt, mask, target, connectivity)
        }
      }
      expect_identical(feasible, seq.int(min(feasible), n))
      expect_identical(min(feasible), minimum_k)
    }
  }
})

test_that("all binary initial partitions on a small grid retain connected labels", {
  mask <- exact_k_mask(array(1, c(2, 2, 1)))
  features <- matrix(c(0, 1, 4, 9), ncol = 1L)
  for (bits in 0:(2^4 - 1L)) {
    initial <- as.integer(as.logical(intToBits(bits)[1:4])) + 1L
    for (connectivity in c(6L, 18L, 26L)) {
      for (target in 1:4) {
        result <- force_exact_k(initial, features, target, mask, connectivity)
        expect_exact_k_topology(result, mask, target, connectivity)
      }
    }
  }
})

test_that("gap, thin-mask, disconnected, under-target, over-target, and K=N cases are safe", {
  line_mask <- exact_k_mask(array(1, c(7, 1, 1)))
  line_features <- matrix(c(0, 8, 0, 8, 0, 8, 0), ncol = 1L)

  under <- force_exact_k(rep(1L, 7), line_features, 4, line_mask, 6)
  expect_exact_k_topology(under, line_mask, 4L, 6L)

  over <- force_exact_k(seq_len(7), line_features, 2, line_mask, 6)
  expect_exact_k_topology(over, line_mask, 2L, 6L)

  singleton <- force_exact_k(rep(1L, 7), line_features, 7, line_mask, 6)
  expect_exact_k_topology(singleton, line_mask, 7L, 6L)

  thin_values <- array(0, c(5, 3, 1))
  thin_values[, 2, 1] <- 1
  thin <- exact_k_mask(thin_values)
  thin_labels <- force_exact_k(rep(1L, 5), matrix(1:5, ncol = 1), 3, thin, 6)
  expect_exact_k_topology(thin_labels, thin, 3L, 6L)

  gap_values <- array(0, c(5, 1, 1))
  gap_values[c(1, 2, 4, 5), 1, 1] <- 1
  gap <- exact_k_mask(gap_values)
  condition <- tryCatch(
    force_exact_k(c(1, 1, 2, 2), matrix(1:4, ncol = 1), 1, gap, 6),
    cluster4d_exact_k_infeasible = identity
  )
  expect_s3_class(condition, "cluster4d_exact_k_infeasible")
  expect_identical(condition$reason, "disconnected_mask_components")
  expect_identical(condition$minimum_k, 2L)

  gap_result <- force_exact_k(
    c(1, 1, 1, 1), matrix(1:4, ncol = 1), 3, gap, 6
  )
  expect_exact_k_topology(gap_result, gap, 3L, 6L)
})

test_that("6, 18, and 26 connectivity have distinct diagonal feasibility", {
  edge_values <- array(0, c(2, 2, 1))
  edge_values[c(1, 4)] <- 1
  edge_mask <- exact_k_mask(edge_values)
  edge_features <- matrix(c(0, 1), ncol = 1)
  edge_condition <- tryCatch(
    force_exact_k(c(1, 1), edge_features, 1, edge_mask, 6),
    cluster4d_exact_k_infeasible = identity
  )
  expect_s3_class(edge_condition, "cluster4d_exact_k_infeasible")
  expect_exact_k_topology(
    force_exact_k(c(1, 1), edge_features, 1, edge_mask, 18),
    edge_mask, 1L, 18L
  )

  corner_values <- array(0, c(2, 2, 2))
  corner_values[c(1, 8)] <- 1
  corner_mask <- exact_k_mask(corner_values)
  corner_condition <- tryCatch(
    force_exact_k(c(1, 1), edge_features, 1, corner_mask, 18),
    cluster4d_exact_k_infeasible = identity
  )
  expect_s3_class(corner_condition, "cluster4d_exact_k_infeasible")
  expect_exact_k_topology(
    force_exact_k(c(1, 1), edge_features, 1, corner_mask, 26),
    corner_mask, 1L, 26L
  )
})

slow_reference_merge <- function(labels, features, edges) {
  candidates <- unique(t(apply(edges, 1, function(edge) {
    sort(labels[edge])
  })))
  candidates <- candidates[candidates[, 1] != candidates[, 2], , drop = FALSE]
  candidates <- unique(candidates[order(candidates[, 1], candidates[, 2]), , drop = FALSE])
  costs <- apply(candidates, 1, function(pair) {
    left <- features[labels == pair[1], , drop = FALSE]
    right <- features[labels == pair[2], , drop = FALSE]
    nrow(left) * nrow(right) / (nrow(left) + nrow(right)) *
      sum((colMeans(left) - colMeans(right))^2)
  })
  selected <- order(costs, candidates[, 1], candidates[, 2])[1]
  list(
    left = as.integer(candidates[selected, 1]),
    right = as.integer(candidates[selected, 2]),
    cost = as.numeric(costs[selected])
  )
}

test_that("adjacent Ward merge costs and selected operations match a slow reference", {
  mask <- exact_k_mask(array(1, c(6, 1, 1)))
  graph <- .exact_k_graph(mask, 6)
  labels <- c(1L, 1L, 2L, 3L, 3L, 4L)
  features <- cbind(c(0, 0.2, 1, 4, 4.2, 9), c(1, 1, 2, 0, 0, 3))
  fast <- .exact_k_select_merge(labels, features, graph$edges)
  slow <- slow_reference_merge(labels, features, graph$edges)
  expect_identical(fast[c("left", "right")], slow[c("left", "right")])
  expect_equal(fast$cost, slow$cost, tolerance = 1e-12)
})

test_that("exact-K repair is deterministic and does not touch RNG state", {
  mask <- exact_k_mask(array(1, c(8, 1, 1)))
  features <- matrix(c(0, 9, 0, 8, 1, 7, 2, 6), ncol = 1L)
  set.seed(42)
  before <- .Random.seed
  first <- force_exact_k(rep(1L, 8), features, 5, mask, 6)
  after <- .Random.seed
  second <- force_exact_k(rep(1L, 8), features, 5, mask, 6)
  expect_identical(first, second)
  expect_identical(after, before)
})

slow_merge_to_target <- function(labels, features, target, edges) {
  while (max(labels) > target) {
    merge <- .exact_k_select_merge(labels, features, edges)
    labels[labels == merge$right] <- merge$left
    labels <- .exact_k_relabel(labels)
  }
  labels
}

slow_corrected_force_exact_k <- function(labels, features, target, mask,
                                         connectivity,
                                         min_cluster_size = 1L) {
  graph <- .exact_k_graph(mask, connectivity)
  n <- length(labels)
  minimum <- length(unique(graph$components))
  component_sizes <- tabulate(graph$components, nbins = minimum)
  maximum <- sum(component_sizes %/% min_cluster_size)
  labels <- .exact_k_connected_labels(labels, graph$graph, graph$edges)
  labels <- .exact_k_prepare_minimum_size(
    labels, features, graph$edges, target, min_cluster_size,
    minimum, maximum, n
  )
  if (max(labels) > target) {
    labels <- .exact_k_merge_to_target(
      labels, features, target, graph$edges, minimum, n,
      maximum, min_cluster_size
    )
  }
  if (max(labels) < target) {
    labels <- .exact_k_consolidate_split_capacity(
      labels, features, graph$edges, target, min_cluster_size,
      minimum, maximum, n
    )
  }
  while (max(labels) < target) {
    split <- .exact_k_select_split(
      labels, features, graph$edges, min_cluster_size, target
    )
    stopifnot(!is.null(split))
    labels[split$nodes] <- max(labels) + 1L
  }
  .exact_k_relabel(labels)
}

test_that("cached exact-K splits are bit-identical to the corrected slow oracle", {
  set.seed(61227)
  for (trial in seq_len(60L)) {
    dims <- sample(2:5, 3L, replace = TRUE)
    included <- sample(
      c(FALSE, TRUE), prod(dims), replace = TRUE,
      prob = c(0.25, 0.75)
    )
    if (sum(included) < 3L) included[seq_len(3L)] <- TRUE
    mask <- exact_k_mask(array(as.integer(included), dims))
    connectivity <- sample(c(6L, 18L, 26L), 1L)
    graph <- .exact_k_graph(mask, connectivity)
    n <- length(graph$mask_idx)
    features <- matrix(rnorm(n * 4L), ncol = 4L)
    initial <- sample(seq_len(min(n, 7L)), n, replace = TRUE)
    normalized <- .exact_k_connected_labels(initial, graph$graph, graph$edges)
    minimum <- length(unique(graph$components))
    current <- max(normalized)
    if (trial %% 2L && current < n) {
      target <- sample(seq.int(current, min(n, current + 6L)), 1L)
    } else {
      target <- sample(seq.int(minimum, current), 1L)
    }

    expected <- slow_corrected_force_exact_k(
      initial, features, target, mask, connectivity
    )
    actual <- force_exact_k(
      initial, features, target, mask, connectivity,
      min_cluster_size = 1L
    )
    expect_identical(actual, expected, info = paste("trial", trial))
  }

  for (trial in seq_len(30L)) {
    n <- sample(12:30, 1L)
    mask <- exact_k_mask(array(1, c(n, 1L, 1L)))
    features <- matrix(rnorm(n * 3L), ncol = 3L)
    minimum <- sample(2:4, 1L)
    target <- sample(seq_len(n %/% minimum), 1L)
    initial <- rep(1L, n)
    expected <- slow_corrected_force_exact_k(
      initial, features, target, mask, 6L, minimum
    )
    actual <- force_exact_k(
      initial, features, target, mask, 6L,
      min_cluster_size = minimum
    )
    expect_identical(actual, expected, info = paste("minimum-size trial", trial))
  }
})

test_that("incremental exact-K merging matches the rebuilding oracle", {
  set.seed(2112)
  for (trial in seq_len(60L)) {
    dims <- sample(2:5, 3L, replace = TRUE)
    included <- sample(c(FALSE, TRUE), prod(dims), replace = TRUE,
                       prob = c(0.35, 0.65))
    if (sum(included) < 2L) included[seq_len(2L)] <- TRUE
    mask <- exact_k_mask(array(as.integer(included), dims))
    connectivity <- sample(c(6L, 18L, 26L), 1L)
    graph <- .exact_k_graph(mask, connectivity)
    n <- length(graph$mask_idx)
    features <- matrix(rnorm(n * 4L), n, 4L)
    initial <- sample(seq_len(min(n, 8L)), n, replace = TRUE)
    labels <- .exact_k_connected_labels(initial, graph$graph, graph$edges)
    minimum <- length(unique(graph$components))
    target <- sample(seq.int(minimum, max(labels)), 1L)

    expected <- slow_merge_to_target(labels, features, target, graph$edges)
    actual <- .exact_k_merge_to_target(
      labels, features, target, graph$edges, minimum, n
    )
    expect_identical(actual, expected, info = paste("trial", trial))
  }
})

test_that("prebuilt exact-K graphs are bound to mask identity", {
  dims <- c(4L, 1L, 1L)
  space <- neuroim2::NeuroSpace(dims)
  connected <- exact_k_mask(array(c(1, 1, 0, 0), dims))
  disconnected <- neuroim2::NeuroVol(array(c(1, 0, 0, 1), dims), space)
  receipt <- .exact_k_graph(connected, 6L)

  expect_error(
    force_exact_k(
      c(1L, 1L), matrix(c(0, 1), ncol = 1L), 1L,
      disconnected, 6L, graph_info = receipt
    ),
    "graph_info mask does not match"
  )

  reshaped <- neuroim2::NeuroVol(
    array(c(1, 1, 0, 0), c(2L, 2L, 1L)),
    neuroim2::NeuroSpace(c(2L, 2L, 1L))
  )
  expect_error(
    force_exact_k(
      c(1L, 1L), matrix(c(0, 1), ncol = 1L), 1L,
      reshaped, 6L, graph_info = receipt
    ),
    "graph_info dimensions do not match"
  )
})

test_that("ACSC, MCL, and ReNA exact-K callers preserve topology", {
  set.seed(7)
  dims <- c(4L, 4L, 2L)
  sp <- neuroim2::NeuroSpace(dims)
  mask <- neuroim2::NeuroVol(array(1, dims), sp)
  vec <- do.call(neuroim2::concat, lapply(1:6, function(i) {
    neuroim2::NeuroVol(array(stats::rnorm(prod(dims)), dims), sp)
  }))

  acsc_result <- suppressWarnings(acsc(
    vec, mask, K = 4, block_size = 2, refine = FALSE
  ))
  expect_exact_k_topology(acsc_result$cluster, mask, 4L, 6L)

  mcl_result <- suppressWarnings(cluster4d_mcl(
    vec, mask, n_clusters = 4, max_iterations = 2,
    connectivity = 6, exact_k = TRUE
  ))
  expect_exact_k_topology(mcl_result$labels, mask, 4L, 6L)

  rena_result <- suppressWarnings(cluster4d_rena(
    vec, mask, n_clusters = 4, max_iterations = 3,
    connectivity = 6, exact_k = TRUE
  ))
  expect_exact_k_topology(rena_result$labels, mask, 4L, 6L)
})
