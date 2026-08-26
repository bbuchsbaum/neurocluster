topology_oracle_edges <- function(mask_values, connectivity) {
  dims <- dim(mask_values)
  retained <- which(mask_values > 0L)
  if (length(retained) < 2L) return(matrix(integer(), ncol = 2L))
  grid <- arrayInd(retained, dims)
  pairs <- utils::combn(seq_along(retained), 2L)
  delta <- abs(grid[pairs[1L, ], , drop = FALSE] -
                 grid[pairs[2L, ], , drop = FALSE])
  manhattan <- rowSums(delta)
  keep <- apply(delta, 1L, max) <= 1L & switch(
    as.character(connectivity),
    `6` = manhattan == 1L,
    `18` = manhattan >= 1L & manhattan <= 2L,
    `26` = manhattan >= 1L
  )
  if (!any(keep)) return(matrix(integer(), ncol = 2L))
  t(pairs[, keep, drop = FALSE])
}

topology_oracle_component_count <- function(mask_values, connectivity) {
  n_voxels <- sum(mask_values > 0L)
  if (n_voxels == 0L) return(0L)

  edges <- topology_oracle_edges(mask_values, connectivity)
  graph <- igraph::make_empty_graph(n_voxels, directed = FALSE)
  if (nrow(edges)) {
    graph <- igraph::add_edges(graph, as.vector(t(edges)))
  }
  as.integer(igraph::components(graph)$no)
}

topology_native_edges <- function(mask, connectivity) {
  idx <- which(as.array(mask) > 0)
  adjacency <- neurocluster:::build_grid_adjacency(mask, idx, connectivity)
  summary <- Matrix::summary(adjacency)
  keep <- summary$i < summary$j
  edges <- cbind(summary$i[keep], summary$j[keep])
  if (!nrow(edges)) matrix(integer(), ncol = 2L) else
    edges[order(edges[, 1L], edges[, 2L]), , drop = FALSE]
}

topology_labels_connected <- function(labels, mask, connectivity) {
  graph_info <- neurocluster:::.exact_k_graph(mask, connectivity)
  normalized <- neurocluster:::.exact_k_connected_labels(
    labels, graph_info$graph, graph_info$edges
  )
  length(unique(normalized)) == length(unique(labels))
}

topology_slow_dijkstra <- function(feature_mat, seeds, neighbors, distances,
                                   alpha, compactness) {
  n <- ncol(feature_mat)
  best_cost <- rep(Inf, n)
  best_label <- integer(n)
  final <- logical(n)
  labels <- integer(n)
  best_cost[seeds] <- 0
  best_label[seeds] <- seq_along(seeds)

  feature_distance <- function(i, j) {
    if (nrow(feature_mat) > 1L) {
      dot <- sum(feature_mat[, i] * feature_mat[, j])
      return(1 - max(-1, min(1, dot)))
    }
    abs(feature_mat[1L, i] - feature_mat[1L, j])
  }

  repeat {
    candidates <- which(!final & is.finite(best_cost))
    if (!length(candidates)) break
    selected <- candidates[order(
      best_cost[candidates], best_label[candidates], candidates
    )[1L]]
    final[selected] <- TRUE
    labels[selected] <- best_label[selected]
    for (slot in seq_len(ncol(neighbors))) {
      neighbor <- neighbors[selected, slot]
      if (neighbor <= 0L || final[neighbor]) next
      edge <- alpha * feature_distance(selected, neighbor) +
        (1 - alpha) * distances[selected, slot] / compactness
      candidate <- best_cost[selected] + edge
      label <- labels[selected]
      if (candidate < best_cost[neighbor] - 1e-12 ||
          (abs(candidate - best_cost[neighbor]) <= 1e-12 &&
           label < best_label[neighbor])) {
        best_cost[neighbor] <- candidate
        best_label[neighbor] <- label
      }
    }
  }
  labels
}

topology_fixture <- function(spacing = c(1, 1, 1), origin = c(0, 0, 0)) {
  dims <- c(5L, 4L, 1L)
  n_time <- 8L
  set.seed(826)
  values <- array(rnorm(prod(dims) * n_time, sd = 0.15), c(dims, n_time))
  pattern_a <- sin(seq(0, 2 * pi, length.out = n_time))
  pattern_b <- cos(seq(0, 2 * pi, length.out = n_time))
  values[1:2, , 1, ] <- values[1:2, , 1, ] +
    array(rep(pattern_a, each = 8L), c(2L, 4L, n_time))
  values[3:5, , 1, ] <- values[3:5, , 1, ] +
    array(rep(pattern_b, each = 12L), c(3L, 4L, n_time))
  space <- neuroim2::NeuroSpace(dims, spacing = spacing, origin = origin)
  frames <- lapply(seq_len(n_time), function(time) {
    neuroim2::NeuroVol(
      array(values[, , , time], dims), space
    )
  })
  list(
    vec = do.call(neuroim2::concat, frames),
    mask = neuroim2::NeuroVol(array(1, dims), space),
    values = values,
    dims = dims
  )
}

test_that("masked grid adjacency matches an exhaustive independent oracle", {
  dims <- c(2L, 2L, 2L)
  for (bits in 1:(2^prod(dims) - 1L)) {
    included <- as.integer(intToBits(bits))[seq_len(prod(dims))]
    mask_values <- array(included, dims)
    mask <- neuroim2::NeuroVol(
      mask_values, neuroim2::NeuroSpace(dims)
    )
    for (connectivity in c(6L, 18L, 26L)) {
      graph_info <- neurocluster:::.exact_k_graph(mask, connectivity)
      native_edges <- graph_info$edges
      if (nrow(native_edges)) {
        native_edges <- native_edges[
          order(native_edges[, 1L], native_edges[, 2L]), , drop = FALSE
        ]
      }
      expect_equal(
        topology_native_edges(mask, connectivity),
        topology_oracle_edges(mask_values, connectivity),
        info = paste("mask", bits, "connectivity", connectivity)
      )
      expect_equal(
        native_edges,
        topology_oracle_edges(mask_values, connectivity),
        info = paste("direct graph mask", bits, "connectivity", connectivity)
      )
      expect_identical(
        length(unique(graph_info$components)),
        topology_oracle_component_count(mask_values, connectivity)
      )
    }
  }

  gap <- neuroim2::NeuroVol(
    array(c(1, 0, 1), c(3, 1, 1)),
    neuroim2::NeuroSpace(c(3, 1, 1))
  )
  expect_equal(sum(neurocluster:::rena_build_connectivity(
    gap, which(as.array(gap) > 0), 26
  )), 0)
})

test_that("direct grid graph crosses the former 65536-voxel boundary", {
  for (n in c(65535L, 65536L)) {
    graph <- build_grid_edges_cpp(seq_len(n), c(n, 1L, 1L), 6L)
    expect_identical(nrow(graph$edges), n - 1L, info = paste("n", n))
    expect_identical(graph$neighbor_ptr[[n + 1L]], 2L * (n - 1L))
    expect_identical(unique(graph$components), 1L)
  }
})

test_that("G3S native propagation equals a slow multi-source Dijkstra oracle", {
  dims <- c(3L, 2L, 1L)
  mask <- neuroim2::NeuroVol(array(1, dims), neuroim2::NeuroSpace(dims))
  graph <- neurocluster:::build_grid_neighbors_g3s(
    mask, seq_len(prod(dims)), 6L
  )
  feature_mat <- rbind(
    c(1, 0.95, 0.7, 0.2, 0.05, 0),
    c(0, 0.31, 0.71, 0.98, 1, 1)
  )
  feature_mat <- apply(feature_mat, 2L, function(x) x / sqrt(sum(x^2)))
  seeds <- c(1L, 6L)

  for (alpha in c(0, 0.35, 1)) {
    expected <- topology_slow_dijkstra(
      feature_mat, seeds, graph$nn.index, graph$nn.dist,
      alpha = alpha, compactness = 1.7
    )
    actual <- g3s_propagate_cpp(
      feature_mat, seeds, graph$nn.index, graph$nn.dist,
      alpha = alpha, compactness = 1.7
    )
    expect_identical(as.integer(actual), as.integer(expected))
  }
})

test_that("G3S compression expands rank until its retained-variance contract holds", {
  set.seed(2086)
  feature_mat <- matrix(rnorm(48L * 24L), 48L, 24L)
  requested <- 3L
  threshold <- 0.9

  compressed <- suppressWarnings(compress_features_svd(
    feature_mat,
    n_components = requested,
    variance_threshold = threshold,
    use_irlba = FALSE,
    use_rsvd = FALSE
  ))

  scaled <- base::scale(feature_mat)
  singular_values <- base::svd(scaled, nu = 0L, nv = 0L)$d
  oracle_ratio <- sum(singular_values[seq_len(compressed$n_components)]^2) /
    sum(singular_values^2)

  expect_gt(compressed$n_components, requested)
  expect_equal(ncol(compressed$features), compressed$n_components)
  expect_gte(compressed$variance_explained, threshold - 1e-12)
  expect_equal(compressed$variance_explained, oracle_ratio, tolerance = 1e-12)

  full_rank <- suppressWarnings(compress_features_svd(
    feature_mat,
    n_components = 1L,
    variance_threshold = 1,
    use_irlba = FALSE,
    use_rsvd = FALSE
  ))
  full_oracle <- cumsum(singular_values^2) / sum(singular_values^2)
  expected_full_rank <- which(
    full_oracle + 100 * .Machine$double.eps >= 1
  )[[1L]]
  expect_equal(full_rank$n_components, expected_full_rank)
  expect_equal(full_rank$variance_explained, 1, tolerance = 1e-12)

  zero_variance <- suppressWarnings(compress_features_svd(
    matrix(4, 8L, 5L),
    n_components = 2L,
    variance_threshold = 0.95,
    use_irlba = FALSE,
    use_rsvd = FALSE
  ))
  expect_equal(zero_variance$variance_explained, 1)
  expect_true(all(is.finite(zero_variance$features)))
})

test_that("randomized G3S compression is deterministic and RNG-neutral", {
  skip_if_not_installed("rsvd")
  set.seed(2087)
  feature_mat <- matrix(rnorm(10001L * 12L), 10001L, 12L)
  before <- .Random.seed
  first <- suppressMessages(compress_features_svd(
    feature_mat, n_components = 4L, variance_threshold = 0,
    use_irlba = FALSE, use_rsvd = TRUE
  ))
  expect_identical(.Random.seed, before)
  second <- suppressMessages(compress_features_svd(
    feature_mat, n_components = 4L, variance_threshold = 0,
    use_irlba = FALSE, use_rsvd = TRUE
  ))
  expect_identical(second$features, first$features)
  expect_identical(.Random.seed, before)

  saved <- .Random.seed
  on.exit(assign(".Random.seed", saved, envir = .GlobalEnv), add = TRUE)
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(suppressMessages(compress_features_svd(
    feature_mat, n_components = 4L, variance_threshold = 0,
    use_irlba = FALSE, use_rsvd = TRUE
  )))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("G3S out-of-sample projection reproduces its training scores", {
  set.seed(904)
  training <- matrix(rnorm(36L * 14L), 36L, 14L)
  compression <- compress_features_svd(
    training,
    n_components = 6L,
    variance_threshold = 0,
    use_irlba = FALSE,
    use_rsvd = FALSE
  )

  transformed <- transform_new_data_svd(training, compression)
  scaled <- base::scale(
    training,
    center = compression$center,
    scale = compression$scale
  )
  oracle <- scaled %*% compression$rotation
  oracle_norms <- sqrt(rowSums(oracle^2))
  oracle_norms[oracle_norms == 0] <- 1
  oracle <- oracle / oracle_norms

  expect_equal(transformed, oracle, tolerance = 1e-12)
  expect_equal(transformed, compression$features, tolerance = 1e-12)
})

test_that("ReNA and G3S fail closed or cover disconnected mask components", {
  dims <- c(5L, 2L, 1L)
  mask_values <- array(0, dims)
  mask_values[1:2, , 1] <- 1
  mask_values[4:5, , 1] <- 1
  space <- neuroim2::NeuroSpace(dims)
  mask <- neuroim2::NeuroVol(mask_values, space)
  set.seed(33)
  frames <- lapply(seq_len(6L), function(i) {
    neuroim2::NeuroVol(array(rnorm(prod(dims)), dims), space)
  })
  vec <- do.call(neuroim2::concat, frames)

  expect_error(rena(vec, mask, K = 1L, connectivity = 6L),
               class = "cluster4d_exact_k_infeasible")
  expect_error(
    cluster4d_g3s(
      vec, mask, K = 1L, n_components = 3L, variance_threshold = 0,
      use_irlba = FALSE, use_rsvd = FALSE, max_refinement_iter = 0,
      connectivity = 6L
    ),
    class = "cluster4d_exact_k_infeasible"
  )

  rena_fit <- rena(vec, mask, K = 2L, connectivity = 6L)
  g3s_fit <- cluster4d_g3s(
    vec, mask, K = 2L, n_components = 3L, variance_threshold = 0,
    use_irlba = FALSE, use_rsvd = FALSE, max_refinement_iter = 1,
    connectivity = 6L
  )
  expect_equal(length(unique(rena_fit$cluster)), 2L)
  expect_equal(length(unique(g3s_fit$cluster)), 2L)
  expect_true(topology_labels_connected(rena_fit$cluster, mask, 6L))
  expect_true(topology_labels_connected(g3s_fit$cluster, mask, 6L))
})

test_that("G3S compactness is positive in 2-D and equivariant in physical units", {
  base <- topology_fixture(spacing = c(1, 1, 1), origin = c(0, 0, 0))
  scaled <- topology_fixture(spacing = c(3, 3, 3), origin = c(10, -4, 7))
  args <- list(
    K = 4L, n_components = 4L, variance_threshold = 0,
    use_irlba = FALSE, use_rsvd = FALSE,
    max_refinement_iter = 1L, connectivity = 6L
  )
  fit_base <- do.call(cluster4d_g3s, c(list(base$vec, base$mask), args))
  fit_scaled <- do.call(cluster4d_g3s, c(list(scaled$vec, scaled$mask), args))

  expect_gt(fit_base$parameters$compactness, 0)
  expect_equal(fit_base$metadata$spatial$effective_dimension, 2L)
  expect_equal(
    fit_scaled$parameters$compactness / fit_base$parameters$compactness,
    3,
    tolerance = 1e-12
  )
  expect_identical(fit_base$cluster, fit_scaled$cluster)

  explicit_base <- do.call(
    cluster4d_g3s,
    c(list(base$vec, base$mask, compactness = 2), args)
  )
  explicit_scaled <- do.call(
    cluster4d_g3s,
    c(list(scaled$vec, scaled$mask, compactness = 6), args)
  )
  expect_identical(explicit_base$cluster, explicit_scaled$cluster)

  for (bad in list(0, -1, Inf, NA_real_)) {
    expect_error(
      do.call(cluster4d_g3s, c(list(base$vec, base$mask, compactness = bad), args)),
      "compactness"
    )
  }
})

test_that("ReNA is grid-invariant while reporting physical coordinates", {
  base <- topology_fixture(spacing = c(1, 1, 1), origin = c(0, 0, 0))
  scaled <- topology_fixture(spacing = c(3, 3, 3), origin = c(10, -4, 7))
  fit_base <- rena(base$vec, base$mask, K = 4L, connectivity = 18L)
  fit_scaled <- rena(scaled$vec, scaled$mask, K = 4L, connectivity = 18L)
  expect_identical(fit_base$cluster, fit_scaled$cluster)
  expect_equal(
    as.numeric(stats::dist(fit_scaled$coord_centers)),
    3 * as.numeric(stats::dist(fit_base$coord_centers)),
    tolerance = 1e-10
  )
  expect_false(neurocluster:::cluster4d_method_contract("rena")$spatial_weight)
  expect_error(
    cluster4d(
      base$vec, base$mask, n_clusters = 4L, method = "rena",
      spatial_weight = 0.5, connectivity = 18L
    ),
    "spatial_weight is not supported"
  )
})

test_that("G3S refinement retains the declared spatial term", {
  labels <- c(1L, 1L, 2L, 2L, 2L)
  features <- rbind(rep(1, 5L), rep(0, 5L))
  coords <- cbind(0:4, 0, 0)
  neighbors <- matrix(
    c(0, 2, 1, 3, 2, 4, 3, 5, 4, 0),
    nrow = 5L, byrow = TRUE
  )
  spatial <- refine_boundaries_g3s_cpp(
    labels, features, coords, neighbors,
    alpha = 0, compactness = 1, max_iter = 2L
  )
  feature_only <- refine_boundaries_g3s_cpp(
    labels, features, coords, neighbors,
    alpha = 1, compactness = 1, max_iter = 2L
  )
  expect_identical(as.integer(spatial), labels)
  expect_false(identical(as.integer(spatial), as.integer(feature_only)))
})

test_that("G3S connectivity repair absorbs fragmented labels as units", {
  dims <- c(7L, 1L, 1L)
  mask <- neuroim2::NeuroVol(array(1, dims), neuroim2::NeuroSpace(dims))
  graph <- build_grid_neighbors_g3s(mask, seq_len(prod(dims)), 6L)
  labels <- c(1L, 1L, 2L, 1L, 2L, 2L, 2L)
  features <- rbind(c(0, 0, 0.2, 0.8, 1, 1, 1))

  repaired <- enforce_label_connectivity_cpp(
    labels, features, graph$coords, graph$nn.index,
    alpha = 1, compactness = 1
  )
  normalized <- .exact_k_connected_labels(
    repaired, graph$graph, graph$edges
  )

  expect_identical(length(unique(repaired)), 2L)
  expect_identical(length(unique(normalized)), 2L)
  expect_identical(as.integer(repaired), c(1L, 1L, 1L, 2L, 2L, 2L, 2L))
})

test_that("G3S tied-gradient seeds cover equal distant lobes", {
  n_each <- 1000L
  coords <- cbind(
    c(seq_len(n_each), 3000L + seq_len(n_each)),
    0,
    0
  )
  features <- matrix(1, nrow = 2L * n_each, ncol = 1L)
  seeds <- find_gradient_seeds_g3s(
    features, coords, K = 60L, k_neighbors = 26L,
    distance = "euclidean", spatial_scale = 2L * n_each / 60L
  )
  per_lobe <- tabulate(1L + (seeds > n_each), nbins = 2L)

  expect_identical(length(seeds), 60L)
  expect_true(all(per_lobe >= 20L), info = paste(per_lobe, collapse = "/"))
})

test_that("final ReNA and refined G3S parcels are flood-fill connected", {
  fixture <- topology_fixture()
  rena_fit <- rena(fixture$vec, fixture$mask, K = 5L, connectivity = 6L)
  g3s_fit <- cluster4d_g3s(
    fixture$vec, fixture$mask,
    K = 5L, n_components = 4L, variance_threshold = 0,
    use_irlba = FALSE, use_rsvd = FALSE,
    max_refinement_iter = 3L, connectivity = 6L
  )
  expect_true(topology_labels_connected(rena_fit$cluster, fixture$mask, 6L))
  expect_true(topology_labels_connected(g3s_fit$cluster, fixture$mask, 6L))
  expect_equal(length(unique(rena_fit$cluster)), 5L)
  expect_equal(length(unique(g3s_fit$cluster)), 5L)
})
