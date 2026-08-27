slice_msf_oracle_detrend_zscore <- function(y) {
  stopifnot(is.numeric(y), length(y) >= 2L, all(is.finite(y)))
  n_time <- length(y)
  time <- 0:(n_time - 1L)
  mean_time <- 0.5 * (n_time - 1L)
  var_time <- (n_time * n_time - 1) / 12
  mean_y <- mean(y)
  beta <- (mean(time * y) - mean_time * mean_y) / var_time
  alpha <- mean_y - beta * mean_time
  residual <- y - (alpha + beta * time)
  centered <- residual - mean(residual)
  sum_squares <- sum(centered * centered)
  if (sum_squares <= 0) sum_squares <- 1
  standard_deviation <- sqrt(sum_squares / max(1L, n_time - 1L))
  if (standard_deviation <= 0) standard_deviation <- 1
  centered / standard_deviation
}

slice_msf_oracle_adjacent_correlation <- function(z) {
  n_pairs <- length(z) %/% 2L
  if (n_pairs < 2L) return(0)
  left <- z[2L * seq_len(n_pairs) - 1L]
  right <- z[2L * seq_len(n_pairs)]
  if (sum((left - mean(left))^2) <= 0 ||
      sum((right - mean(right))^2) <= 0) {
    return(0)
  }
  value <- stats::cor(left, right)
  if (!is.finite(value)) 0 else max(-1, min(1, value))
}

slice_msf_oracle_sketch <- function(y, frequencies, weights = NULL) {
  z <- slice_msf_oracle_detrend_zscore(y)
  n_time <- length(z)
  if (is.null(weights)) weights <- rep(1, length(frequencies))
  basis <- vapply(frequencies, function(frequency) {
    sqrt(2 / n_time) * cos(
      pi * ((0:(n_time - 1L)) + 0.5) * frequency / n_time
    )
  }, numeric(n_time))
  coefficients <- as.vector(crossprod(basis, z)) * weights
  norm <- sqrt(sum(coefficients * coefficients))
  coefficients / max(1e-6, norm)
}

slice_msf_oracle_edges <- function(mask, sketch, nbhd = 8L,
                                   stitch_z = TRUE) {
  dims <- dim(mask)
  included <- as.vector(as.array(mask)) > 0
  gids <- which(included)
  global_to_local <- integer(prod(dims))
  global_to_local[gids] <- seq_along(gids)
  offsets <- rbind(c(1L, 0L, 0L), c(0L, 1L, 0L), c(0L, 0L, 1L))
  if (nbhd == 8L) {
    extra <- list()
    for (dz in 0:1) for (dy in -1:1) for (dx in -1:1) {
      if (dx == 0L && dy == 0L && dz == 0L) next
      if (dz == 0L && dy == 0L && dx <= 0L) next
      if (dz == 0L && abs(dx) + abs(dy) <= 1L) next
      if (dz == 1L && dx == 0L && dy == 0L) next
      extra[[length(extra) + 1L]] <- c(dx, dy, dz)
    }
    offsets <- rbind(offsets, do.call(rbind, extra))
  }
  result <- vector("list", 0L)
  for (global in gids) {
    coordinate <- arrayInd(global, dims)
    for (row in seq_len(nrow(offsets))) {
      offset <- offsets[row, ]
      if (!stitch_z && offset[[3L]] != 0L) next
      neighbor_coordinate <- coordinate + offset
      if (any(neighbor_coordinate < 1L) ||
          any(neighbor_coordinate > dims)) next
      neighbor_global <- neighbor_coordinate[[1L]] +
        dims[[1L]] * (neighbor_coordinate[[2L]] - 1L) +
        dims[[1L]] * dims[[2L]] * (neighbor_coordinate[[3L]] - 1L)
      neighbor_local <- global_to_local[[neighbor_global]]
      if (neighbor_local == 0L) next
      local <- global_to_local[[global]]
      similarity <- sum(sketch[, global] * sketch[, neighbor_global])
      distance <- 1 - max(-1, min(1, similarity))
      result[[length(result) + 1L]] <- c(local, neighbor_local, distance)
    }
  }
  if (!length(result)) {
    return(data.frame(left = integer(), right = integer(), distance = numeric()))
  }
  values <- do.call(rbind, result)
  data.frame(
    left = as.integer(values[, 1L]),
    right = as.integer(values[, 2L]),
    distance = as.numeric(values[, 3L])
  )
}

slice_msf_oracle_fh <- function(n_vertices, edges, fh_scale, min_size) {
  parent <- seq_len(n_vertices)
  rank <- integer(n_vertices)
  size <- rep(1L, n_vertices)
  internal <- numeric(n_vertices)
  find_root <- function(node) {
    while (node != parent[[node]]) {
      parent[[node]] <<- parent[[parent[[node]]]]
      node <- parent[[node]]
    }
    node
  }
  unite <- function(left, right, weight) {
    left <- find_root(left)
    right <- find_root(right)
    if (left == right) return(left)
    if (rank[[left]] < rank[[right]]) {
      swap <- left
      left <- right
      right <- swap
    }
    parent[[right]] <<- left
    size[[left]] <<- size[[left]] + size[[right]]
    internal[[left]] <<- max(internal[[left]], internal[[right]], weight)
    if (rank[[left]] == rank[[right]]) rank[[left]] <<- rank[[left]] + 1L
    left
  }

  quantized <- as.integer(pmax(0, pmin(2, edges$distance)) * 32767.5)
  edges <- edges[order(quantized, seq_len(nrow(edges))), , drop = FALSE]
  for (row in seq_len(nrow(edges))) {
    left <- find_root(edges$left[[row]])
    right <- find_root(edges$right[[row]])
    if (left == right) next
    left_threshold <- internal[[left]] + fh_scale / size[[left]]
    right_threshold <- internal[[right]] + fh_scale / size[[right]]
    if (edges$distance[[row]] <= left_threshold &&
        edges$distance[[row]] <= right_threshold) {
      unite(left, right, edges$distance[[row]])
    }
  }

  for (iteration in seq_len(32L)) {
    best_weight <- rep(Inf, n_vertices)
    best_neighbor <- integer(n_vertices)
    for (row in seq_len(nrow(edges))) {
      left <- find_root(edges$left[[row]])
      right <- find_root(edges$right[[row]])
      if (left == right) next
      weight <- edges$distance[[row]]
      if (size[[left]] < min_size && weight < best_weight[[left]]) {
        best_weight[[left]] <- weight
        best_neighbor[[left]] <- right
      }
      if (size[[right]] < min_size && weight < best_weight[[right]]) {
        best_weight[[right]] <- weight
        best_neighbor[[right]] <- left
      }
    }
    changed <- FALSE
    for (node in seq_len(n_vertices)) {
      root <- find_root(node)
      if (root != node) next
      if (size[[root]] < min_size && best_neighbor[[root]] > 0L) {
        unite(root, best_neighbor[[root]], best_weight[[root]])
        changed <- TRUE
      }
    }
    if (!changed) break
  }
  roots <- vapply(seq_len(n_vertices), find_root, integer(1L))
  as.integer(match(roots, unique(roots)))
}

slice_msf_truth_is_connected <- function(fixture) {
  graph <- neurocluster:::.exact_k_graph(
    fixture$mask, fixture$contract$connectivity
  )
  split <- neurocluster:::.exact_k_connected_labels(
    fixture$truth, graph$graph, graph$edges
  )
  length(unique(split)) == length(unique(fixture$truth))
}

slice_msf_all_connected <- function(cluster_map, connectivity = 6L) {
  offsets <- as.matrix(expand.grid(-1L:1L, -1L:1L, -1L:1L))
  offsets <- offsets[rowSums(abs(offsets)) > 0L, , drop = FALSE]
  if (connectivity == 6L) {
    offsets <- offsets[rowSums(abs(offsets)) == 1L, , drop = FALSE]
  } else if (connectivity == 18L) {
    offsets <- offsets[rowSums(abs(offsets)) <= 2L, , drop = FALSE]
  }
  dims <- dim(cluster_map)
  all(vapply(sort(unique(cluster_map[cluster_map > 0L])), function(label) {
    members <- which(cluster_map == label, arr.ind = TRUE)
    seen <- array(FALSE, dims)
    queue <- matrix(members[1L, ], nrow = 1L)
    visited <- 0L
    while (nrow(queue) > 0L) {
      current <- queue[1L, ]
      queue <- queue[-1L, , drop = FALSE]
      if (seen[current[1L], current[2L], current[3L]]) next
      if (cluster_map[current[1L], current[2L], current[3L]] != label) next
      seen[current[1L], current[2L], current[3L]] <- TRUE
      visited <- visited + 1L
      candidates <- sweep(offsets, 2L, current, `+`)
      valid <- rowSums(candidates < 1L) == 0L &
        candidates[, 1L] <= dims[[1L]] &
        candidates[, 2L] <= dims[[2L]] &
        candidates[, 3L] <= dims[[3L]]
      if (any(valid)) queue <- rbind(queue, candidates[valid, , drop = FALSE])
    }
    visited == nrow(members)
  }, logical(1L)))
}

slice_msf_quality_diagnostic <- function(fit, truth, threshold, fixture,
                                         estimand, seed) {
  observed <- clustering_accuracy(fit$cluster, truth)$ari
  sizes <- sort(tabulate(fit$cluster))
  repair <- fit$metadata$exact_k_repair
  natural_k <- if (is.null(repair)) {
    length(unique(fit$cluster))
  } else {
    repair$normalized_k
  }
  requested_k <- if (is.null(repair)) {
    NA_integer_
  } else {
    repair$requested_k
  }
  min_size <- fit$parameters$min_size
  if (is.null(min_size)) min_size <- NA_integer_
  sprintf(
    paste0(
      "estimand=%s; fixture=%s; seed=%d; observed_ari=%.6f; ",
      "threshold=%.2f; natural_k=%s; repaired_k=%s; min_size=%s; ",
      "smallest=%d; largest=%d"
    ),
    estimand, fixture, as.integer(seed), observed, threshold,
    ifelse(is.na(natural_k), "NA", as.character(natural_k)),
    ifelse(is.na(requested_k), "NA", as.character(requested_k)),
    ifelse(is.na(min_size), "NA", as.character(min_size)),
    min(sizes), max(sizes)
  )
}

slice_msf_make_spherical_fixture <- function() {
  dims <- c(15L, 15L, 15L)
  n_clusters <- 8L
  n_time <- 100L
  set.seed(42)
  centers <- matrix(0, n_clusters, 3L)
  centers[1L, ] <- c(runif(1L, 3, 12), runif(1L, 3, 12), runif(1L, 3, 12))
  minimum_distance <- min(dims) / (n_clusters^(1 / 3) + 1)
  for (cluster in 2:n_clusters) {
    best <- NULL
    best_distance <- 0
    for (attempt in seq_len(100L)) {
      candidate <- c(runif(1L, 3, 12), runif(1L, 3, 12), runif(1L, 3, 12))
      distances <- sqrt(rowSums(
        (centers[seq_len(cluster - 1L), , drop = FALSE] -
           matrix(candidate, cluster - 1L, 3L, byrow = TRUE))^2
      ))
      if (min(distances) > best_distance) {
        best_distance <- min(distances)
        best <- candidate
      }
      if (min(distances) > minimum_distance) break
    }
    centers[cluster, ] <- best
  }
  base_radius <- min(dims) / (2 * n_clusters^(1 / 3))
  radii <- base_radius * (1 + 0.5 * (runif(n_clusters) - 0.5))
  coordinates <- expand.grid(x = seq_len(dims[[1L]]), y = seq_len(dims[[2L]]),
                             z = seq_len(dims[[3L]]))
  truth_full <- integer(nrow(coordinates))
  for (voxel in seq_len(nrow(coordinates))) {
    coordinate <- as.numeric(coordinates[voxel, ])
    distances <- sqrt(rowSums((centers - matrix(
      coordinate, n_clusters, 3L, byrow = TRUE
    ))^2))
    nearest <- which.min(distances)
    if (distances[[nearest]] <= radii[[nearest]] * 1.5) truth_full[[voxel]] <- nearest
  }
  included <- truth_full > 0L
  signals <- qr.Q(qr(matrix(rnorm(n_time * n_clusters), n_time, n_clusters)))
  values <- array(0, c(dims, n_time))
  for (voxel in which(included)) {
    coordinate <- as.integer(coordinates[voxel, ])
    values[coordinate[[1L]], coordinate[[2L]], coordinate[[3L]], ] <-
      signals[, truth_full[[voxel]]] + rnorm(n_time, sd = 0.3)
  }
  mask_values <- array(as.integer(included), dims)
  list(
    vec = neuroim2::NeuroVec(values, neuroim2::NeuroSpace(c(dims, n_time))),
    mask = neuroim2::NeuroVol(mask_values, neuroim2::NeuroSpace(dims)),
    truth = truth_full[included],
    contract = list(
      estimand = "eight connected spherical/Voronoi parcels",
      connectivity = 26L,
      seed = 42L,
      invariances = c("offset", "linear_trend", "positive_scale")
    )
  )
}

slice_msf_make_block_fixture <- function(dims = c(12L, 12L, 3L),
                                         n_time = 48L, seed = 701L) {
  set.seed(seed)
  time <- 0:(n_time - 1L)
  basis <- vapply(seq_len(16L), function(frequency) {
    cos(pi * (time + 0.5) * frequency / n_time)
  }, numeric(n_time))
  truth_array <- array(0L, dims)
  values <- array(0, c(dims, n_time))
  for (z in seq_len(dims[[3L]])) for (y in seq_len(dims[[2L]])) {
    for (x in seq_len(dims[[1L]])) {
      region <- 1L + (x - 1L) %/% 3L + 4L * ((y - 1L) %/% 6L)
      coefficients <- rnorm(16L, sd = 0.04)
      coefficients[[region]] <- 1
      coefficients[[region + 8L]] <- 0.4
      values[x, y, z, ] <- as.vector(basis %*% coefficients) +
        rnorm(n_time, sd = 0.02)
      truth_array[x, y, z] <- region
    }
  }
  mask <- neuroim2::NeuroVol(array(1, dims), neuroim2::NeuroSpace(dims))
  list(
    vec = neuroim2::NeuroVec(values, neuroim2::NeuroSpace(c(dims, n_time))),
    mask = mask,
    truth = as.integer(truth_array),
    contract = list(
      estimand = "eight connected axis-aligned DCT parcels",
      connectivity = 26L,
      seed = as.integer(seed),
      invariances = c("offset", "linear_trend", "positive_scale")
    )
  )
}

slice_msf_make_cube_fixture <- function() {
  fixture <- slice_msf_make_block_fixture(c(9L, 9L, 3L), 48L, 811L)
  fixture$contract$estimand <- "six connected cuboid parcels"
  fixture
}

slice_msf_make_gradient_fixture <- function() {
  fixture <- slice_msf_make_block_fixture(c(12L, 12L, 4L), 64L, 913L)
  values <- as.array(fixture$vec)
  for (x in seq_len(dim(values)[[1L]])) {
    alpha <- (x - 1L) / max(1L, dim(values)[[1L]] - 1L)
    values[x, , , ] <- (1 - 0.25 * alpha) * values[x, , , ] +
      0.25 * alpha * values[pmin(x + 1L, dim(values)[[1L]]), , , ]
  }
  fixture$vec <- neuroim2::NeuroVec(
    values, neuroim2::NeuroSpace(dim(values))
  )
  fixture$contract$estimand <- "eight connected parcels with graded x-boundaries"
  fixture
}
