library(testthat)
library(neuroim2)

.run_snic_core <- function(features, seeds, spatial_mix = 0.5, s = 100,
                           coords = cbind(seq_len(ncol(features)) - 1L, 0L, 0L)) {
  n <- ncol(features)
  dims <- c(n, 1L, 1L)
  labels <- array(0L, dims)
  mask <- array(1L, dims)
  norm_coords <- matrix(as.numeric(coords), ncol = 3L)
  lookup <- array(seq_len(n) - 1L, dims)

  snic_main_optimized(
    labels, mask, norm_coords[seeds + 1L, , drop = FALSE],
    as.integer(seeds), matrix(as.integer(coords), ncol = 3L), norm_coords,
    as.matrix(features), length(seeds), s, spatial_mix, lookup
  )
}

.snic_reference <- function(features, seeds, spatial_mix, s, coords) {
  features <- as.matrix(features)
  coords <- as.matrix(coords)
  n <- ncol(features)
  k <- length(seeds)
  labels <- integer(n)
  counts <- rep.int(1L, k)
  sum_features <- features[, seeds + 1L, drop = FALSE]
  sum_coords <- t(coords[seeds + 1L, , drop = FALSE])
  labels[seeds + 1L] <- seq_len(k)
  best <- rep(Inf, n)
  best[seeds + 1L] <- 0
  queue <- data.frame(
    index = seeds,
    label = seq_len(k),
    distance = rep(0, k),
    is_seed = rep(TRUE, k)
  )
  assignment_order <- seeds + 1L
  assignment_labels <- seq_len(k)
  seeds_expanded <- 0L

  while (nrow(queue) &&
         (length(assignment_order) < n || seeds_expanded < k)) {
    take <- order(
      queue$distance, -as.integer(queue$is_seed),
      queue$label, queue$index
    )[1L]
    event <- queue[take, , drop = FALSE]
    queue <- queue[-take, , drop = FALSE]

    expand <- FALSE
    if (event$is_seed) {
      seeds_expanded <- seeds_expanded + 1L
      expand <- TRUE
    } else if (labels[event$index + 1L] == 0L) {
      cluster <- event$label
      labels[event$index + 1L] <- cluster
      assignment_order <- c(assignment_order, event$index + 1L)
      assignment_labels <- c(assignment_labels, cluster)
      counts[cluster] <- counts[cluster] + 1L
      sum_features[, cluster] <-
        sum_features[, cluster] + features[, event$index + 1L]
      sum_coords[, cluster] <-
        sum_coords[, cluster] + coords[event$index + 1L, ]
      expand <- TRUE
    }

    if (!expand) next
    cluster <- event$label
    feature_center <- sum_features[, cluster] / counts[cluster]
    if (nrow(features) > 1L) {
      feature_norm <- sqrt(sum(feature_center^2))
      if (feature_norm > 0) feature_center <- feature_center / feature_norm
    }
    coord_center <- sum_coords[, cluster] / counts[cluster]
    delta <- abs(sweep(coords, 2L, coords[event$index + 1L, ], "-"))
    neighbors <- which(rowSums(delta > 1) == 0 & rowSums(delta) > 0) - 1L

    for (neighbor in neighbors) {
      if (labels[neighbor + 1L] != 0L) next
      feature_distance <- sum(
        (feature_center - features[, neighbor + 1L])^2
      )
      spatial_distance <- sum(
        (coord_center - coords[neighbor + 1L, ])^2
      ) / s
      distance <- spatial_mix * spatial_distance +
        (1 - spatial_mix) * feature_distance
      if (distance < best[neighbor + 1L]) {
        best[neighbor + 1L] <- distance
        queue <- rbind(
          queue,
          data.frame(
            index = neighbor, label = cluster, distance = distance,
            is_seed = FALSE
          )
        )
      }
    }
  }

  list(
    labels = labels,
    counts = counts,
    assignment_order = assignment_order,
    assignment_labels = assignment_labels
  )
}

test_that("SNIC seeds own their voxels and are counted exactly once", {
  output <- .run_snic_core(
    matrix(0, nrow = 1L, ncol = 3L), seeds = c(0L, 2L)
  )

  expect_identical(as.integer(output), c(1L, 1L, 2L))
  expect_identical(as.integer(output)[c(1L, 3L)], c(1L, 2L))
  expect_identical(
    attr(output, "snic_centroid_counts"),
    tabulate(as.integer(output), nbins = 2L)
  )
  expect_identical(attr(output, "snic_assignment_order"), c(1L, 3L, 2L))
  expect_identical(attr(output, "snic_assignment_labels"), c(1L, 2L, 1L))
  expect_identical(sort(attr(output, "snic_assignment_order")), 1:3)
})

test_that("optimized SNIC agrees step by step with a slow reference", {
  coords <- cbind(0:6, 0L, 0L)
  features <- rbind(
    c(-1, -0.8, -0.2, 0.1, 0.3, 0.8, 1),
    c(0.2, 0.4, 0.9, 1, 0.8, 0.3, 0.1),
    c(1, 0.8, 0.4, 0, -0.3, -0.7, -1)
  )
  features <- neurocluster:::.snic_normalize_features(features)
  seeds <- c(0L, 3L, 6L)
  expected <- .snic_reference(features, seeds, 0.35, 4, coords)
  observed <- .run_snic_core(features, seeds, 0.35, 4, coords)

  expect_identical(as.integer(observed), expected$labels)
  expect_identical(attr(observed, "snic_centroid_counts"), expected$counts)
  expect_identical(
    attr(observed, "snic_assignment_order"), expected$assignment_order
  )
  expect_identical(
    attr(observed, "snic_assignment_labels"), expected$assignment_labels
  )
})

test_that("SNIC feature normalization is coherent by data modality", {
  multi <- matrix(c(1, 2, 4, 2, 3, 8, 5, 5, 5), nrow = 3L)
  normalized_multi <- neurocluster:::.snic_normalize_features(multi)
  expect_equal(colMeans(normalized_multi), rep(0, 3), tolerance = 1e-14)
  expect_equal(sqrt(colSums(normalized_multi^2)), c(1, 1, 0), tolerance = 1e-14)
  expect_identical(
    attr(normalized_multi, "snic_normalization"),
    "per_voxel_centered_unit_l2"
  )

  structural <- neurocluster:::.snic_normalize_features(matrix(c(0, 0, 10), nrow = 1L))
  expect_equal(mean(structural), 0, tolerance = 1e-14)
  expect_equal(stats::sd(as.numeric(structural)), 1, tolerance = 1e-14)
  expect_identical(
    attr(structural, "snic_normalization"),
    "global_zscore_across_voxels"
  )
})

test_that("single-frame intensity remains an active SNIC feature", {
  low_middle <- neurocluster:::.snic_normalize_features(
    matrix(c(0, 0, 10), nrow = 1L)
  )
  high_middle <- neurocluster:::.snic_normalize_features(
    matrix(c(0, 10, 10), nrow = 1L)
  )

  low_result <- .run_snic_core(low_middle, c(0L, 2L), spatial_mix = 0)
  high_result <- .run_snic_core(high_middle, c(0L, 2L), spatial_mix = 0)
  expect_identical(as.integer(low_result)[2L], 1L)
  expect_identical(as.integer(high_result)[2L], 2L)
})

test_that("public SNIC reports requested and actual K with core diagnostics", {
  dims <- c(5L, 5L, 2L)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  values <- array(seq_len(prod(dims)), dims)
  vec <- NeuroVol(values, NeuroSpace(dims))

  result <- snic(vec, mask, K = 5L, compactness = 5)
  expect_identical(result$parameters$requested_K, 5L)
  expect_identical(result$parameters$actual_K, 5L)
  expect_identical(result$actual_k, 5L)
  expect_identical(result$parameters$spatial_mix, 0.5)
  expect_identical(
    result$metadata$snic$centroid_counts,
    tabulate(result$labels, nbins = result$actual_k)
  )
  expect_identical(
    sort(result$metadata$snic$assignment_order),
    seq_along(result$labels)
  )
})

test_that("SNIC rejects its inactive iteration argument and invalid mixtures", {
  dims <- c(3L, 3L, 2L)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  vec <- NeuroVol(array(seq_len(prod(dims)), dims), NeuroSpace(dims))

  expect_error(snic(vec, mask, K = 2L, max_iter = 1L), "not supported")
  expect_error(snic(vec, mask, K = 2L, compactness = -1), "between 0 and 10")
  expect_error(snic(vec, mask, K = 2L, compactness = 11), "between 0 and 10")
  expect_error(
    cluster4d_snic(vec, mask, n_clusters = 2L, max_iterations = 1L),
    "not supported"
  )
  expect_error(
    neurocluster:::map_cluster4d_params("snic", max_iter = 1L),
    "not supported"
  )
})
