test_that("empty supervoxel labels are distinctly and deterministically reseeded", {
  feature_mat <- matrix(c(10, 10, 10, 10, 4, 4), nrow = 1L)
  coords <- cbind(c(0, 10, 10, 10, 20, 21), 0, 0)
  labels <- c(1L, 1L, 1L, 1L, 2L, 2L)

  state1 <- neurocluster:::.supervoxel_center_state(
    labels, feature_mat, coords, K = 4L,
    alpha = 0, sigma1 = 1, sigma2 = 5,
    parallel = FALSE
  )
  state2 <- neurocluster:::.supervoxel_center_state(
    labels, feature_mat, coords, K = 4L,
    alpha = 0, sigma1 = 1, sigma2 = 5,
    parallel = TRUE
  )

  expect_identical(state1$empty_labels, c(3L, 4L))
  expect_length(unique(state1$reseeded_voxels), 2L)
  expect_true(all(state1$counts > 0L))
  expect_identical(sort(unique(state1$labels)), 1:4)
  expect_identical(state1$labels, state2$labels)
  expect_equal(state1$centers, state2$centers, tolerance = 1e-12)
  expect_equal(state1$coord_centers, state2$coord_centers, tolerance = 1e-12)
})


test_that("native assignment rejects an unoccupied phantom center", {
  data <- matrix(rep(10, 4), nrow = 1L)
  coords <- rbind(c(0, 10, 10, 10), c(0, 0, 0, 0), c(0, 0, 0, 0))
  centers <- compute_centroids_parallel_fast(
    c(0L, 0L, 0L, 0L), data, coords, 2L
  )
  nn_index <- matrix(c(1L, 0L, 0L, 0L), ncol = 1L)
  nn_dist <- matrix(rep(1, 4), ncol = 1L)

  expect_error(
    fused_assignment_binned(
      nn_index, nn_dist, rep(0L, 4L), coords,
      centers$centers, centers$coord_centers, centers$counts, data,
      20, 1, 5, 0, window_factor = 100, bin_expand = 100L
    ),
    "every center is non-empty"
  )
  expect_error(
    fused_assignment_parallel_binned(
      nn_index, nn_dist, rep(0L, 4L), coords,
      centers$centers, centers$coord_centers, centers$counts, data,
      20, 1, 5, 0, grain_size = 1L,
      window_factor = 100, bin_expand = 100L
    ),
    "every center is non-empty"
  )
})


test_that("instrumented iterations keep K labels and K finite centers", {
  n <- 24L
  # Coincident coordinates plus identical features make every score tie. The
  # native deterministic tie rule collapses assignments to the first label,
  # forcing the exact-K reseed path on every non-converged iteration.
  coords <- matrix(0, nrow = n, ncol = 3L)
  feature_mat <- matrix(0, nrow = 2L, ncol = n)
  initclus <- rep(1:4, each = n / 4L)

  fit <- suppressWarnings(neurocluster:::supervoxel_cluster_fit(
    feature_mat, coords, K = 4L, sigma1 = 1, sigma2 = 2,
    alpha = 1, iterations = 4L, connectivity = 2L,
    initclus = initclus, use_gradient = FALSE, parallel = FALSE,
    converge_thresh = 1e-12, trace = TRUE
  ))

  expect_identical(sort(unique(fit$clusters)), 1:4)
  expect_identical(length(fit$counts), 4L)
  expect_true(all(fit$counts > 0L))
  expect_true(length(fit$reseed_events) >= 1L)
  expect_true(length(fit$iteration_trace) >= 1L)
  for (entry in fit$iteration_trace) {
    expect_length(entry$counts, 4L)
    expect_true(all(entry$counts > 0L))
    expect_identical(entry$n_feature_centers, 4L)
    expect_identical(entry$n_coord_centers, 4L)
    expect_length(unique(entry$reseeded_voxels), length(entry$reseeded_voxels))
  }
})


test_that("parallel and sequential fits obey an exact deterministic contract", {
  dims <- c(16L, 16L, 5L)
  coords <- arrayInd(seq_len(prod(dims)), dims)
  n <- nrow(coords)
  feature_mat <- rbind(
    sin(seq_len(n) / 11),
    cos(seq_len(n) / 17),
    (coords[, 1] - mean(coords[, 1])) / max(coords[, 1])
  )
  initclus <- as.integer(cut(
    coords[, 1] + 16 * (coords[, 2] - 1L),
    breaks = 8L, labels = FALSE
  ))

  args <- list(
    feature_mat = feature_mat, coords = coords, K = 8L,
    sigma1 = 1.5, sigma2 = 3, alpha = 0.55,
    iterations = 3L, connectivity = 6L, initclus = initclus,
    use_gradient = FALSE, grain_size = 32L, converge_thresh = 1e-12,
    trace = TRUE
  )
  sequential <- do.call(
    neurocluster:::supervoxel_cluster_fit,
    c(args, list(parallel = FALSE))
  )
  parallel <- do.call(
    neurocluster:::supervoxel_cluster_fit,
    c(args, list(parallel = TRUE))
  )

  expect_false(sequential$parallel_used)
  expect_true(parallel$parallel_used)
  expect_identical(parallel$clusters, sequential$clusters)
  expect_equal(parallel$centers, sequential$centers, tolerance = 1e-12)
  expect_equal(
    parallel$coord_centers, sequential$coord_centers, tolerance = 1e-12
  )
  expect_identical(parallel$counts, sequential$counts)
})
