slice_msf_public_fit <- function(fixture, seed = 1L, num_runs = 1L,
                                 consensus = FALSE) {
  slice_msf(
    fixture$vec, fixture$mask,
    target_k_global = length(unique(fixture$truth)),
    r = 16L, compactness = 4, min_size = 2L,
    num_runs = num_runs, consensus = consensus,
    stitch_z = TRUE, nbhd = 8L, seed = seed
  )
}

slice_msf_vec_from_array <- function(values, spatial_spacing = NULL) {
  space <- if (is.null(spatial_spacing)) {
    neuroim2::NeuroSpace(dim(values))
  } else {
    neuroim2::NeuroSpace(dim(values), spacing = as.numeric(spatial_spacing))
  }
  neuroim2::NeuroVec(values, space)
}

test_that("predeclared public Slice-MSF quality gates pass", {
  slice_fixture <- slice_msf_make_block_fixture(
    dims = c(12L, 12L, 1L), n_time = 48L, seed = 701L
  )
  direct <- slice_msf_public_fit(slice_fixture)
  wrapped <- cluster4d(
    slice_fixture$vec, slice_fixture$mask,
    n_clusters = 8L, method = "slice_msf",
    spatial_weight = 0.4, min_size = 2L, r = 16L,
    num_runs = 1L, consensus = FALSE, verbose = FALSE
  )
  for (fit_name in c("direct", "wrapper")) {
    fit <- if (fit_name == "direct") direct else wrapped
    observed <- clustering_accuracy(fit$cluster, slice_fixture$truth)$ari
    expect_true(
      observed >= 0.85,
      info = slice_msf_quality_diagnostic(
        fit, slice_fixture$truth, 0.85,
        paste0("slice_friendly_", fit_name),
        slice_fixture$contract$estimand, 701L
      )
    )
  }

  cube_fixture <- make_synthetic_clusters(
    n_time = 100L, noise_sd = 0, n_clusters = 27L, seed = 42L
  )
  cube_fit <- cluster4d(
    cube_fixture$vec, cube_fixture$mask,
    n_clusters = 27L, method = "slice_msf", verbose = FALSE
  )
  cube_ari <- clustering_accuracy(
    cube_fit$cluster, cube_fixture$true_labels
  )$ari
  expect_true(
    cube_ari >= 0.60,
    info = slice_msf_quality_diagnostic(
      cube_fit, cube_fixture$true_labels, 0.60,
      "cube_fixture_1", "27 separable 3x3x3 cubes", 42L
    )
  )

  ensemble_fixture <- slice_msf_make_block_fixture()
  ensemble <- slice_msf_public_fit(
    ensemble_fixture, seed = 1L, num_runs = 3L, consensus = TRUE
  )
  observed <- clustering_accuracy(
    ensemble$cluster, ensemble_fixture$truth
  )$ari
  expect_true(
    observed >= 0.90,
    info = slice_msf_quality_diagnostic(
      ensemble, ensemble_fixture$truth, 0.90,
      "dct_ensemble_consensus", ensemble_fixture$contract$estimand,
      ensemble_fixture$contract$seed
    )
  )
  run_signatures <- vapply(
    ensemble$runs, function(run) paste(run$labels, collapse = ","),
    character(1L)
  )
  expect_identical(length(unique(run_signatures)), 3L)
})

test_that("declared temporal invariances and spatial permutation are metamorphic", {
  fixture <- slice_msf_make_block_fixture(
    dims = c(12L, 12L, 1L), n_time = 48L, seed = 701L
  )
  base <- slice_msf_public_fit(fixture)
  values <- as.array(fixture$vec)
  dims <- dim(values)[1:3]
  n_time <- dim(values)[4L]
  coordinates <- arrayInd(seq_len(prod(dims)), dims)

  offset <- array(2 + coordinates[, 1L] / dims[[1L]], dims)
  offset_values <- values + array(
    rep(as.vector(offset), n_time), dim(values)
  )
  slope <- array((coordinates[, 2L] - 1L) / max(1L, dims[[2L]] - 1L), dims)
  trend_values <- values + array(
    as.vector(outer(as.vector(slope), 0:(n_time - 1L))), dim(values)
  )
  scale <- array(0.5 + coordinates[, 1L] / dims[[1L]], dims)
  scale_values <- values * array(
    rep(as.vector(scale), n_time), dim(values)
  )

  transformed <- list(
    offset = offset_values,
    linear_trend = trend_values,
    positive_scale = scale_values
  )
  for (name in names(transformed)) {
    changed_fixture <- fixture
    changed_fixture$vec <- slice_msf_vec_from_array(transformed[[name]])
    changed <- slice_msf_public_fit(changed_fixture)
    expect_identical(
      changed$cluster, base$cluster,
      info = paste(name, "is declared label-identical after detrend/z-score")
    )
  }

  reflected_values <- values[rev(seq_len(dims[[1L]])), , , , drop = FALSE]
  reflected_mask <- neuroim2::NeuroVol(
    as.array(fixture$mask)[rev(seq_len(dims[[1L]])), , , drop = FALSE],
    neuroim2::NeuroSpace(dims)
  )
  reflected_fixture <- fixture
  reflected_fixture$vec <- slice_msf_vec_from_array(reflected_values)
  reflected_fixture$mask <- reflected_mask
  reflected <- slice_msf_public_fit(reflected_fixture)
  inverse_map <- reflected$cluster_map[
    rev(seq_len(dims[[1L]])), , , drop = FALSE
  ]
  expect_equal(
    clustering_accuracy(as.integer(inverse_map), base$cluster)$ari,
    1, tolerance = 1e-12,
    info = "x-reflection is ARI-equivalent after the inverse voxel permutation"
  )

  relabel_map <- rev(seq_len(length(unique(base$cluster))))
  relabeled <- relabel_map[base$cluster]
  expect_equal(
    clustering_accuracy(relabeled, base$cluster)$ari,
    1, tolerance = 1e-12,
    info = "stable label renaming is partition-equivalent"
  )
})

test_that("adversarial temporal series remain finite and non-finite input fails", {
  dims <- c(6L, 2L, 1L)
  n_time <- 24L
  time <- seq(0, 2 * pi, length.out = n_time)
  series <- rbind(
    constant = rep(3, n_time),
    near_constant = 3 + 1e-10 * sin(time),
    zero_adjacent = rep(c(-1, 0, 1, 0), length.out = n_time),
    negative_adjacent = rep(c(-1, 1), length.out = n_time),
    alternating_phase = rep(c(1, -1, -1, 1), length.out = n_time),
    low_frequency = sin(time)
  )
  values <- array(0, c(dims, n_time))
  for (voxel in seq_len(prod(dims))) {
    values[arrayInd(voxel, dims)[1L], arrayInd(voxel, dims)[2L], 1L, ] <-
      series[1L + (voxel - 1L) %% nrow(series), ]
  }
  vec <- slice_msf_vec_from_array(values)
  mask <- neuroim2::NeuroVol(array(1, dims), neuroim2::NeuroSpace(dims))
  native <- slice_msf_single(
    vec, mask, r = 8L, k = 0.25, min_size = 1L,
    nbhd = 8L, stitch_z = TRUE, gamma = 0, z_mult = 0
  )
  fit <- slice_msf(
    vec, mask, target_k_global = 4L,
    r = 8L, compactness = 4, min_size = 1L,
    num_runs = 1L, consensus = FALSE
  )
  expect_true(all(is.finite(native$sketch)))
  expect_true(all(is.finite(native$temporal_smoothness)))
  expect_true(all(fit$cluster > 0L))
  expect_identical(length(unique(fit$cluster)), 4L)
  expect_lte(native$temporal_smoothness[[3L]], 0)
  expect_lt(native$temporal_smoothness[[5L]], 0)

  bad_values <- values
  bad_values[1L, 1L, 1L, 3L] <- NA_real_
  expect_error(
    slice_msf(
      slice_msf_vec_from_array(bad_values), mask,
      r = 8L, min_size = 1L, num_runs = 1L, consensus = FALSE
    ),
    "finite|non-finite"
  )
})

test_that("holes, a thin bridge, disconnected input, and anisotropy are explicit", {
  dims <- c(12L, 9L, 3L)
  n_time <- 32L
  mask_values <- array(0L, dims)
  mask_values[1:5, 1:7, ] <- 1L
  mask_values[3:4, 3:5, 2L] <- 0L
  mask_values[6:7, 4L, 2L] <- 1L
  mask_values[8:10, 2:7, ] <- 1L
  mask_values[12L, 8:9, 1:2] <- 1L
  space <- neuroim2::NeuroSpace(dims, spacing = c(1, 2, 4))
  mask <- neuroim2::NeuroVol(mask_values, space)
  time <- 0:(n_time - 1L)
  values <- array(0, c(dims, n_time))
  for (voxel in which(mask_values > 0L)) {
    coordinate <- arrayInd(voxel, dims)
    frequency <- 1L + (coordinate[[1L]] - 1L) %% 6L
    values[coordinate[[1L]], coordinate[[2L]], coordinate[[3L]], ] <-
      cos(pi * (time + 0.5) * frequency / n_time)
  }
  vec <- slice_msf_vec_from_array(values, c(1, 2, 4))
  natural <- slice_msf(
    vec, mask, r = 8L, compactness = 4, min_size = 2L,
    num_runs = 1L, consensus = FALSE, nbhd = 4L
  )
  repaired <- slice_msf(
    vec, mask, target_k_global = 6L,
    r = 8L, compactness = 4, min_size = 2L,
    num_runs = 1L, consensus = FALSE, nbhd = 4L
  )

  expect_true(slice_msf_all_connected(natural$cluster_map, 6L))
  expect_true(slice_msf_all_connected(repaired$cluster_map, 6L))
  expect_identical(length(unique(repaired$cluster)), 6L)
  expect_true(all(tabulate(repaired$cluster) >= 2L))
  expect_true(all(repaired$cluster_map[mask_values == 0L] == 0L))
  expect_equal(neuroim2::spacing(mask), c(1, 2, 4))

  condition <- tryCatch(
    slice_msf(
      vec, mask, target_k_global = 1L,
      r = 8L, compactness = 4, min_size = 2L,
      num_runs = 1L, consensus = FALSE, nbhd = 4L
    ),
    cluster4d_exact_k_infeasible = identity
  )
  expect_s3_class(condition, "cluster4d_exact_k_infeasible")
  expect_identical(condition$reason, "disconnected_mask_components")
})

test_that("seeded ensembles replay and preserve the caller RNG state", {
  fixture <- slice_msf_make_block_fixture(
    dims = c(8L, 8L, 2L), n_time = 32L, seed = 733L
  )
  set.seed(99173)
  before <- .Random.seed
  first <- slice_msf_public_fit(
    fixture, seed = 17L, num_runs = 3L, consensus = TRUE
  )
  after <- .Random.seed
  second <- slice_msf_public_fit(
    fixture, seed = 17L, num_runs = 3L, consensus = TRUE
  )
  expect_identical(after, before)
  expect_identical(first$cluster, second$cluster)
  expect_identical(first$metadata$ensemble$runs,
                   second$metadata$ensemble$runs)

  different <- slice_msf_public_fit(
    fixture, seed = 18L, num_runs = 3L, consensus = TRUE
  )
  expect_false(identical(
    first$metadata$ensemble$runs,
    different$metadata$ensemble$runs
  ))
})
