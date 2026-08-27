test_that("Slice-MSF fixtures declare connected estimands and invariances", {
  fixtures <- list(
    spherical = slice_msf_make_spherical_fixture(),
    gradient = slice_msf_make_gradient_fixture(),
    cube = slice_msf_make_cube_fixture(),
    dct_ensemble = slice_msf_make_block_fixture()
  )
  for (name in names(fixtures)) {
    fixture <- fixtures[[name]]
    expect_named(
      fixture$contract,
      c("estimand", "connectivity", "seed", "invariances"),
      info = name
    )
    expect_true(nzchar(fixture$contract$estimand), info = name)
    expect_true(slice_msf_truth_is_connected(fixture), info = name)
    expect_setequal(
      fixture$contract$invariances,
      c("offset", "linear_trend", "positive_scale")
    )
  }
})

test_that("independent R preprocessing and FH oracle match native gamma-zero path", {
  dims <- c(3L, 2L, 1L)
  n_time <- 12L
  time <- 0:(n_time - 1L)
  signals <- rbind(
    cos(pi * (time + 0.5) / n_time),
    cos(3 * pi * (time + 0.5) / n_time),
    -cos(pi * (time + 0.5) / n_time)
  )
  values <- array(0, c(dims, n_time))
  for (voxel in seq_len(prod(dims))) {
    coordinate <- arrayInd(voxel, dims)
    signal <- signals[1L + (voxel - 1L) %% 3L, ]
    values[coordinate[[1L]], coordinate[[2L]], coordinate[[3L]], ] <-
      3 + 0.2 * time + signal
  }
  mask <- neuroim2::NeuroVol(array(1, dims), neuroim2::NeuroSpace(dims))
  vec <- neuroim2::NeuroVec(values, neuroim2::NeuroSpace(c(dims, n_time)))
  frequencies <- c(1L, 2L, 3L, 5L)
  mode_weights <- c(1, 0.5, 1.5, 0.75)
  native <- slice_msf_single(
    vec, mask, r = length(frequencies), k = 0.25, min_size = 1L,
    nbhd = 8L, stitch_z = TRUE, gamma = 0, z_mult = 0,
    dct_frequencies = frequencies, dct_weights = mode_weights
  )
  expected_sketch <- vapply(seq_len(prod(dims)), function(voxel) {
    coordinate <- arrayInd(voxel, dims)
    slice_msf_oracle_sketch(
      values[coordinate[[1L]], coordinate[[2L]], coordinate[[3L]], ],
      frequencies, mode_weights
    )
  }, numeric(length(frequencies)))
  expect_equal(native$sketch, expected_sketch, tolerance = 1e-11)
  expected_smoothness <- vapply(seq_len(prod(dims)), function(voxel) {
    coordinate <- arrayInd(voxel, dims)
    slice_msf_oracle_adjacent_correlation(
      slice_msf_oracle_detrend_zscore(
        values[coordinate[[1L]], coordinate[[2L]], coordinate[[3L]], ]
      )
    )
  }, numeric(1L))
  expect_equal(native$temporal_smoothness, expected_smoothness, tolerance = 1e-12)
  expect_equal(native$params$feature_distance, "1-cosine")

  edges <- slice_msf_oracle_edges(mask, expected_sketch, 8L, TRUE)
  expected_labels <- slice_msf_oracle_fh(
    prod(dims), edges, fh_scale = 0.25, min_size = 1L
  )
  expect_equal(
    clustering_accuracy(native$labels, expected_labels)$ari,
    1,
    tolerance = 1e-12
  )
})

test_that("supported Slice-MSF defaults recover the spherical estimand", {
  expect_identical(formals(slice_msf)$gamma, 0)
  expect_identical(formals(slice_msf)$z_mult, 0.0)
  expect_identical(formals(cluster4d_slice_msf)$gamma, 0)
  expect_identical(formals(cluster4d_slice_msf)$z_mult, 0)
  fixture <- slice_msf_make_spherical_fixture()
  direct <- slice_msf(
    fixture$vec, fixture$mask,
    target_k_global = -1L, r = 16L, compactness = 4,
    min_size = 2L, num_runs = 1L, consensus = FALSE,
    stitch_z = TRUE, nbhd = 8L
  )
  wrapped <- cluster4d(
    fixture$vec, fixture$mask, n_clusters = 8L, method = "slice_msf",
    spatial_weight = 0.4, min_size = 2L, r = 16L,
    num_runs = 1L, consensus = FALSE, verbose = FALSE
  )
  ablation <- cluster4d(
    fixture$vec, fixture$mask, n_clusters = 8L, method = "slice_msf",
    spatial_weight = 0.4, min_size = 2L, r = 16L,
    num_runs = 1L, consensus = FALSE, gamma = 0, z_mult = 0,
    verbose = FALSE
  )

  expect_gte(length(unique(direct$cluster)), 2L)
  expect_gte(clustering_accuracy(wrapped$cluster, fixture$truth)$ari, 0.80)
  expect_gte(clustering_accuracy(ablation$cluster, fixture$truth)$ari, 0.80)
  expect_gte(min(tabulate(wrapped$cluster)), 2L)
  expect_equal(wrapped$parameters$spatial_weight, 0.4, tolerance = 1e-12)
  expect_equal(wrapped$parameters$compactness, 4, tolerance = 1e-12)
  expect_error(
    slice_msf(
      fixture$vec, fixture$mask, num_runs = 1L, consensus = FALSE,
      min_size = 2L, gamma = 1e-6
    ),
    "gamma must be zero",
    class = "slice_msf_unsupported_gamma"
  )
  expect_error(
    cluster4d(
      fixture$vec, fixture$mask, n_clusters = 8L, method = "slice_msf",
      gamma = 1, verbose = FALSE
    ),
    "gamma must be zero",
    class = "slice_msf_unsupported_gamma"
  )
})

test_that("exact-K repair preserves spherical and gradient parcel structure", {
  fixtures <- list(
    spherical = slice_msf_make_spherical_fixture(),
    gradient = slice_msf_make_gradient_fixture()
  )
  minimum_ari <- c(spherical = 0.80, gradient = 0.95)

  for (name in names(fixtures)) {
    fixture <- fixtures[[name]]
    fit <- cluster4d(
      fixture$vec, fixture$mask, n_clusters = 8L, method = "slice_msf",
      spatial_weight = 0.4, min_size = 2L, r = 16L,
      num_runs = 1L, consensus = FALSE, verbose = FALSE
    )
    sizes <- tabulate(fit$cluster)
    repair <- fit$metadata$exact_k_repair

    expect_identical(length(sizes), 8L, info = name)
    expect_true(min(sizes) >= 2L, info = name)
    expect_true(max(sizes) < 0.5 * length(fit$cluster), info = name)
    observed <- clustering_accuracy(fit$cluster, fixture$truth)$ari
    expect_true(
      observed >= minimum_ari[[name]],
      info = slice_msf_quality_diagnostic(
        fit, fixture$truth, minimum_ari[[name]], name,
        fixture$contract$estimand, fixture$contract$seed
      )
    )
    expect_identical(repair$requested_k, 8L, info = name)
    expect_identical(repair$min_cluster_size, 2L, info = name)
    expect_identical(repair$final_sizes, as.integer(sizes), info = name)
    expect_identical(sum(repair$input_sizes), length(fit$cluster), info = name)
    expect_true(is.list(repair$operations), info = name)
  }
})

test_that("axially smoothed sketches are renormalized before cosine distance", {
  fixture <- slice_msf_make_block_fixture()
  fraction <- 0.25
  raw <- slice_msf_single(
    fixture$vec, fixture$mask,
    r = 8L, k = 0.4, min_size = 2L,
    nbhd = 8L, stitch_z = TRUE, gamma = 0, z_mult = 0
  )
  observed <- slice_msf_single(
    fixture$vec, fixture$mask,
    r = 8L, k = 0.4, min_size = 2L,
    nbhd = 8L, stitch_z = TRUE, gamma = 0, z_mult = fraction
  )

  dims <- dim(fixture$mask)
  included <- as.vector(as.array(fixture$mask)) > 0
  expected <- raw$sketch
  for (mode in seq_len(nrow(expected))) {
    for (y in seq_len(dims[[2L]])) {
      for (x in seq_len(dims[[1L]])) {
        global <- x + dims[[1L]] * (y - 1L) +
          dims[[1L]] * dims[[2L]] * (seq_len(dims[[3L]]) - 1L)
        column <- ifelse(included[global], expected[mode, global], 0)
        for (z in seq_len(dims[[3L]])) {
          if (!included[[global[[z]]]]) next
          left <- if (z > 1L) column[[z - 1L]] else column[[z]]
          right <- if (z < dims[[3L]]) column[[z + 1L]] else column[[z]]
          expected[mode, global[[z]]] <-
            (1 - fraction) * column[[z]] + 0.5 * fraction * (left + right)
        }
      }
    }
  }
  expected_norms <- sqrt(colSums(expected[, included, drop = FALSE]^2))
  expected[, included] <- sweep(
    expected[, included, drop = FALSE], 2L,
    pmax(1e-6, expected_norms), `/`
  )
  observed_norms <- sqrt(colSums(observed$sketch[, included, drop = FALSE]^2))

  expect_equal(observed$sketch, expected, tolerance = 1e-12)
  expect_equal(observed_norms, rep(1, sum(included)), tolerance = 1e-12)
  expect_identical(observed$params$feature_distance, "1-cosine")
})

test_that("stitched nbhd 8 uses the same 26-neighbor topology throughout", {
  dims <- c(2L, 2L, 2L)
  n_time <- 12L
  mask_values <- array(0, dims)
  mask_values[1L, 1L, 1L] <- 1
  mask_values[2L, 2L, 2L] <- 1
  signal <- sin(seq(0, 2 * pi, length.out = n_time))
  values <- array(0, c(dims, n_time))
  values[1L, 1L, 1L, ] <- signal
  values[2L, 2L, 2L, ] <- signal
  mask <- neuroim2::NeuroVol(mask_values, neuroim2::NeuroSpace(dims))
  vec <- neuroim2::NeuroVec(values, neuroim2::NeuroSpace(c(dims, n_time)))

  diagonal <- slice_msf(
    vec, mask, target_k_global = -1L,
    r = 4L, compactness = 0, min_size = 1L,
    num_runs = 1L, consensus = FALSE, stitch_z = TRUE, nbhd = 8L
  )
  axial <- slice_msf(
    vec, mask, target_k_global = -1L,
    r = 4L, compactness = 0, min_size = 1L,
    num_runs = 1L, consensus = FALSE, stitch_z = TRUE, nbhd = 4L
  )

  expect_identical(diagonal$actual_k, 1L)
  expect_identical(diagonal$metadata$native$n_components_fh, 1L)
  expect_identical(diagonal$metadata$native$n_components_final, 1L)
  expect_identical(diagonal$metadata$exact_k_repair$natural_k, 1L)
  expect_identical(axial$actual_k, 2L)
})

test_that("adjacent-correlation diagnostics cannot define feature distance", {
  alternating <- rep(c(-1, 1), 8L)
  shifted <- c(alternating[-1L], alternating[[1L]])
  zero_diagnostic <- rep(c(-1, 0, 1, 0), 4L)
  left <- slice_msf_oracle_sketch(alternating, c(1L, 3L, 5L))
  right <- slice_msf_oracle_sketch(shifted, c(1L, 3L, 5L))
  distance <- 1 - max(-1, min(1, sum(left * right)))

  expect_lte(slice_msf_oracle_adjacent_correlation(alternating), 0)
  expect_lte(slice_msf_oracle_adjacent_correlation(zero_diagnostic), 0)
  expect_gt(distance, 0)

  dims <- c(2L, 1L, 1L)
  n_time <- 20L
  phase <- seq(0, 2 * pi, length.out = n_time / 2L)
  first <- as.vector(rbind(sin(phase), -sin(phase)))
  second <- as.vector(rbind(cos(phase), -cos(phase)))
  values <- array(0, c(dims, n_time))
  values[1L, 1L, 1L, ] <- first
  values[2L, 1L, 1L, ] <- second
  mask <- neuroim2::NeuroVol(array(1, dims), neuroim2::NeuroSpace(dims))
  vec <- neuroim2::NeuroVec(values, neuroim2::NeuroSpace(c(dims, n_time)))
  native <- slice_msf_single(
    vec, mask, r = 6L, k = 1e-6, min_size = 1L,
    nbhd = 4L, stitch_z = TRUE, gamma = 0, z_mult = 0
  )
  expect_true(all(native$temporal_smoothness < 0))
  expect_equal(length(unique(native$labels)), 2L)
})
