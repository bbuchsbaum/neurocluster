slice_msf_contract_fixture <- function(dims = c(8L, 8L, 2L), n_time = 32L,
                                       seed = 7L, negative_mask = FALSE) {
  set.seed(seed)
  tt <- 0:(n_time - 1L)
  basis <- vapply(
    seq_len(15L),
    function(k) cos(pi * (tt + 0.5) * k / n_time),
    numeric(n_time)
  )
  values <- array(0, c(dims, n_time))
  for (z in seq_len(dims[3L])) {
    for (y in seq_len(dims[2L])) {
      for (x in seq_len(dims[1L])) {
        region <- 1L + (x - 1L) %/% 2L + 4L * ((y - 1L) %/% 4L)
        coefficients <- stats::rnorm(15L, sd = 0.08)
        coefficients[region] <- 1
        coefficients[min(15L, region + 4L)] <- 0.5
        values[x, y, z, ] <-
          as.vector(basis %*% coefficients) + stats::rnorm(n_time, sd = 0.02)
      }
    }
  }
  mask_values <- array(1, dims)
  if (negative_mask) mask_values[1L, 1L, 1L] <- -1
  list(
    vec = neuroim2::NeuroVec(
      values, neuroim2::NeuroSpace(c(dims, n_time))
    ),
    mask = neuroim2::NeuroVol(
      mask_values, neuroim2::NeuroSpace(dims)
    ),
    dims = dims,
    n_time = n_time
  )
}

slice_msf_all_connected <- function(cluster_map, connectivity = 6L) {
  offsets <- rbind(
    c(-1L, 0L, 0L), c(1L, 0L, 0L),
    c(0L, -1L, 0L), c(0L, 1L, 0L),
    c(0L, 0L, -1L), c(0L, 0L, 1L)
  )
  if (connectivity == 26L) {
    offsets <- as.matrix(expand.grid(-1L:1L, -1L:1L, -1L:1L))
    offsets <- offsets[rowSums(abs(offsets)) > 0L, , drop = FALSE]
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
        candidates[, 1L] <= dims[1L] &
        candidates[, 2L] <= dims[2L] &
        candidates[, 3L] <= dims[3L]
      if (any(valid)) queue <- rbind(queue, candidates[valid, , drop = FALSE])
    }
    visited == nrow(members)
  }, logical(1L)))
}

test_that("Slice-MSF ensemble is diverse, seeded, and gamma fails closed", {
  fixture <- slice_msf_contract_fixture()
  fit <- slice_msf(
    fixture$vec, fixture$mask,
    r = 6L, compactness = 10, min_size = 2L,
    num_runs = 3L, consensus = TRUE, use_features = TRUE,
    seed = 11L, ensemble_fraction = 0.6
  )
  replay <- slice_msf(
    fixture$vec, fixture$mask,
    r = 6L, compactness = 10, min_size = 2L,
    num_runs = 3L, consensus = TRUE, use_features = TRUE,
    seed = 11L, ensemble_fraction = 0.6
  )
  signatures <- vapply(
    fit$runs, function(run) paste(run$labels, collapse = ","), character(1L)
  )

  expect_gte(length(unique(signatures)), 3L)
  expect_identical(fit$cluster, replay$cluster)
  expect_identical(fit$metadata$ensemble$runs, replay$metadata$ensemble$runs)

  gamma_zero <- slice_msf(
    fixture$vec, fixture$mask,
    r = 8L, compactness = 10, min_size = 2L,
    num_runs = 1L, consensus = FALSE, gamma = 0
  )
  expect_equal(gamma_zero$metadata$native$gamma, 0)
  expect_error(
    slice_msf(
      fixture$vec, fixture$mask,
      r = 8L, compactness = 10, min_size = 2L,
      num_runs = 1L, consensus = FALSE, gamma = 4
    ),
    "gamma must be zero",
    class = "slice_msf_unsupported_gamma"
  )
})

test_that("Slice-MSF exact K splits, merges, and preserves 3-D connectivity", {
  fixture <- slice_msf_contract_fixture()
  natural <- slice_msf(
    fixture$vec, fixture$mask,
    r = 8L, compactness = 10, min_size = 2L,
    num_runs = 1L, consensus = FALSE
  )
  merged <- slice_msf(
    fixture$vec, fixture$mask, target_k_global = 3L,
    r = 8L, compactness = 10, min_size = 2L,
    num_runs = 1L, consensus = FALSE
  )
  expect_gt(length(unique(natural$cluster)), 3L)
  expect_equal(length(unique(merged$cluster)), 3L)
  expect_true(slice_msf_all_connected(merged$cluster_map, 26L))

  dims <- c(6L, 6L, 1L)
  n_time <- 16L
  signal <- sin(seq(0, 2 * pi, length.out = n_time))
  constant_values <- array(rep(signal, each = prod(dims)), c(dims, n_time))
  constant_vec <- neuroim2::NeuroVec(
    constant_values, neuroim2::NeuroSpace(c(dims, n_time))
  )
  constant_mask <- neuroim2::NeuroVol(
    array(1, dims), neuroim2::NeuroSpace(dims)
  )
  under <- slice_msf(
    constant_vec, constant_mask,
    r = 4L, compactness = 0, min_size = 1L,
    num_runs = 1L, consensus = FALSE
  )
  split <- slice_msf(
    constant_vec, constant_mask, target_k_global = 5L,
    r = 4L, compactness = 0, min_size = 1L,
    num_runs = 1L, consensus = FALSE
  )
  expect_lt(length(unique(under$cluster)), 5L)
  expect_equal(length(unique(split$cluster)), 5L)
  expect_true(slice_msf_all_connected(split$cluster_map, 26L))
})

test_that("Slice-MSF applies mask and result contracts to consensus", {
  fixture <- slice_msf_contract_fixture(negative_mask = TRUE)
  fit <- slice_msf(
    fixture$vec, fixture$mask, target_k_global = 7L,
    r = 6L, compactness = 10, min_size = 2L,
    num_runs = 3L, consensus = TRUE, use_features = TRUE,
    seed = 29L, ensemble_fraction = 0.6
  )

  expect_equal(length(fit$cluster), prod(fixture$dims) - 1L)
  expect_equal(fit$cluster_map[1L, 1L, 1L], 0L)
  expect_equal(dim(fit$centers), c(7L, fixture$n_time))
  expect_equal(dim(fit$coord_centers), c(7L, 3L))
  expect_equal(fit$parameters$seed, 29L)
  expect_equal(fit$metadata$ensemble$subspace_fraction, 0.6)
  expect_length(fit$metadata$ensemble$runs, 3L)
  expect_equal(fit$metadata$consensus$stitch_z, TRUE)
  expect_equal(fit$metadata$topology$exact_k_engine,
               "shared_adjacency_preserving")
  expect_true(slice_msf_all_connected(fit$cluster_map, 26L))
})

test_that("retained Slice-MSF controls reach their declared intermediates", {
  fixture <- slice_msf_contract_fixture()
  base <- slice_msf_single(
    fixture$vec, fixture$mask,
    r = 4L, k = 0.2, min_size = 2L, nbhd = 6L,
    stitch_z = TRUE, gamma = 0, z_mult = 0
  )
  changed <- slice_msf_single(
    fixture$vec, fixture$mask,
    r = 4L, k = 0.4, min_size = 3L, nbhd = 8L,
    stitch_z = FALSE, gamma = 0, z_mult = 0.25,
    dct_frequencies = c(1L, 3L, 5L, 7L),
    dct_weights = c(1, 2, 1, 0.5)
  )
  expect_equal(nrow(base$sketch), 4L)
  expect_equal(changed$params$fh_scale, 0.4)
  expect_equal(changed$params$min_size, 3L)
  expect_equal(changed$params$nbhd, 8L)
  expect_false(changed$params$stitch_z)
  expect_equal(changed$params$gamma, 0)
  expect_equal(changed$params$z_mult, 0.25)
  expect_equal(changed$params$dct_frequencies, c(1L, 3L, 5L, 7L))
  expect_false(isTRUE(all.equal(base$sketch, changed$sketch)))

  fused <- slice_msf_consensus(
    list(base, changed), fixture$mask,
    nbhd = 6L, k_fuse = 0.25, min_size_fuse = 2L,
    use_features = TRUE, lambda = 0.35, stitch_z = TRUE
  )
  expect_equal(
    unname(fused$params[c("nbhd", "fh_scale", "min_size", "use_features",
                          "lambda", "stitch_z")]),
    list(4L, 0.25, 2L, TRUE, 0.35, TRUE)
  )
  expect_equal(fused$params$run_weighting, "uniform")
  altered_diagnostics <- list(base, changed)
  altered_diagnostics[[1L]]$temporal_smoothness[] <- -100
  altered_diagnostics[[2L]]$temporal_smoothness[] <- 100
  fused_altered <- slice_msf_consensus(
    altered_diagnostics, fixture$mask,
    nbhd = 6L, k_fuse = 0.25, min_size_fuse = 2L,
    use_features = TRUE, lambda = 0.35, stitch_z = TRUE
  )
  expect_identical(fused$labels, fused_altered$labels)
  expect_error(
    slice_msf_consensus(
      list(base, changed), fixture$mask, min_size_fuse = 2L, gamma = 1
    ),
    "gamma must be zero",
    class = "slice_msf_unsupported_gamma"
  )

  seed_a <- slice_msf(
    fixture$vec, fixture$mask,
    r = 6L, compactness = 10, min_size = 2L,
    num_runs = 2L, consensus = TRUE, seed = 1L,
    ensemble_fraction = 1
  )
  seed_b <- slice_msf(
    fixture$vec, fixture$mask,
    r = 6L, compactness = 10, min_size = 2L,
    num_runs = 2L, consensus = TRUE, seed = 2L,
    ensemble_fraction = 0.5
  )
  expect_false(identical(
    seed_a$metadata$ensemble$runs,
    seed_b$metadata$ensemble$runs
  ))
  expect_error(
    slice_msf(
      fixture$vec, fixture$mask,
      num_runs = 2L, consensus = FALSE, min_size = 2L
    ),
    "requires consensus"
  )
  expect_error(
    slice_msf(
      fixture$vec, fixture$mask, target_k_global = 2L,
      num_runs = 1L, consensus = FALSE, stitch_z = FALSE,
      min_size = 2L
    ),
    "requires stitch_z"
  )
})
