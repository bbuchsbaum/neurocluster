library(testthat)
library(neuroim2)

.flash3d_slow_score_labels <- function(core, dims, mask_idx, lambda,
                                       vox_scale, barrier) {
  coords <- arrayInd(mask_idx, dims) - 1
  n <- length(mask_idx)
  k <- nrow(core$assignment_site_coords)
  bits <- ncol(core$voxel_hash_bits)
  scaled <- matrix(vox_scale, nrow = n, ncol = 3L, byrow = TRUE)
  scores <- matrix(NA_real_, nrow = n, ncol = k)

  for (cluster in seq_len(k)) {
    site_bits <- matrix(
      core$assignment_site_hash_bits[cluster, ],
      nrow = n, ncol = bits, byrow = TRUE
    )
    hamming <- rowSums(abs(core$voxel_hash_bits - site_bits)) / bits
    spatial <- rowSums(
      (sweep(coords, 2L, core$assignment_site_coords[cluster, ], "-") * scaled)^2
    ) / core$S2
    barrier_cost <- abs(
      barrier[mask_idx] - barrier[core$assignment_site_seed_grid[cluster]]
    )
    scores[, cluster] <- core$lambda_s_final * spatial +
      lambda[[2L]] * hamming + lambda[[3L]] * barrier_cost
  }
  max.col(-scores, ties.method = "first")
}

test_that("FLASH-3D boundary field is candidate-dependent", {
  dims <- c(9L, 1L, 1L)
  ts <- rbind(
    seq_len(9L), rev(seq_len(9L)),
    rep(c(0, 1, 0), 3L), rep(c(1, 0, 1), 3L)
  )
  lambda <- c(0.2, 1, 50)
  zero <- flash3d_supervoxels_cpp(
    ts, 1:9, dims, 2L, lambda, 1L, 64L, 4L,
    c(1, 1, 1), rep(0, 9), FALSE
  )
  boundary <- flash3d_supervoxels_cpp(
    ts, 1:9, dims, 2L, lambda, 1L, 64L, 4L,
    c(1, 1, 1), c(rep(0, 4), rep(10, 5)), FALSE
  )

  expect_identical(as.integer(zero$labels), c(1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L))
  expect_identical(as.integer(boundary$labels), c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L))
  expect_false(identical(zero$labels, boundary$labels))
})

test_that("one-frame FLASH-3D completes in a subprocess with expected labels", {
  skip_if_not_installed("processx")
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  package_path <- system.file(package = "neurocluster")
  installed_test <- file.exists(file.path(package_path, "Meta", "package.rds"))
  load_code <- if (installed_test) {
    paths <- paste(
      encodeString(.libPaths(), quote = '"'), collapse = ","
    )
    paste0(".libPaths(c(", paths, ")); library(neurocluster); ")
  } else {
    paste0(
      "devtools::load_all(", encodeString(root, quote = '"'),
      ", quiet=TRUE, recompile=FALSE); "
    )
  }
  code <- paste0(
    load_code,
    "out <- neurocluster:::flash3d_supervoxels_cpp(matrix(1:9,nrow=1),1:9,",
    "c(9L,1L,1L),2L,c(.5,1,0),1L,64L,4L,c(1,1,1),NULL,FALSE); ",
    "stopifnot(identical(as.integer(out$labels),c(1L,1L,1L,2L,2L,2L,2L,2L,2L))); ",
    "cat('FLASH_ONE_FRAME_OK')"
  )
  process <- processx::run(
    file.path(R.home("bin"), "Rscript"), c("-e", code),
    timeout = 5000, error_on_status = FALSE,
    env = c(LC_ALL = "en_US.UTF-8", LANG = "en_US.UTF-8")
  )
  expect_false(process$timeout)
  expect_identical(process$status, 0L, info = process$stderr)
  expect_match(process$stdout, "FLASH_ONE_FRAME_OK", info = process$stderr)
})

test_that("FLASH-3D reseeds empty clusters with distinct valid donors", {
  repaired <- flash3d_reseed_empty_cpp(
    c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L),
    seq_len(8L), K = 4L
  )
  expect_identical(sort(unique(repaired$labels)), 1:4)
  expect_length(repaired$reseed_indices, 2L)
  expect_length(unique(repaired$reseed_indices), 2L)
  expect_true(all(repaired$reseed_indices %in% seq_len(8L)))
  expect_true(all(tabulate(repaired$labels, nbins = 4L) > 0L))

  set.seed(1)
  core <- flash3d_supervoxels_cpp(
    matrix(rnorm(6 * 16), nrow = 6L), 1:16, c(4L, 4L, 1L),
    12L, c(0, 5, 0), 1L, 64L, 6L, c(1, 1, 1), NULL, FALSE, TRUE
  )
  expect_identical(sort(unique(core$labels)), 1:12)
  expect_true(all(tabulate(core$labels, nbins = 12L) > 0L))
})

test_that("FLASH-3D scoring agrees with an exhaustive slow oracle", {
  dims <- c(7L, 5L, 1L)
  n <- prod(dims)
  set.seed(42)
  ts <- matrix(rnorm(10L * n), nrow = 10L)
  barrier <- rep(c(0, 2), length.out = n)
  lambda <- c(0.35, 0.65, 1.2)
  vox_scale <- c(1, 1.5, 2)

  for (bits in c(64L, 128L)) {
    core <- flash3d_supervoxels_cpp(
      ts, seq_len(n), dims, 4L, lambda, 1L, bits, 8L,
      vox_scale, barrier, FALSE, TRUE
    )
    oracle <- .flash3d_slow_score_labels(
      core, dims, seq_len(n), lambda, vox_scale, barrier
    )
    expect_identical(as.integer(core$labels), as.integer(oracle), info = bits)
  }
})

test_that("public FLASH-3D centers are final-label means in physical space", {
  dims <- c(5L, 4L, 2L)
  spatial <- NeuroSpace(dims, spacing = c(2, 3, 4))
  mask <- NeuroVol(array(1, dims), spatial)
  values <- array(seq_len(prod(dims)), dims)
  vec <- NeuroVol(values, spatial)
  result <- supervoxels_flash3d(vec, mask, K = 5L, rounds = 1L)

  mask_idx <- which(mask > 0)
  coords <- as.matrix(index_to_coord(mask, mask_idx))
  features <- matrix(as.numeric(values), ncol = 1L)
  expected_coords <- rowsum(coords, result$labels) /
    as.numeric(tabulate(result$labels, nbins = 5L))
  expected_centers <- rowsum(features, result$labels) /
    as.numeric(tabulate(result$labels, nbins = 5L))

  expect_equal(result$coord_centers, unname(expected_coords), tolerance = 1e-12)
  expect_equal(result$centers, unname(expected_centers), tolerance = 1e-12)
  expect_identical(result$actual_k, 5L)
  expect_identical(result$parameters$actual_K, 5L)
  expect_false(isTRUE(all.equal(
    result$coord_centers, result$metadata$cpp_coord_centers_voxel
  )))
})

test_that("FLASH-3D validates every native scoring parameter", {
  dims <- c(3L, 3L, 1L)
  spatial <- NeuroSpace(dims)
  mask <- NeuroVol(array(1, dims), spatial)
  vec <- NeuroVol(array(seq_len(prod(dims)), dims), spatial)

  expect_error(supervoxels_flash3d(vec, mask, 2L, lambda_s = Inf), "finite")
  expect_error(supervoxels_flash3d(vec, mask, 2L, lambda_t = -1), "non-negative")
  expect_error(supervoxels_flash3d(vec, mask, 2L, lambda_g = 1), "requires barrier")
  expect_error(supervoxels_flash3d(vec, mask, 2L, rounds = 1.5), "integer")
  expect_error(supervoxels_flash3d(vec, mask, 2L, bits = 64.5), "integer")
  expect_error(supervoxels_flash3d(vec, mask, 2L, dctM = NA), "finite")
  expect_error(supervoxels_flash3d(vec, mask, 2L, vox_scale = c(1, 0, 1)), "positive")
  expect_error(
    supervoxels_flash3d(vec, mask, 2L, lambda_g = 1,
                        barrier = array(-1, dims)),
    "non-negative"
  )
})

test_that("unified FLASH-3D preserves spatial-weight endpoints", {
  dims <- c(4L, 4L, 1L)
  spatial <- NeuroSpace(dims)
  mask <- NeuroVol(array(1, dims), spatial)
  vec <- NeuroVol(array(seq_len(prod(dims)), dims), spatial)

  feature_only <- cluster4d_flash3d(
    vec, mask, n_clusters = 3L, spatial_weight = 0,
    max_iterations = 1L, bits = 64L, dctM = 4L
  )
  spatial_only <- cluster4d_flash3d(
    vec, mask, n_clusters = 3L, spatial_weight = 1,
    max_iterations = 1L, bits = 64L, dctM = 4L
  )
  expect_identical(feature_only$parameters$lambda_s, 0)
  expect_identical(feature_only$parameters$lambda_t, 1)
  expect_identical(spatial_only$parameters$lambda_s, 1)
  expect_identical(spatial_only$parameters$lambda_t, 0)
  expect_true(feature_only$parameters$normalized_weight_mixture)
  expect_true(spatial_only$parameters$normalized_weight_mixture)
})
