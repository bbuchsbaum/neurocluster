slic_contract_fixture <- function(dims = c(6L, 6L, 2L), n_time = 3L, seed = 1L) {
  set.seed(seed)
  space <- neuroim2::NeuroSpace(dims)
  list(
    mask = neuroim2::NeuroVol(array(1, dims), space),
    vec = neuroim2::NeuroVec(
      array(rnorm(prod(dims) * n_time), c(dims, n_time)),
      neuroim2::NeuroSpace(c(dims, n_time))
    )
  )
}

expect_slic_labels_connected <- function(labels, mask, connectivity) {
  graph <- neurocluster:::.exact_k_graph(mask, connectivity)
  normalized <- neurocluster:::.exact_k_connected_labels(
    labels, graph$graph, graph$edges
  )
  expect_identical(
    length(unique(normalized)), length(unique(labels))
  )
}

test_that("native SLIC centers summarize its final live labels", {
  fixture <- slic_contract_fixture()
  mask_idx <- which(as.array(fixture$mask) > 0)
  features <- t(as.matrix(neuroim2::series(fixture$vec, mask_idx)))
  features <- scale(features)
  features[!is.finite(features)] <- 0
  coords <- as.matrix(neuroim2::index_to_coord(fixture$mask, mask_idx))

  core <- neurocluster:::slic4d_core(
    features, coords, as.integer(mask_idx - 1L), c(6L, 6L, 2L),
    c(1, 1, 1), 8L, 0, 1L, 1, 1L, "mask_grid",
    TRUE, 0L, 6L, TRUE, FALSE, 0L, numeric(), 0L, FALSE
  )
  counts <- tabulate(core$labels, nbins = core$actual_k)
  expected_features <- rowsum(features, core$labels, reorder = TRUE) / counts
  expected_coords <- rowsum(coords, core$labels, reorder = TRUE) / counts

  expect_true(core$connectivity_changed)
  expect_identical(sort(unique(core$labels)), seq_len(core$actual_k))
  expect_identical(nrow(core$center_feats), core$actual_k)
  expect_identical(nrow(core$center_coords), core$actual_k)
  expect_equal(core$center_feats, unname(expected_features), tolerance = 1e-12)
  expect_equal(core$center_coords, unname(expected_coords), tolerance = 1e-12)
  expect_true(all(core$connectivity_elapsed_ms >= 0))
  expect_true(all(core$center_recompute_elapsed_ms >= 0))
})

test_that("direct SLIC centers and classes use final public labels", {
  fixture <- slic_contract_fixture()
  result <- slic4d_supervoxels(
    fixture$vec, fixture$mask, K = 8, compactness = 0,
    max_iter = 1, n_threads = 1, seed_method = "mask_grid",
    seed_relocate = "none", connectivity = 6,
    strict_connectivity = TRUE, preserve_k = FALSE
  )
  mask_idx <- which(as.array(fixture$mask) > 0)
  features <- t(as.matrix(neuroim2::series(fixture$vec, mask_idx)))
  coords <- as.matrix(neuroim2::index_to_coord(fixture$mask, mask_idx))
  expected <- neurocluster:::compute_cluster_centers(
    result$labels, features, coords
  )

  expect_s3_class(result, "cluster4d_result")
  expect_s3_class(result, "cluster_result")
  expect_equal(result$centers, expected$centers, tolerance = 1e-12)
  expect_equal(result$coord_centers, expected$coord_centers, tolerance = 1e-12)
  expect_identical(nrow(result$centers), result$actual_k)
  expect_identical(nrow(result$coord_centers), result$actual_k)
  expect_slic_labels_connected(result$labels, fixture$mask, 6L)
})

test_that("checkerboard refinement preserves K and strict connectivity", {
  dims <- c(6L, 6L, 1L)
  index <- arrayInd(seq_len(prod(dims)), dims)
  checkerboard <- array((rowSums(index) %% 2L) * 10, dims)
  mask <- neuroim2::NeuroVol(
    array(1, dims), neuroim2::NeuroSpace(dims)
  )
  image <- neuroim2::NeuroVol(checkerboard, neuroim2::space(mask))

  result <- slic4d_supervoxels(
    image, mask, K = 2, compactness = 0, max_iter = 3,
    n_threads = 1, seed_method = "mask_grid", seed_relocate = "none",
    connectivity = 6, strict_connectivity = TRUE, preserve_k = TRUE,
    topup_iters = 2
  )

  expect_true(result$metadata$native_connectivity_changed)
  expect_gt(result$metadata$native_actual_k, 2L)
  expect_identical(result$actual_k, 2L)
  expect_true(result$metadata$preserve_k_feasible)
  expect_slic_labels_connected(result$labels, mask, 6L)
})

test_that("fragmented SLIC fixture is repaired deterministically", {
  fixture <- slic_contract_fixture(c(8L, 8L, 2L), 2L)
  args <- list(
    fixture$vec, fixture$mask, K = 12, compactness = 0,
    max_iter = 1, n_threads = 1, seed_method = "mask_grid",
    seed_relocate = "none", connectivity = 6,
    strict_connectivity = TRUE, preserve_k = TRUE, topup_iters = 2
  )
  first <- do.call(slic4d_supervoxels, args)
  second <- do.call(slic4d_supervoxels, args)

  expect_identical(first$labels, second$labels)
  expect_identical(first$actual_k, 12L)
  expect_gt(first$metadata$native_actual_k, 12L)
  expect_slic_labels_connected(first$labels, fixture$mask, 6L)
})

test_that("strict connectivity takes precedence when preserve-K is infeasible", {
  dims <- c(3L, 1L, 1L)
  mask_values <- array(c(1, 0, 1), dims)
  space <- neuroim2::NeuroSpace(dims)
  mask <- neuroim2::NeuroVol(mask_values, space)
  image <- neuroim2::NeuroVol(array(c(1, 0, 2), dims), space)

  expect_warning(
    result <- slic4d_supervoxels(
      image, mask, K = 1, max_iter = 1, n_threads = 1,
      seed_relocate = "none", connectivity = 6,
      strict_connectivity = TRUE, preserve_k = TRUE
    ),
    "strict connectivity takes precedence"
  )
  expect_false(result$metadata$preserve_k_feasible)
  expect_identical(result$actual_k, 2L)
  expect_slic_labels_connected(result$labels, mask, 6L)

  unconstrained <- suppressWarnings(slic4d_supervoxels(
    image, mask, K = 1, max_iter = 1, n_threads = 1,
    seed_relocate = "none", connectivity = 6,
    strict_connectivity = FALSE, preserve_k = TRUE
  ))
  expect_identical(unconstrained$actual_k, 1L)
})

test_that("SLIC validates finite method and native parameters", {
  fixture <- slic_contract_fixture(c(3L, 3L, 2L), 2L)
  base <- list(fixture$vec, fixture$mask, K = 2, seed_relocate = "none")
  expect_error(do.call(slic4d_supervoxels, c(base, compactness = -1)), "non-negative")
  expect_error(do.call(slic4d_supervoxels, c(base, compactness = Inf)), "finite")
  expect_error(do.call(slic4d_supervoxels, c(base, max_iter = 0)), "positive")
  expect_error(do.call(slic4d_supervoxels, c(base, max_iter = 1.5)), "integer")
  expect_error(do.call(slic4d_supervoxels, c(base, n_threads = -1)), "non-negative")
  expect_error(do.call(slic4d_supervoxels, c(base, step_mm = 0)), "positive")
  expect_error(do.call(slic4d_supervoxels, c(base, n_components = -1)), "non-negative")
  expect_error(do.call(slic4d_supervoxels, c(base, connectivity = 18)), "6 or 26")
  expect_error(do.call(slic4d_supervoxels, c(base, topup_iters = -1)), "non-negative")
  expect_error(do.call(slic4d_supervoxels, c(base, strict_connectivity = NA)), "TRUE or FALSE")

  X <- matrix(seq_len(16), nrow = 8L)
  coords <- cbind(rep(0:1, 4), rep(rep(0:1, each = 2), 2), rep(0:1, each = 4))
  call_core <- function(...) neurocluster:::slic4d_core(
    X, coords, 0:7, c(2L, 2L, 2L), c(1, 1, 1), 2L, ...
  )
  expect_error(call_core(compactness = -1), "compactness")
  expect_error(call_core(max_iter = 0L), "max_iter")
  expect_error(call_core(connectivity = 18L), "connectivity")
  expect_error(call_core(grad_masked = rep(1, 7L)), "grad_masked")
  bad <- X
  bad[1L, 1L] <- NaN
  expect_error(
    neurocluster:::slic4d_core(
      bad, coords, 0:7, c(2L, 2L, 2L), c(1, 1, 1), 2L
    ),
    "finite"
  )
})

test_that("unified SLIC retains the canonical result class", {
  fixture <- slic_contract_fixture(c(4L, 4L, 2L), 3L)
  result <- cluster4d(
    fixture$vec, fixture$mask, n_clusters = 4, method = "slic",
    max_iterations = 1, parallel = FALSE, preserve_k = TRUE,
    seed_relocate = "none", connectivity = 6
  )
  expect_s3_class(result, "cluster4d_result")
  expect_identical(result$actual_k, 4L)
  expect_slic_labels_connected(result$labels, fixture$mask, 6L)
})

test_that("large unified SLIC activates auto parallelism with exact serial parity", {
  fixture <- slic_contract_fixture(c(16L, 16L, 8L), 4L, seed = 20260825L)
  common <- list(
    vec = fixture$vec,
    mask = fixture$mask,
    n_clusters = 32L,
    spatial_weight = 0.4,
    max_iterations = 2L,
    connectivity = 6L,
    seed_method = "mask_grid",
    seed_relocate = "none",
    strict_connectivity = FALSE,
    enforce_connectivity = FALSE,
    preserve_k = FALSE,
    topup_iters = 0L
  )

  thread_env <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)

  sequential <- do.call(cluster4d_slic, c(common, parallel = FALSE))
  expect_identical(
    Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_),
    thread_env
  )
  automatic <- do.call(cluster4d_slic, c(common, parallel = TRUE))
  expect_identical(
    Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_),
    thread_env
  )

  expect_false(sequential$metadata$assignment_parallel_requested)
  expect_false(sequential$metadata$assignment_parallel_used)
  expect_identical(sequential$metadata$native_n_threads, 1L)
  expect_true(automatic$metadata$assignment_parallel_requested)
  expect_true(automatic$metadata$assignment_parallel_used)
  expect_identical(automatic$metadata$native_n_threads, 0L)
  expect_identical(automatic$labels, sequential$labels)
  expect_equal(automatic$centers, sequential$centers, tolerance = 0)
  expect_equal(automatic$coord_centers, sequential$coord_centers, tolerance = 0)
})

test_that("SLIC thread selection never mutates the caller's thread environment", {
  fixture <- slic_contract_fixture(c(5L, 5L, 2L), 3L, seed = 20260826L)
  old_threads <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)
  on.exit({
    if (is.na(old_threads)) {
      Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS")
    } else {
      Sys.setenv(RCPP_PARALLEL_NUM_THREADS = old_threads)
    }
  }, add = TRUE)

  for (initial in list(NA_character_, "3")) {
    if (is.na(initial)) {
      Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS")
    } else {
      Sys.setenv(RCPP_PARALLEL_NUM_THREADS = initial)
    }
    before <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)

    invisible(slic4d_supervoxels(
      fixture$vec, fixture$mask, K = 5L, compactness = 2,
      max_iter = 1L, n_threads = 1L, seed_method = "mask_grid",
      seed_relocate = "none", strict_connectivity = FALSE,
      enforce_connectivity = FALSE
    ))

    expect_identical(
      Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_),
      before
    )
  }
})
