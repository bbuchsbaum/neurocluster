result_contract_fixture <- function(dims = c(3L, 3L, 2L), n_time = 6L) {
  sp <- neuroim2::NeuroSpace(
    dims, spacing = c(1.5, 2, 2.5), origin = c(10, -4, 7)
  )
  voxel <- seq_len(prod(dims))
  frames <- lapply(seq_len(n_time), function(i) {
    values <- array(
      sin(voxel * (i + 1) * 0.23) + cos(voxel * 0.17 - i * 0.31),
      dims
    )
    neuroim2::NeuroVol(values, sp)
  })
  list(
    vec = do.call(neuroim2::concat, frames),
    mask = neuroim2::NeuroVol(array(1, dims), sp)
  )
}

result_contract_args <- list(
  supervoxels = list(max_iterations = 1, parallel = FALSE),
  snic = list(),
  slic = list(
    max_iterations = 1, parallel = FALSE, seed_relocate = "none",
    preserve_k = TRUE
  ),
  corr_slic = list(
    max_iterations = 1, parallel = FALSE, embedding_dim = 8,
    refine_exact_iters = 0
  ),
  brs_slic = list(
    max_iterations = 1, parallel = FALSE, embedding_dim = 8,
    boundary_passes = 0, global_passes = 0
  ),
  slice_msf = list(
    num_runs = 1, consensus = FALSE, stitch_z = TRUE, r = 2
  ),
  flash3d = list(max_iterations = 1, dctM = 4),
  g3s = list(
    max_iterations = 1, n_components = 2,
    use_irlba = FALSE, use_rsvd = FALSE
  ),
  rena = list(max_iterations = 2, connectivity = 6),
  rena_plus = list(max_iterations = 2, connectivity = 6, r = 2),
  mcl = list(max_iterations = 2, connectivity = 6),
  acsc = list(max_iterations = 1, block_size = 2, refine = FALSE),
  commute = list(ncomp = 2)
)

direct_result_contract_args <- result_contract_args
direct_result_contract_args$g3s <- list(
  max_refinement_iter = 1, n_components = 2,
  use_irlba = FALSE, use_rsvd = FALSE
)

direct_result_contract_functions <- list(
  supervoxels = cluster4d_supervoxels,
  snic = cluster4d_snic,
  slic = cluster4d_slic,
  corr_slic = cluster4d_corrslic,
  brs_slic = cluster4d_brsslic,
  slice_msf = cluster4d_slice_msf,
  flash3d = cluster4d_flash3d,
  g3s = cluster4d_g3s,
  rena = cluster4d_rena,
  rena_plus = cluster4d_rena_plus,
  mcl = cluster4d_mcl,
  acsc = cluster4d_acsc,
  commute = cluster4d_commute
)

run_result_contract_methods <- function(fixture) {
  lapply(names(result_contract_args), function(method) {
    args <- c(
      list(
        vec = fixture$vec,
        mask = fixture$mask,
        n_clusters = 2,
        method = method,
        verbose = FALSE
      ),
      result_contract_args[[method]]
    )
    suppressWarnings(do.call(cluster4d, args))
  }) |>
    stats::setNames(names(result_contract_args))
}

run_direct_result_contract_methods <- function(fixture) {
  lapply(names(direct_result_contract_functions), function(method) {
    args <- c(
      list(fixture$vec, fixture$mask, 2),
      direct_result_contract_args[[method]]
    )
    suppressWarnings(do.call(direct_result_contract_functions[[method]], args))
  }) |>
    stats::setNames(names(direct_result_contract_functions))
}

independent_result_centers <- function(labels, vec, mask) {
  mask_idx <- which(as.array(mask) > 0)
  features <- t(as.matrix(neuroim2::series(vec, mask_idx)))
  coords <- as.matrix(neuroim2::index_to_coord(mask, mask_idx))
  ids <- sort(unique(labels))
  list(
    centers = do.call(rbind, lapply(ids, function(id) {
      colMeans(features[labels == id, , drop = FALSE])
    })),
    coord_centers = do.call(rbind, lapply(ids, function(id) {
      colMeans(coords[labels == id, , drop = FALSE])
    }))
  )
}

test_that("every unified method returns the typed cluster4d schema", {
  fixture <- result_contract_fixture()
  results <- run_result_contract_methods(fixture)
  required <- c(
    "labels", "cluster", "clusvol", "centers", "coord_centers",
    "actual_k", "n_clusters", "label_ids", "method", "parameters",
    "provenance"
  )

  for (method in names(results)) {
    result <- results[[method]]
    expect_true(all(required %in% names(result)), info = method)
    expect_identical(result$labels, result$cluster, info = method)
    expect_identical(result$actual_k, result$n_clusters, info = method)
    expect_identical(result$label_ids, seq_len(result$actual_k), info = method)
    expect_identical(sort(unique(result$labels)), result$label_ids, info = method)
    expect_identical(result$method, method, info = method)
    expect_identical(result$provenance$feature_space$representation, "original")
    expect_identical(result$provenance$feature_space$summary, "mean")
    expect_identical(result$provenance$coordinate_space$units, "mm")
    expect_true(validate_cluster4d(result, fixture$vec, fixture$mask)$valid)
  }
})

test_that("every exported cluster4d wrapper returns the same typed schema", {
  fixture <- result_contract_fixture()
  results <- run_direct_result_contract_methods(fixture)
  required <- c(
    "labels", "cluster", "clusvol", "centers", "coord_centers",
    "actual_k", "n_clusters", "label_ids", "method", "parameters",
    "provenance"
  )

  for (method in names(results)) {
    result <- results[[method]]
    expect_true(all(required %in% names(result)), info = method)
    expect_identical(result$method, method, info = method)
    expect_true(
      validate_cluster4d(result, fixture$vec, fixture$mask)$valid,
      info = method
    )
  }
})

test_that("centers and coordinate centers independently summarize final labels", {
  fixture <- result_contract_fixture()
  results <- run_result_contract_methods(fixture)

  for (method in names(results)) {
    result <- results[[method]]
    oracle <- independent_result_centers(
      result$labels, fixture$vec, fixture$mask
    )
    expect_equal(result$centers, oracle$centers, tolerance = 1e-10, info = method)
    expect_equal(
      result$coord_centers, oracle$coord_centers,
      tolerance = 1e-10, info = method
    )
    expect_identical(dim(result$centers), c(result$actual_k, 6L), info = method)
    expect_identical(
      dim(result$coord_centers), c(result$actual_k, 3L), info = method
    )
  }
})

test_that("cluster vector and ClusteredNeuroVol round-trip exactly", {
  fixture <- result_contract_fixture()
  results <- run_result_contract_methods(fixture)

  for (method in names(results)) {
    result <- results[[method]]
    expect_identical(as.integer(result$clusvol@clusters), result$labels, info = method)
    expect_identical(
      which(as.array(result$clusvol@mask) > 0),
      which(as.array(fixture$mask) > 0),
      info = method
    )
    expect_true(
      .cluster4d_same_numeric(
        .cluster4d_geometry(result$clusvol)$affine,
        .cluster4d_geometry(fixture$mask)$affine
      ),
      info = method
    )
  }
})

test_that("validate_cluster4d fails closed on malformed and stale results", {
  fixture <- result_contract_fixture()
  result <- run_result_contract_methods(fixture)$snic

  expect_false(validate_cluster4d(result)$valid)

  malformed <- result
  malformed$provenance <- NULL
  expect_false(validate_cluster4d(malformed, fixture$vec, fixture$mask)$valid)

  malformed <- result
  malformed$centers[1, 1] <- Inf
  expect_false(validate_cluster4d(malformed, fixture$vec, fixture$mask)$valid)

  malformed <- result
  malformed$coord_centers[1, 1] <- NaN
  expect_false(validate_cluster4d(malformed, fixture$vec, fixture$mask)$valid)

  malformed <- result
  malformed$centers <- malformed$centers[, -1, drop = FALSE]
  expect_false(validate_cluster4d(malformed, fixture$vec, fixture$mask)$valid)

  malformed <- result
  malformed$coord_centers <- malformed$coord_centers[, -1, drop = FALSE]
  expect_false(validate_cluster4d(malformed, fixture$vec, fixture$mask)$valid)

  stale <- result
  stale$centers[1, 1] <- stale$centers[1, 1] + 0.25
  expect_false(validate_cluster4d(stale, fixture$vec, fixture$mask)$valid)

  stale <- result
  stale$coord_centers[1, 1] <- stale$coord_centers[1, 1] + 0.25
  expect_false(validate_cluster4d(stale, fixture$vec, fixture$mask)$valid)

  malformed <- result
  malformed$labels[malformed$labels == 2L] <- 3L
  malformed$cluster <- malformed$labels
  expect_false(validate_cluster4d(malformed, fixture$vec, fixture$mask)$valid)

  malformed <- result
  malformed$cluster[1] <- if (malformed$cluster[1] == 1L) 2L else 1L
  expect_false(validate_cluster4d(malformed, fixture$vec, fixture$mask)$valid)

  malformed <- result
  malformed$clusvol@clusters[1] <- if (malformed$labels[1] == 1L) 2L else 1L
  expect_false(validate_cluster4d(malformed, fixture$vec, fixture$mask)$valid)

  malformed <- result
  malformed$provenance$coordinate_space$geometry$affine[1, 4] <-
    malformed$provenance$coordinate_space$geometry$affine[1, 4] + 1
  expect_false(validate_cluster4d(malformed, fixture$vec, fixture$mask)$valid)

  wrong_space <- neuroim2::NeuroSpace(
    dim(fixture$mask), spacing = c(1.5, 2, 2.5), origin = c(11, -4, 7)
  )
  wrong_mask <- neuroim2::NeuroVol(array(1, dim(fixture$mask)), wrong_space)
  expect_false(validate_cluster4d(result, fixture$vec, wrong_mask)$valid)
})
