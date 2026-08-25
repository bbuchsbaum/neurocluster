library(testthat)
library(neuroim2)

corr_brs_fixture <- function(seed = 1901L) {
  set.seed(seed)
  dims <- c(5L, 4L, 3L)
  n <- prod(dims)
  timepoints <- 12L
  features <- matrix(rnorm(n * timepoints), nrow = n, ncol = timepoints)
  mask0 <- 0:(n - 1L)
  coords1 <- cbind(
    mask0 %% dims[1L],
    (mask0 %/% dims[1L]) %% dims[2L],
    mask0 %/% (dims[1L] * dims[2L])
  ) + 1
  list(
    dims = dims,
    n = n,
    timepoints = timepoints,
    features = features,
    mask0 = mask0,
    coords1 = coords1
  )
}

expect_corr_brs_native_summary <- function(core, features, coords1) {
  labels <- as.integer(core$labels)
  ids <- seq_len(as.integer(core$actual_k))
  expect_identical(as.integer(core$label_ids), ids)
  expect_identical(sort(unique(labels)), ids)
  expect_identical(as.integer(core$counts), tabulate(labels, nbins = length(ids)))
  expect_identical(sum(as.integer(core$counts)), nrow(features))

  counts <- as.numeric(tabulate(labels, nbins = length(ids)))
  feature_oracle <- rowsum(features, labels, reorder = TRUE) / counts
  coord_oracle <- rowsum(coords1, labels, reorder = TRUE) / counts
  dimnames(feature_oracle) <- NULL
  dimnames(coord_oracle) <- NULL

  expect_equal(unname(core$original_centers), feature_oracle, tolerance = 1e-12)
  expect_equal(unname(core$centers), feature_oracle, tolerance = 1e-12)
  expect_equal(unname(core$centers_xyz), coord_oracle, tolerance = 1e-12)
  expect_identical(dim(core$original_centers), c(length(ids), ncol(features)))
  expect_identical(dim(core$centers_xyz), c(length(ids), 3L))
}

run_corr_contract_core <- function(fixture, ...) {
  args <- c(list(
    feat = fixture$features,
    mask_lin_idx = fixture$mask0,
    dims = fixture$dims,
    K = 5L,
    d = 8L,
    sketch_repeats = 1L,
    alpha = 0.35,
    max_iter = 4L,
    seed = 31L,
    assign_stride = 1L,
    quantize_assign = FALSE,
    embed_basis = "hash",
    whiten_embed = FALSE,
    refine_exact_iters = 0L,
    refine_boundary_only = TRUE,
    refine_stride = 1L,
    refine_alpha = -1,
    connectivity = 6L,
    min_size = 1L,
    n_threads = 1L,
    verbose = FALSE
  ), list(...))
  do.call(corrslic_core, args[!duplicated(names(args), fromLast = TRUE)])
}

run_brs_contract_core <- function(fixture, ...) {
  args <- c(list(
    feat = fixture$features,
    mask_lin_idx = fixture$mask0,
    dims = fixture$dims,
    K = 5L,
    d = 8L,
    sketch_repeats = 1L,
    alpha = 0.35,
    coarse_iter = 4L,
    boundary_passes = 0L,
    global_passes = 0L,
    refine_spatial = 0.35,
    refine_l2 = 0,
    refine_stride = 1L,
    seed = 31L,
    connectivity = 6L,
    min_size = 1L,
    n_threads = 1L,
    verbose = FALSE
  ), list(...))
  do.call(brs_slic_core, args[!duplicated(names(args), fromLast = TRUE)])
}

test_that("Corr-SLIC native summaries follow final labels in every approximation mode", {
  fixture <- corr_brs_fixture()
  modes <- list(
    coarse = list(),
    refined = list(refine_exact_iters = 2L, refine_boundary_only = FALSE),
    quantized = list(quantize_assign = TRUE),
    strided = list(assign_stride = 3L, refine_exact_iters = 1L, refine_stride = 2L)
  )

  for (mode in names(modes)) {
    core <- do.call(run_corr_contract_core, c(list(fixture = fixture), modes[[mode]]))
    expect_corr_brs_native_summary(core, fixture$features, fixture$coords1)
    expect_identical(ncol(core$embedding_centers), 8L)
    expect_identical(nrow(core$embedding_centers), as.integer(core$actual_k))
    expect_identical(core$params$supports_weight_endpoints, TRUE)
  }
})

test_that("BRS-SLIC native summaries follow coarse and refined final labels", {
  fixture <- corr_brs_fixture(1902L)
  modes <- list(
    coarse = list(),
    boundary_refined = list(boundary_passes = 2L),
    global_refined = list(boundary_passes = 1L, global_passes = 2L),
    strided = list(boundary_passes = 1L, global_passes = 1L, refine_stride = 3L)
  )

  for (mode in names(modes)) {
    core <- do.call(run_brs_contract_core, c(list(fixture = fixture), modes[[mode]]))
    expect_corr_brs_native_summary(core, fixture$features, fixture$coords1)
    expect_identical(ncol(core$coarse_embedding$centers), 8L)
    expect_identical(core$coarse_embedding$stage, "pre_refinement")
    expect_identical(core$params$supports_weight_endpoints, TRUE)
  }
})

test_that("public Corr-SLIC and BRS-SLIC separate original and embedding centers", {
  fixture <- corr_brs_fixture(1903L)
  arr <- array(
    fixture$features,
    dim = c(fixture$dims, fixture$timepoints)
  )
  vec <- NeuroVec(arr, NeuroSpace(c(fixture$dims, fixture$timepoints)))
  mask <- NeuroVol(array(1, dim = fixture$dims), NeuroSpace(fixture$dims))

  corr <- cluster4d_corrslic(
    vec, mask, n_clusters = 5L, embedding_dim = 8L,
    max_iterations = 4L, parallel = FALSE, min_size = 1L,
    refine_exact_iters = 1L, seed = 31L
  )
  brs <- cluster4d_brsslic(
    vec, mask, n_clusters = 5L, embedding_dim = 8L,
    max_iterations = 4L, parallel = FALSE, min_size = 1L,
    boundary_passes = 1L, global_passes = 1L, seed = 31L
  )

  expect_identical(ncol(corr$centers), 12L)
  expect_identical(ncol(corr$metadata$approximation$embedding$centers), 8L)
  expect_identical(ncol(brs$centers), 12L)
  expect_identical(ncol(brs$metadata$approximation$coarse_embedding$centers), 8L)
  expect_identical(corr$metadata$approximation$embedding$row_to_label, corr$label_ids)
  expect_identical(brs$metadata$approximation$native_final$label_ids, brs$label_ids)

  corr_oracle <- rowsum(fixture$features, corr$cluster, reorder = TRUE) /
    as.numeric(tabulate(corr$cluster, nbins = corr$actual_k))
  brs_oracle <- rowsum(fixture$features, brs$cluster, reorder = TRUE) /
    as.numeric(tabulate(brs$cluster, nbins = brs$actual_k))
  dimnames(corr_oracle) <- NULL
  dimnames(brs_oracle) <- NULL
  expect_equal(corr$centers, corr_oracle, tolerance = 1e-12)
  expect_equal(brs$centers, brs_oracle, tolerance = 1e-12)
})

test_that("exact refinement scores match a direct Pearson oracle", {
  fixture <- corr_brs_fixture(1904L)
  features <- fixture$features[1:7, , drop = FALSE]
  prototypes <- fixture$features[8:10, , drop = FALSE]
  coords <- fixture$coords1[1:7, , drop = FALSE]
  center_coords <- fixture$coords1[8:10, , drop = FALSE]
  weight <- 0.37
  scale <- 2.4
  stride <- 2L
  l2_weight <- 0.11

  observed <- corrslic_exact_scores_cpp(
    features, prototypes, coords, center_coords,
    spatial_weight = weight, spatial_scale = scale,
    stride = stride, l2_weight = l2_weight
  )
  keep <- seq.int(1L, ncol(features), by = stride)
  expected <- matrix(NA_real_, nrow(features), nrow(prototypes))
  for (v in seq_len(nrow(features))) {
    for (k in seq_len(nrow(prototypes))) {
      feature_distance <- 1 - stats::cor(features[v, keep], prototypes[k, keep])
      l2 <- mean((features[v, keep] - prototypes[k, keep])^2)
      spatial <- sum((coords[v, ] - center_coords[k, ])^2) / scale^2
      expected[v, k] <- (1 - weight) * feature_distance + l2_weight * l2 + weight * spatial
    }
  }
  expect_equal(observed, expected, tolerance = 1e-8)
})

test_that("spatial blend endpoints disable the opposite score term", {
  fixture <- corr_brs_fixture(1905L)
  x <- fixture$features[1:5, , drop = FALSE]
  p <- fixture$features[6:8, , drop = FALSE]
  xyz <- fixture$coords1[1:5, , drop = FALSE]
  cxyz <- fixture$coords1[6:8, , drop = FALSE]

  feature_only_a <- corrslic_exact_scores_cpp(x, p, xyz, cxyz, 0, 2)
  feature_only_b <- corrslic_exact_scores_cpp(x, p, xyz * 20, cxyz - 11, 0, 2)
  spatial_only_a <- corrslic_exact_scores_cpp(x, p, xyz, cxyz, 1, 2)
  spatial_only_b <- corrslic_exact_scores_cpp(x, p * -3 + 7, xyz, cxyz, 1, 2)
  expect_equal(feature_only_a, feature_only_b, tolerance = 0)
  expect_equal(spatial_only_a, spatial_only_b, tolerance = 0)

  independent <- fixture
  set.seed(777)
  independent$features <- matrix(
    rnorm(length(fixture$features)), nrow(fixture$features), ncol(fixture$features)
  )
  corr_a <- run_corr_contract_core(fixture, alpha = 1, refine_alpha = 1)
  corr_b <- run_corr_contract_core(independent, alpha = 1, refine_alpha = 1)
  brs_a <- run_brs_contract_core(
    fixture, alpha = 1, refine_spatial = 1, refine_l2 = 0,
    boundary_passes = 2L, global_passes = 1L
  )
  brs_b <- run_brs_contract_core(
    independent, alpha = 1, refine_spatial = 1, refine_l2 = 0,
    boundary_passes = 2L, global_passes = 1L
  )
  expect_identical(corr_a$labels, corr_b$labels)
  expect_identical(brs_a$labels, brs_b$labels)
})

test_that("Corr-SLIC and BRS-SLIC reject lossy or out-of-domain parameters", {
  fixture <- corr_brs_fixture(1906L)
  bad_features <- fixture$features
  bad_features[2, 3] <- Inf

  expect_error(run_corr_contract_core(fixture, alpha = 1.01), "0, 1")
  expect_error(run_corr_contract_core(fixture, refine_alpha = 1.01), "0, 1")
  expect_error(run_corr_contract_core(fixture, feat = bad_features), "finite")
  expect_error(
    run_corr_contract_core(fixture, mask_lin_idx = c(fixture$mask0[-1], fixture$mask0[2])),
    "duplicate"
  )
  expect_error(run_brs_contract_core(fixture, refine_spatial = 1.01), "0, 1")
  expect_error(run_brs_contract_core(fixture, refine_l2 = NaN), "finite")

  arr <- array(fixture$features, dim = c(fixture$dims, fixture$timepoints))
  vec <- NeuroVec(arr, NeuroSpace(c(fixture$dims, fixture$timepoints)))
  mask <- NeuroVol(array(1, dim = fixture$dims), NeuroSpace(fixture$dims))
  expect_error(
    cluster4d_corrslic(vec, mask, n_clusters = 5, alpha = 1.01),
    "alpha must be in"
  )
  expect_error(
    cluster4d_brsslic(vec, mask, n_clusters = 5, boundary_passes = 1.5),
    "integer-valued"
  )
})
