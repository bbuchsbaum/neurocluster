library(testthat)
library(neuroim2)

test_that("cluster4d corr_slic returns valid structure", {
  set.seed(123)

  arr <- array(rnorm(6 * 6 * 6 * 12), dim = c(6, 6, 6, 12))
  vec <- NeuroVec(arr, NeuroSpace(c(6, 6, 6, 12)))
  mask <- NeuroVol(array(1, dim = c(6, 6, 6)), NeuroSpace(c(6, 6, 6)))

  result <- cluster4d(
    vec,
    mask,
    n_clusters = 8,
    method = "corr_slic",
    max_iterations = 3,
    embedding_dim = 32,
    seed = 42,
    parallel = FALSE,
    verbose = FALSE
  )

  expect_s3_class(result, "cluster4d_result")
  expect_equal(result$method, "corr_slic")
  expect_length(result$cluster, sum(mask > 0))
  expect_false(any(is.na(result$cluster)))
  expect_gte(result$n_clusters, 1)
  expect_lte(result$n_clusters, 8)
})

test_that("cluster4d corr_slic is deterministic with fixed seed", {
  set.seed(321)

  arr <- array(rnorm(5 * 5 * 5 * 10), dim = c(5, 5, 5, 10))
  vec <- NeuroVec(arr, NeuroSpace(c(5, 5, 5, 10)))
  mask <- NeuroVol(array(1, dim = c(5, 5, 5)), NeuroSpace(c(5, 5, 5)))

  res1 <- cluster4d(
    vec,
    mask,
    n_clusters = 6,
    method = "corr_slic",
    max_iterations = 3,
    embedding_dim = 32,
    seed = 7,
    parallel = FALSE,
    verbose = FALSE
  )

  res2 <- cluster4d(
    vec,
    mask,
    n_clusters = 6,
    method = "corr_slic",
    max_iterations = 3,
    embedding_dim = 32,
    seed = 7,
    parallel = FALSE,
    verbose = FALSE
  )

  expect_equal(res1$cluster, res2$cluster)
})

test_that("cluster4d corr_slic supports assignment subsampling stride", {
  set.seed(99)

  arr <- array(rnorm(6 * 6 * 5 * 16), dim = c(6, 6, 5, 16))
  vec <- NeuroVec(arr, NeuroSpace(c(6, 6, 5, 16)))
  mask <- NeuroVol(array(1, dim = c(6, 6, 5)), NeuroSpace(c(6, 6, 5)))

  base <- cluster4d(
    vec,
    mask,
    n_clusters = 7,
    method = "corr_slic",
    max_iterations = 4,
    embedding_dim = 32,
    assign_stride = 1,
    seed = 11,
    parallel = FALSE,
    verbose = FALSE
  )

  stride <- cluster4d(
    vec,
    mask,
    n_clusters = 7,
    method = "corr_slic",
    max_iterations = 4,
    embedding_dim = 32,
    assign_stride = 2,
    seed = 11,
    parallel = FALSE,
    verbose = FALSE
  )

  expect_length(stride$cluster, sum(mask > 0))
  expect_false(any(is.na(stride$cluster)))
  expect_true(stride$n_clusters >= 1)
  expect_true(stride$n_clusters <= 7)

  # Check that subsampling does not catastrophically diverge from baseline.
  ari <- clustering_accuracy(stride$cluster, base$cluster)$ari
  expect_gt(ari, 0.05)
})

test_that("cluster4d corr_slic supports quantized assignment path", {
  set.seed(1234)

  arr <- array(rnorm(7 * 7 * 5 * 18), dim = c(7, 7, 5, 18))
  vec <- NeuroVec(arr, NeuroSpace(c(7, 7, 5, 18)))
  mask <- NeuroVol(array(1, dim = c(7, 7, 5)), NeuroSpace(c(7, 7, 5)))

  base <- cluster4d(
    vec,
    mask,
    n_clusters = 8,
    method = "corr_slic",
    max_iterations = 4,
    embedding_dim = 32,
    quantize_assign = FALSE,
    seed = 21,
    parallel = FALSE,
    verbose = FALSE
  )

  quant <- cluster4d(
    vec,
    mask,
    n_clusters = 8,
    method = "corr_slic",
    max_iterations = 4,
    embedding_dim = 32,
    quantize_assign = TRUE,
    seed = 21,
    parallel = FALSE,
    verbose = FALSE
  )

  expect_length(quant$cluster, sum(mask > 0))
  expect_false(any(is.na(quant$cluster)))
  expect_true(quant$n_clusters >= 1)
  expect_true(quant$n_clusters <= 8)

  # Quantized assignment should remain in the same broad regime as float assignment.
  ari <- clustering_accuracy(quant$cluster, base$cluster)$ari
  expect_gt(ari, 0.05)
})

test_that("cluster4d corr_slic supports adaptive embedding and exact refinement", {
  set.seed(4242)

  arr <- array(rnorm(5 * 5 * 4 * 14), dim = c(5, 5, 4, 14))
  vec <- NeuroVec(arr, NeuroSpace(c(5, 5, 4, 14)))
  mask <- NeuroVol(array(1, dim = c(5, 5, 4)), NeuroSpace(c(5, 5, 4)))

  res <- cluster4d(
    vec,
    mask,
    n_clusters = 6,
    method = "corr_slic",
    max_iterations = 3,
    embedding_dim = "auto",
    adaptive_embedding = TRUE,
    embedding_basis = "dct",
    embedding_whiten = TRUE,
    quantize_assign = TRUE,
    refine_exact_iters = 1,
    refine_boundary_only = TRUE,
    seed = 17,
    parallel = FALSE,
    verbose = FALSE
  )

  expect_length(res$cluster, sum(mask > 0))
  expect_false(any(is.na(res$cluster)))
  expect_true(res$parameters$embedding_dim >= 8)
  expect_identical(res$parameters$adaptive_embedding, TRUE)
  expect_identical(res$parameters$embedding_basis, "dct")
  expect_identical(res$parameters$embedding_whiten, TRUE)
  expect_identical(res$parameters$refine_exact_iters, 1L)
})
