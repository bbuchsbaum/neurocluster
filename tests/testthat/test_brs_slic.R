library(testthat)
library(neuroim2)

test_that("cluster4d brs_slic returns valid structure", {
  if (!exists("brs_slic_core")) {
    skip("Skipping brs_slic test - C++ core not available")
  }

  set.seed(123)
  arr <- array(rnorm(6 * 6 * 6 * 12), dim = c(6, 6, 6, 12))
  vec <- NeuroVec(arr, NeuroSpace(c(6, 6, 6, 12)))
  mask <- NeuroVol(array(1, dim = c(6, 6, 6)), NeuroSpace(c(6, 6, 6)))

  result <- cluster4d(
    vec, mask,
    n_clusters = 8,
    method = "brs_slic",
    embedding_dim = 32,
    max_iterations = 2,
    boundary_passes = 1,
    parallel = FALSE,
    verbose = FALSE
  )

  expect_s3_class(result, "cluster4d_result")
  expect_equal(result$method, "brs_slic")
  expect_length(result$cluster, sum(mask > 0))
  expect_false(any(is.na(result$cluster)))
  expect_gte(result$n_clusters, 1)
  expect_lte(result$n_clusters, 8)
  expect_equal(nrow(result$coord_centers), result$n_clusters)
})

test_that("cluster4d brs_slic is deterministic with fixed seed", {
  if (!exists("brs_slic_core")) {
    skip("Skipping brs_slic test - C++ core not available")
  }

  set.seed(321)
  arr <- array(rnorm(5 * 5 * 5 * 10), dim = c(5, 5, 5, 10))
  vec <- NeuroVec(arr, NeuroSpace(c(5, 5, 5, 10)))
  mask <- NeuroVol(array(1, dim = c(5, 5, 5)), NeuroSpace(c(5, 5, 5)))

  res1 <- cluster4d(
    vec, mask,
    n_clusters = 6,
    method = "brs_slic",
    embedding_dim = 32,
    max_iterations = 2,
    boundary_passes = 1,
    seed = 17,
    parallel = FALSE,
    verbose = FALSE
  )

  res2 <- cluster4d(
    vec, mask,
    n_clusters = 6,
    method = "brs_slic",
    embedding_dim = 32,
    max_iterations = 2,
    boundary_passes = 1,
    seed = 17,
    parallel = FALSE,
    verbose = FALSE
  )

  expect_equal(res1$cluster, res2$cluster)
})
