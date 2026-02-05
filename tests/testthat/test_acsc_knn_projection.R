test_that("build_acsc_graph random projection does not perturb RNG state", {
  # Internal helper used by acsc(); safe to call directly in package tests.
  if (!exists("build_acsc_graph", mode = "function")) {
    skip("build_acsc_graph not available")
  }

  nb <- 12L
  T <- 40L
  block_summary <- list(
    summaries = matrix(rnorm(nb * T), nrow = nb),
    spatial = matrix(rnorm(nb * 3), nrow = nb)
  )

  set.seed(42)
  expected <- runif(5)

  set.seed(42)
  g <- build_acsc_graph(
    block_summary,
    ann_k = 3,
    alpha = 0.5,
    spatial_weighting = "binary",
    block_size = 2,
    knn_proj_dim = 10L,
    knn_proj_method = "rp",
    knn_proj_seed = 999L
  )
  actual <- runif(5)

  expect_equal(actual, expected)
  expect_true(igraph::is_igraph(g))
  expect_equal(igraph::vcount(g), nb)
})

test_that("build_acsc_graph supports DCT projection when available", {
  if (!exists("build_acsc_graph", mode = "function")) {
    skip("build_acsc_graph not available")
  }
  if (!exists("make_dct_basis", mode = "function")) {
    skip("make_dct_basis not available (C++ not built)")
  }

  nb <- 10L
  T <- 30L
  block_summary <- list(
    summaries = matrix(rnorm(nb * T), nrow = nb),
    spatial = matrix(rnorm(nb * 3), nrow = nb)
  )

  g <- build_acsc_graph(
    block_summary,
    ann_k = 3,
    alpha = 0.5,
    spatial_weighting = "gaussian",
    block_size = 2,
    knn_proj_dim = 8L,
    knn_proj_method = "dct"
  )

  expect_true(igraph::is_igraph(g))
  expect_equal(igraph::vcount(g), nb)
})

