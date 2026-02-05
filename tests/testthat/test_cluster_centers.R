library(testthat)
library(neurocluster)

context("Cluster center computation")

test_that("compute_cluster_centers mean matches manual calculation", {
  set.seed(1)
  labels <- c(1L, 1L, 2L, 2L, 2L, 3L)
  features <- matrix(rnorm(length(labels) * 4), nrow = length(labels), ncol = 4)
  coords <- matrix(rnorm(length(labels) * 3), nrow = length(labels), ncol = 3)

  res <- neurocluster:::compute_cluster_centers(labels, features, coords, method = "mean")

  expect_equal(res$n_clusters, 3)

  manual_centers <- rbind(
    colMeans(features[labels == 1L, , drop = FALSE]),
    colMeans(features[labels == 2L, , drop = FALSE]),
    colMeans(features[labels == 3L, , drop = FALSE])
  )
  manual_coords <- rbind(
    colMeans(coords[labels == 1L, , drop = FALSE]),
    colMeans(coords[labels == 2L, , drop = FALSE]),
    colMeans(coords[labels == 3L, , drop = FALSE])
  )

  expect_equal(res$centers, manual_centers, tolerance = 1e-12)
  expect_equal(res$coord_centers, manual_coords, tolerance = 1e-12)
})

