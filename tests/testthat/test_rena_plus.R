library(testthat)
library(neurocluster)
library(neuroim2)

set.seed(123)

test_that("rena_plus basic functionality", {
  mask <- NeuroVol(array(1, c(4,4,4)), NeuroSpace(c(4,4,4)))
  vec <- NeuroVec(array(rnorm(4*4*4*6), c(4,4,4,6)),
                  NeuroSpace(c(4,4,4,6)))

  res <- rena_plus(vec, mask, K = 5, r = 3, lambda = 0,
                   connectivity = 6, max_iterations = 20, verbose = FALSE)

  expect_s3_class(res, "cluster4d_result")
  expect_equal(length(res$cluster), sum(mask > 0))
  expect_equal(res$method, "rena_plus")
  expect_equal(res$n_clusters, 5)
  expect_true(all(res$cluster >= 1))
  expect_true(all(!is.na(res$cluster)))
})

test_that("rena_plus uses gradient penalty to separate halves", {
  dims <- c(6,6,4)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))

  data_array <- array(0, c(dims, 8))
  # Use opposing temporal profiles. A constant mean shift would be removed by
  # the documented per-voxel time-series standardization.
  profile <- c(-1, -0.5, 0, 0.5, 1, 0.5, 0, -0.5)
  for (time in seq_along(profile)) {
    data_array[1:3, , , time] <- profile[[time]]
    data_array[4:6, , , time] <- -profile[[time]]
  }

  vec <- NeuroVec(data_array, NeuroSpace(c(dims, 8)))

  grad_vol <- array(0, dims)
  grad_vol[1:3,,] <- 0
  grad_vol[4:6,,] <- 5

  res <- rena_plus(vec, mask, K = 2, r = 2, lambda = 5,
                   grad_img = grad_vol, connectivity = 6, verbose = FALSE)

  labs <- array(res$cluster, dims)
  left <- labs[1:3,,]
  right <- labs[4:6,,]

  maj_left <- as.integer(names(which.max(table(left))))
  maj_right <- as.integer(names(which.max(table(right))))

  acc <- (sum(left == maj_left) + sum(right == maj_right)) / length(labs)

  expect_false(maj_left == maj_right)
  expect_gt(acc, 0.8)
})

test_that("cluster4d dispatches rena_plus", {
  mask <- NeuroVol(array(1, c(3,3,3)), NeuroSpace(c(3,3,3)))
  vec <- NeuroVec(array(rnorm(3*3*3*5), c(3,3,3,5)),
                  NeuroSpace(c(3,3,3,5)))

  res <- cluster4d(vec, mask, n_clusters = 3, method = "rena_plus",
                   r = 2, lambda = 0.5, connectivity = 6, verbose = FALSE)

  expect_s3_class(res, "cluster4d_result")
  expect_equal(res$method, "rena_plus")
  expect_true(res$n_clusters <= 3)
  expect_equal(res$parameters$n_clusters_requested, 3)
  expect_equal(res$parameters$r, 2)
})

test_that("coarsening stops at K_prime and ward returns exact K", {
  mask <- NeuroVol(array(1, c(4,4,3)), NeuroSpace(c(4,4,3)))
  vec <- NeuroVec(array(rnorm(4*4*3*6), c(4,4,3,6)),
                  NeuroSpace(c(4,4,3,6)))

  res <- rena_plus(vec, mask, K = 4, r = 2, lambda = 0,
                   connectivity = 6, max_iterations = 30, verbose = FALSE)

  K_prime <- ceiling(4 * 2)
  # A whole reciprocal-neighbor level is not applied when it would jump below
  # K_prime; Ward then performs exact one-at-a-time adjacent merges.
  expect_gte(res$metadata$coarse_n_clusters, K_prime)
  expect_equal(res$n_clusters, 4)
})
