# Regression tests for issue #3: rena_rnn_coarse_cpp walked each row of the
# column-compressed adjacency with an O(N) row iterator, making ReNA++
# coarsening quadratic in the voxel count. The reference labels below were
# produced by the pre-fix implementation (neurocluster f034a57) and pin that
# the faster neighbour walk returns bit-identical results.

rena_plus_reference_fixture <- function() {
  dims <- c(7L, 6L, 5L)
  g <- as.matrix(expand.grid(seq_len(dims[1]), seq_len(dims[2]), seq_len(dims[3])))
  keep <- rowSums((g - 3.5)^2) <= 14
  n_time <- 6L
  arr <- array(0, c(dims, n_time))
  for (t in seq_len(n_time)) {
    arr[, , , t] <- sin(g[, 1] * t / 3) + cos(g[, 2] * (t + 1) / 4) +
      0.1 * g[, 3] * t +
      0.05 * ((g[, 1] * 7 + g[, 2] * 3 + g[, 3] * 5 + t) %% 11)
  }
  list(
    mask = neuroim2::NeuroVol(
      array(as.numeric(keep), dims),
      neuroim2::NeuroSpace(dims, spacing = c(1, 1, 1))
    ),
    vec = neuroim2::NeuroVec(
      arr, neuroim2::NeuroSpace(c(dims, n_time), spacing = c(1, 1, 1))
    ),
    grad = as.numeric((g[, 1] + 2 * g[, 2] + 3 * g[, 3]) %% 7) / 7
  )
}

test_that("rena_plus reproduces pre-fix reference labels", {
  fixture <- rena_plus_reference_fixture()
  expect_equal(sum(fixture$mask > 0), 164)

  r1 <- rena_plus(fixture$vec, fixture$mask, K = 6, connectivity = 26,
                  lambda = 0)
  expect_identical(as.integer(r1$cluster), c(
    1L, 2L, 1L, 1L, 2L, 3L, 4L, 1L, 2L, 2L, 4L, 4L, 4L,
    4L, 4L, 4L, 4L, 4L, 5L, 5L, 4L, 6L, 5L, 5L, 1L, 1L, 2L, 3L, 4L,
    1L, 1L, 2L, 3L, 3L, 4L, 2L, 2L, 2L, 4L, 4L, 4L, 4L, 4L, 4L, 4L,
    4L, 4L, 5L, 5L, 4L, 6L, 6L, 5L, 5L, 5L, 6L, 4L, 1L, 1L, 2L, 3L,
    3L, 4L, 1L, 1L, 2L, 3L, 3L, 4L, 4L, 2L, 2L, 4L, 4L, 4L, 4L, 4L,
    4L, 4L, 4L, 4L, 4L, 4L, 5L, 5L, 4L, 6L, 6L, 4L, 5L, 5L, 5L, 6L,
    6L, 4L, 1L, 1L, 2L, 4L, 3L, 4L, 1L, 1L, 2L, 3L, 3L, 4L, 4L, 2L,
    2L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 5L, 5L, 4L, 6L,
    6L, 4L, 5L, 5L, 5L, 6L, 6L, 1L, 1L, 2L, 4L, 4L, 1L, 1L, 2L, 4L,
    3L, 4L, 4L, 2L, 2L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 5L, 5L,
    4L, 6L, 6L, 5L, 5L, 5L, 6L
  ))

  r2 <- rena_plus(fixture$vec, fixture$mask, K = 9, connectivity = 6,
                  lambda = 1.5, grad_img = fixture$grad, r = 3)
  expect_identical(as.integer(r2$cluster), c(
    1L, 2L, 3L, 1L, 2L, 4L, 5L,
    3L, 2L, 2L, 4L, 6L, 5L, 5L, 5L, 5L, 4L, 6L, 7L, 7L, 5L, 8L, 7L,
    7L, 3L, 1L, 2L, 4L, 5L, 3L, 1L, 2L, 4L, 6L, 5L, 2L, 2L, 2L, 4L,
    6L, 5L, 5L, 5L, 5L, 4L, 6L, 9L, 7L, 7L, 5L, 8L, 8L, 7L, 7L, 5L,
    8L, 5L, 3L, 1L, 2L, 4L, 6L, 5L, 3L, 1L, 2L, 4L, 6L, 5L, 5L, 2L,
    2L, 4L, 6L, 6L, 5L, 5L, 5L, 5L, 4L, 6L, 6L, 9L, 9L, 7L, 5L, 8L,
    8L, 9L, 7L, 7L, 5L, 8L, 8L, 5L, 3L, 1L, 2L, 4L, 6L, 5L, 3L, 1L,
    2L, 4L, 6L, 5L, 5L, 2L, 2L, 4L, 6L, 6L, 5L, 5L, 5L, 5L, 4L, 6L,
    6L, 9L, 9L, 9L, 5L, 8L, 8L, 9L, 7L, 7L, 5L, 8L, 8L, 3L, 1L, 2L,
    4L, 5L, 3L, 1L, 2L, 4L, 6L, 5L, 5L, 2L, 2L, 4L, 6L, 5L, 5L, 5L,
    5L, 4L, 6L, 9L, 9L, 9L, 5L, 8L, 8L, 7L, 7L, 5L, 8L
  ))
})

test_that("rena coarsening searches row neighbours of an asymmetric graph", {
  # Row i of G (not column i) lists the candidates of node i. In this
  # directed graph the two differ: rows give {1,2} and {3,4,5,6}, while
  # columns would chain every node into one cluster.
  X <- cbind(c(0, 1, 3, 6, 10, 15), c(0, 0.5, 0, 0.5, 0, 0.5))
  G <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4, 5, 6, 2, 6, 1),
    j = c(2, 3, 4, 5, 6, 1, 1, 4, 5),
    x = 1, dims = c(6, 6)
  )
  G <- methods::as(G, "dgCMatrix")
  coarse <- neurocluster:::rena_rnn_coarse_cpp(X, G, numeric(0), 2L, 0, 1L)
  expect_identical(as.integer(coarse$labels_coarse), c(1L, 1L, 2L, 2L, 2L, 2L))
  expect_equal(as.matrix(coarse$G_coarse), matrix(c(0, 1, 1, 0), 2L, 2L))
})

test_that("rena_plus coarsening is near-linear on a large grid", {
  skip_on_cran()
  dims <- c(40L, 40L, 15L)
  space <- neuroim2::NeuroSpace(dims, c(2, 2, 2))
  mask <- neuroim2::NeuroVol(array(1, dims), space)
  set.seed(5)
  vec <- neuroim2::NeuroVec(
    array(stats::rnorm(prod(dims) * 5L), c(dims, 5L)),
    neuroim2::NeuroSpace(c(dims, 5L), c(2, 2, 2))
  )
  # 24,000 voxels: the quadratic row walk took over a minute here; the
  # column walk takes well under a second on a laptop.
  elapsed <- system.time(
    result <- cluster4d(vec, mask, n_clusters = 100L, method = "rena_plus")
  )[["elapsed"]]
  expect_identical(result$n_clusters, 100L)
  expect_lt(elapsed, 30)
})
