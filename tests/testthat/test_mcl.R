library(neurocluster)
library(neuroim2)

make_mcl_test_data <- function(dims = c(8, 8, 2), ntime = 24) {
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  coords <- arrayInd(seq_len(nvox), dims)

  t_seq <- seq(0, 2 * pi, length.out = ntime)
  ts_data <- matrix(0, nrow = nvox, ncol = ntime)

  for (v in seq_len(nvox)) {
    x <- coords[v, 1]
    region <- if (x <= dims[1] / 2) 1 else 2
    signal <- if (region == 1) sin(t_seq) else cos(t_seq)
    ts_data[v, ] <- signal + rnorm(ntime, sd = 0.08)
  }

  vols <- lapply(seq_len(ntime), function(t) {
    arr <- array(0, dims)
    arr[mask > 0] <- ts_data[, t]
    NeuroVol(arr, NeuroSpace(dims))
  })

  list(mask = mask, vec = do.call(concat, vols))
}

test_that("cluster4d mcl produces a valid cluster4d_result", {
  set.seed(1)
  dat <- make_mcl_test_data()

  res <- cluster4d(
    dat$vec,
    dat$mask,
    n_clusters = 6,
    method = "mcl",
    max_iterations = 8,
    connectivity = 6,
    inflation = 2.0,
    verbose = FALSE
  )

  expect_s3_class(res, "cluster4d_result")
  expect_equal(res$method, "mcl")
  expect_equal(length(res$cluster), sum(dat$mask > 0))
  expect_true(res$n_clusters >= 2)
  expect_true("iterations" %in% names(res$metadata))
  expect_true("converged" %in% names(res$metadata))
  expect_true(is.logical(res$metadata$converged))
})

test_that("cluster4d mcl can force exact K when requested", {
  set.seed(2)
  dat <- make_mcl_test_data(dims = c(10, 10, 2), ntime = 20)

  target_k <- 5
  res <- cluster4d(
    dat$vec,
    dat$mask,
    n_clusters = target_k,
    method = "mcl",
    exact_k = TRUE,
    max_iterations = 6,
    connectivity = 6,
    inflation = 2.2,
    verbose = FALSE
  )

  expect_equal(res$n_clusters, target_k)
  expect_equal(length(unique(res$cluster)), target_k)
})
