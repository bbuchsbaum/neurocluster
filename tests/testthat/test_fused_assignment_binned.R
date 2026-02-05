test_that("fused_assignment_binned matches brute-force scoring when all centroids are candidates", {
  if (!exists("fused_assignment_binned")) {
    skip("fused_assignment_binned not available (C++ not built)")
  }

  set.seed(1)
  n <- 4L
  K <- 3L
  D <- 6L

  # coords: 3 x n (columns are voxels)
  coords <- matrix(rnorm(3 * n, sd = 0.05), nrow = 3)

  # keep centroid coords very close so binning captures all clusters easily
  coord_centroids <- matrix(rnorm(3 * K, sd = 0.05), nrow = 3)

  # feature data: D x n, centroids: D x K
  data <- matrix(rnorm(D * n), nrow = D)
  data_centroids <- matrix(rnorm(D * K), nrow = D)

  nn_index <- matrix(-1L, nrow = n, ncol = 1)
  nn_dist <- matrix(0.0, nrow = n, ncol = 1)
  curclus <- rep(0L, n)

  dthresh <- 0
  sigma1 <- 1.3
  sigma2 <- 2.2
  alpha <- 0.6

  # Make candidate set effectively include all centroids
  res <- fused_assignment_binned(
    nn_index = nn_index,
    nn_dist = nn_dist,
    curclus = curclus,
    coords = coords,
    data_centroids = data_centroids,
    coord_centroids = coord_centroids,
    data = data,
    dthresh = dthresh,
    sigma1 = sigma1,
    sigma2 = sigma2,
    alpha = alpha,
    window_factor = 100,
    bin_expand = 5L
  )

  brute <- integer(n)
  inv2sigma1 <- 1 / (2 * sigma1^2)
  inv2sigma2 <- 1 / (2 * sigma2^2)
  for (i in seq_len(n)) {
    vx <- coords[, i]
    x <- data[, i]
    scores <- vapply(seq_len(K), function(k) {
      c <- data_centroids[, k]
      cc <- coord_centroids[, k]
      dd2 <- sum((x - c)^2)
      sd2 <- sum((vx - cc)^2)
      alpha * exp(-dd2 * inv2sigma1) + (1 - alpha) * exp(-sd2 * inv2sigma2)
    }, numeric(1))
    brute[i] <- which.max(scores) - 1L
  }

  expect_identical(as.integer(res), as.integer(brute))
})

