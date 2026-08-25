ward_exact_oracle <- function(X, G, sizes, target) {
  X <- as.matrix(X)
  adj <- as.matrix(G) != 0
  adj <- adj | t(adj)
  diag(adj) <- FALSE

  n <- nrow(X)
  active <- rep(TRUE, n)
  means <- lapply(seq_len(n), function(i) X[i, ])
  weights <- as.numeric(sizes)
  members <- lapply(seq_len(n), identity)
  merge_a <- merge_b <- active_counts <- integer()
  merge_cost <- numeric()

  while (sum(active) > target) {
    candidates <- which(
      upper.tri(adj) & adj & outer(active, active, `&`),
      arr.ind = TRUE
    )
    if (!nrow(candidates)) {
      stop("oracle cannot reach target through adjacent merges")
    }

    costs <- apply(candidates, 1L, function(edge) {
      a <- edge[[1L]]
      b <- edge[[2L]]
      weights[[a]] * weights[[b]] / (weights[[a]] + weights[[b]]) *
        sum((means[[a]] - means[[b]])^2)
    })
    pick <- order(costs, candidates[, 1L], candidates[, 2L])[[1L]]
    a <- candidates[pick, 1L]
    b <- candidates[pick, 2L]

    merge_a <- c(merge_a, a)
    merge_b <- c(merge_b, b)
    merge_cost <- c(merge_cost, costs[[pick]])

    total_weight <- weights[[a]] + weights[[b]]
    means[[a]] <- (
      weights[[a]] * means[[a]] + weights[[b]] * means[[b]]
    ) / total_weight
    weights[[a]] <- total_weight
    members[[a]] <- c(members[[a]], members[[b]])
    active[[b]] <- FALSE

    neighbors <- which(active & (adj[a, ] | adj[b, ]))
    neighbors <- setdiff(neighbors, a)
    adj[a, ] <- FALSE
    adj[, a] <- FALSE
    adj[b, ] <- FALSE
    adj[, b] <- FALSE
    adj[a, neighbors] <- TRUE
    adj[neighbors, a] <- TRUE
    active_counts <- c(active_counts, sum(active))
  }

  labels <- integer(n)
  for (cluster in which(active)) {
    labels[members[[cluster]]] <- cluster
  }
  labels <- match(labels, unique(labels))

  list(
    labels_super = labels,
    n_clusters = sum(active),
    merge_a = unname(merge_a),
    merge_b = unname(merge_b),
    merge_cost = unname(merge_cost),
    active_counts = unname(active_counts)
  )
}

ward_candidate <- function(X, G, sizes, target) {
  neurocluster:::ward_on_supervoxels_cpp(
    X_coarse = X,
    G_coarse = G,
    sizes = as.integer(sizes),
    n_clusters = as.integer(target)
  )
}

expect_ward_matches_oracle <- function(X, G, sizes, target) {
  expected <- ward_exact_oracle(X, G, sizes, target)
  observed <- ward_candidate(X, G, sizes, target)

  expect_identical(observed$labels_super, expected$labels_super)
  expect_identical(observed$n_clusters, expected$n_clusters)
  expect_identical(observed$merge_a, expected$merge_a)
  expect_identical(observed$merge_b, expected$merge_b)
  expect_equal(observed$merge_cost, expected$merge_cost, tolerance = 1e-11)
  expect_identical(observed$active_counts, expected$active_counts)
  observed
}

test_that("ReNA++ rejects the six-node stale-Ward counterexample", {
  X <- structure(c(
    -0.841, 1.384, -1.255, 0.070, 1.711, -0.603,
    -0.472, -0.635, -0.286, 0.138, 1.228, -0.802
  ), dim = c(6L, 2L))
  G <- Matrix::sparseMatrix(
    i = c(1:5, 2:6), j = c(2:6, 1:5), x = 1,
    dims = c(6L, 6L)
  )
  sizes <- c(4L, 3L, 2L, 5L, 2L, 2L)

  observed <- expect_ward_matches_oracle(X, G, sizes, target = 2L)

  expect_identical(observed$labels_super, c(1L, 1L, 1L, 1L, 2L, 2L))
  expect_equal(
    observed$merge_cost,
    c(2.764859, 6.908403, 3.337260, 9.475496),
    tolerance = 1e-6
  )
})

test_that("ReNA++ Ward trace agrees with an exact oracle on random graphs", {
  set.seed(20260825)

  for (case in seq_len(60L)) {
    n <- sample(4:9, 1L)
    d <- sample(1:4, 1L)
    X <- matrix(rnorm(n * d), nrow = n)
    sizes <- sample(1:5, n, replace = TRUE)
    target <- sample(seq_len(n), 1L)

    upper <- matrix(FALSE, n, n)
    upper[cbind(seq_len(n - 1L), 2:n)] <- TRUE
    extras <- which(upper.tri(upper, diag = FALSE), arr.ind = TRUE)
    upper[extras] <- upper[extras] | runif(nrow(extras)) < 0.25
    graph <- upper | t(upper)
    G <- methods::as(graph * 1, "dgCMatrix")

    observed <- expect_ward_matches_oracle(X, G, sizes, target)
    expect_true(all(is.finite(observed$merge_cost)))
    expect_true(all(observed$merge_cost >= 0))
    expect_true(all(diff(c(0, cumsum(observed$merge_cost))) >= -1e-12))
    expect_identical(
      observed$active_counts,
      if (target == n) integer() else seq.int(n - 1L, target)
    )
    expect_identical(length(unique(observed$labels_super)), target)
  }
})

test_that("ReNA++ Ward graph handling and ties are deterministic", {
  X <- matrix(c(0, 0, 1, 1, 2, 2, 3, 3), ncol = 2L, byrow = TRUE)
  symmetric <- Matrix::sparseMatrix(
    i = c(1:3, 2:4), j = c(2:4, 1:3), x = 1,
    dims = c(4L, 4L)
  )
  asymmetric_duplicates <- Matrix::sparseMatrix(
    i = rep(1:3, each = 2L), j = rep(2:4, each = 2L), x = 1,
    dims = c(4L, 4L)
  )

  first <- ward_candidate(X, symmetric, rep(1L, 4L), 2L)
  second <- ward_candidate(X, asymmetric_duplicates, rep(1L, 4L), 2L)
  third <- ward_candidate(X, symmetric, rep(1L, 4L), 2L)

  expect_identical(first, second)
  expect_identical(first, third)
  expect_identical(first$merge_a[[1L]], 1L)
  expect_identical(first$merge_b[[1L]], 2L)
})

test_that("ReNA++ Ward actively rejects version-stale queue entries", {
  set.seed(20260825)
  side <- 10L
  index <- matrix(seq_len(side^2), side, side)
  left <- c(as.vector(index[-side, ]), as.vector(index[, -side]))
  right <- c(as.vector(index[-1L, ]), as.vector(index[, -1L]))
  G <- Matrix::sparseMatrix(
    i = c(left, right), j = c(right, left), x = 1,
    dims = rep(side^2, 2L)
  )
  observed <- ward_candidate(
    matrix(rnorm(side^2 * 4L), nrow = side^2),
    G, sample(1:5, side^2, replace = TRUE), 25L
  )

  expect_gt(observed$stale_version_rejections, 0)
  expect_gt(observed$inactive_rejections, 0)
  expect_equal(observed$recomputed_rejections, 0)
  expect_lte(observed$max_queue_size, observed$queue_pushes)
})

test_that("ReNA++ Ward validates its native boundary and fails closed", {
  X <- matrix(seq_len(8), nrow = 4L)
  G <- Matrix::bandSparse(4L, k = c(-1L, 1L), diagonals = list(rep(1, 3), rep(1, 3)))

  expect_error(ward_candidate(X, G, c(1L, 1L), 2L), "sizes")
  expect_error(ward_candidate(X, G, c(1L, 1L, 0L, 1L), 2L), "positive")
  expect_error(ward_candidate(X, G, rep(1L, 4L), 0L), "n_clusters")
  expect_error(ward_candidate(X, G, rep(1L, 4L), 5L), "n_clusters")
  expect_error(ward_candidate(X, G[1:3, 1:3], rep(1L, 4L), 2L), "square")

  bad_X <- X
  bad_X[1L, 1L] <- Inf
  expect_error(ward_candidate(bad_X, G, rep(1L, 4L), 2L), "finite")

  negative_G <- G
  negative_G[1L, 2L] <- -1
  expect_error(ward_candidate(X, negative_G, rep(1L, 4L), 2L), "non-negative")

  disconnected <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L, 4L), j = c(2L, 1L, 4L, 3L), x = 1,
    dims = c(4L, 4L)
  )
  expect_error(
    ward_candidate(X, disconnected, rep(1L, 4L), 1L),
    "cannot reach"
  )
})

test_that("rena_plus validates all method parameters", {
  mask <- neuroim2::NeuroVol(array(1, c(2, 2, 2)), neuroim2::NeuroSpace(c(2, 2, 2)))
  vec <- neuroim2::NeuroVec(
    array(rnorm(2 * 2 * 2 * 3), c(2, 2, 2, 3)),
    neuroim2::NeuroSpace(c(2, 2, 2, 3))
  )

  expect_error(rena_plus(vec, mask, K = 2, r = 0.5), "r must be at least 1")
  expect_error(rena_plus(vec, mask, K = 2, r = Inf), "r must be a finite")
  expect_error(rena_plus(vec, mask, K = 2, lambda = -1), "non-negative")
  expect_error(rena_plus(vec, mask, K = 2, lambda = NA_real_), "lambda must be a finite")
  expect_error(rena_plus(vec, mask, K = 2, connectivity = 27), "6, 18, or 26")
  expect_error(rena_plus(vec, mask, K = 2, connectivity = 6.5), "integer")
  expect_error(rena_plus(vec, mask, K = 2, max_iterations = 0), "positive")
  expect_error(rena_plus(vec, mask, K = 2, max_iterations = 1.5), "integer")
  expect_error(rena_plus(vec, mask, K = 2, grad_img = rep(NA_real_, 8)), "finite")
})
