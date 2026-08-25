make_commute_contract_data <- function(dims = c(4L, 3L, 1L),
                                       n_time = 6L,
                                       constant = FALSE) {
  space <- neuroim2::NeuroSpace(dims, c(2, 3, 4))
  mask <- neuroim2::NeuroVol(array(1, dims), space)
  coords <- arrayInd(seq_len(prod(dims)), dims)
  volumes <- lapply(seq_len(n_time), function(time) {
    values <- if (constant) {
      rep(5, nrow(coords))
    } else {
      sin(coords[, 1] / 2 + time / 3) +
        cos(coords[, 2] + time / 5) + coords[, 1] * time / 100
    }
    neuroim2::NeuroVol(array(values, dims), space)
  })
  list(
    mask = mask,
    vec = do.call(neuroim2::concat, volumes),
    dims = dims
  )
}


test_that("bad voxel policies are deterministic, explicit, and provenanced", {
  constant <- matrix(c(
    1, 1, 1, 1,
    1, 2, 3, 4,
    5, 5, 5, 5
  ), nrow = 3L, byrow = TRUE)
  expect_error(
    neurocluster:::.commute_prepare_features(constant, "error"),
    "zero_constant"
  )
  prepared <- neurocluster:::.commute_prepare_features(
    constant, "zero_constant"
  )
  expect_identical(prepared$provenance$constant_voxels, c(1L, 3L))
  expect_equal(prepared$features[c(1L, 3L), ], matrix(0, 2L, 4L))
  expect_false(prepared$provenance$random_noise_injected)

  for (bad in list(NA_real_, Inf, -Inf)) {
    invalid <- constant
    invalid[2L, 2L] <- bad
    expect_error(
      neurocluster:::.commute_prepare_features(
        invalid, "zero_constant"
      ),
      "non-finite"
    )
  }

  fixture <- make_commute_contract_data(constant = TRUE)
  result <- commute_cluster(
    fixture$vec, fixture$mask, K = 3L, ncomp = 3L,
    connectivity = 3L, alpha = 0,
    bad_voxel_policy = "zero_constant"
  )
  expect_identical(
    result$metadata$data_policy$constant_voxels,
    seq_len(prod(fixture$dims))
  )
  expect_false(result$metadata$data_policy$random_noise_injected)
  expect_true(all(is.finite(result$embedding)))
})


test_that("commute clustering never changes or creates global RNG state", {
  fixture <- make_commute_contract_data()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  no_seed <- commute_cluster(
    fixture$vec, fixture$mask, K = 3L, ncomp = 3L,
    connectivity = 3L
  )
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_true(no_seed$metadata$rng$used)
  expect_true(no_seed$metadata$rng$global_state_preserved)

  set.seed(481)
  before <- .Random.seed
  first <- commute_cluster(
    fixture$vec, fixture$mask, K = 3L, ncomp = 3L,
    connectivity = 3L
  )
  expect_identical(.Random.seed, before)
  second <- commute_cluster(
    fixture$vec, fixture$mask, K = 3L, ncomp = 3L,
    connectivity = 3L
  )
  expect_identical(.Random.seed, before)
  expect_identical(first$labels, second$labels)
  expect_equal(first$embedding, second$embedding, tolerance = 0)

  constants <- make_commute_contract_data(constant = TRUE)
  expect_error(
    commute_cluster(
      constants$vec, constants$mask, K = 2L,
      ncomp = 2L, connectivity = 3L
    ),
    "zero_constant"
  )
  expect_identical(.Random.seed, before)
})


test_that("K and ncomp boundaries fail before eigensolver dispatch", {
  fixture <- make_commute_contract_data(c(3L, 2L, 1L))
  n <- prod(fixture$dims)
  expect_error(commute_cluster(fixture$vec, fixture$mask, K = 0L), "positive")
  expect_error(
    commute_cluster(fixture$vec, fixture$mask, K = n + 1L),
    "Cannot create"
  )
  for (bad in list(0L, n, 1.5, NA_real_, Inf)) {
    expect_error(
      commute_cluster(
        fixture$vec, fixture$mask, K = 2L,
        ncomp = bad, connectivity = n - 1L
      ),
      "ncomp"
    )
  }

  singleton <- commute_cluster(
    fixture$vec, fixture$mask, K = n,
    ncomp = n - 1L, connectivity = n - 1L
  )
  expect_identical(singleton$labels, seq_len(n))
  expect_equal(dim(singleton$embedding), c(n, n - 1L))
  expect_true(all(tabulate(singleton$labels, nbins = n) == 1L))
})


test_that("physical kNN graph contract matches an independent oracle", {
  coords <- cbind(c(0, 1, 4, 10), 0, 0)
  features <- matrix(c(
    -1, 0, 1,
    -1, .1, .9,
    .5, 0, -.5,
    1, 0, -1
  ), nrow = 4L, byrow = TRUE)
  graph <- neurocluster:::.commute_knn_graph(
    coords, features, connectivity = 1L,
    alpha = .4, spatial_sigma = 2, feature_sigma = 1,
    weight_mode = "heat"
  )

  distances <- as.matrix(stats::dist(coords))
  diag(distances) <- Inf
  directed <- cbind(seq_len(4L), max.col(-distances, ties.method = "first"))
  expected <- cbind(
    pmin(directed[, 1L], directed[, 2L]),
    pmax(directed[, 1L], directed[, 2L])
  )
  expected <- unique(expected[order(expected[, 1L], expected[, 2L]), , drop = FALSE])
  expect_identical(graph$edges, expected)
  expect_true(Matrix::isSymmetric(graph$graph, checkDN = FALSE))
  expect_true(all(Matrix::diag(graph$graph) == 0))
  expect_match(graph$contract, "physical-coordinate k-nearest")

  # The 1-to-3 edge crosses a two-unit mask-coordinate gap: this method is
  # explicitly proximity-graph clustering, not masked grid adjacency.
  expect_true(any(apply(graph$edges, 1L, identical, c(2L, 3L))))
})


test_that("spectral embedding matches effective-resistance oracle", {
  graph <- Matrix::Matrix(matrix(c(
    0, 2, 0, 1,
    2, 0, 3, 0,
    0, 3, 0, 4,
    1, 0, 4, 0
  ), 4L, 4L, byrow = TRUE), sparse = TRUE)
  spectral <- neurocluster:::.commute_spectral_embedding(
    graph, ncomp = 3L, eigen_tol = 1e-12
  )

  laplacian <- diag(Matrix::rowSums(graph)) - as.matrix(graph)
  oracle <- eigen(laplacian, symmetric = TRUE)
  positive <- which(oracle$values > 1e-12)
  inverse <- oracle$vectors[, positive, drop = FALSE] %*%
    diag(1 / oracle$values[positive], nrow = length(positive)) %*%
    t(oracle$vectors[, positive, drop = FALSE])
  volume <- sum(Matrix::rowSums(graph))
  resistance <- outer(diag(inverse), diag(inverse), "+") - 2 * inverse
  embedded_d2 <- as.matrix(stats::dist(spectral$embedding))^2
  expect_equal(
    unname(embedded_d2), unname(volume * resistance), tolerance = 1e-10
  )

  oracle_subspace <- oracle$vectors[, positive, drop = FALSE]
  expect_equal(
    spectral$eigenvectors %*% t(spectral$eigenvectors),
    oracle_subspace %*% t(oracle_subspace),
    tolerance = 1e-10
  )
})


test_that("centers and result shapes independently summarize final labels", {
  fixture <- make_commute_contract_data()
  result <- commute_cluster(
    fixture$vec, fixture$mask, K = 3L,
    ncomp = 4L, connectivity = 3L
  )
  original <- neurocluster:::.cluster4d_original_data(
    fixture$vec, fixture$mask, "commute oracle"
  )
  ids <- seq_len(3L)
  center_oracle <- do.call(rbind, lapply(ids, function(id) {
    colMeans(original$features[result$labels == id, , drop = FALSE])
  }))
  coord_oracle <- do.call(rbind, lapply(ids, function(id) {
    colMeans(original$coords[result$labels == id, , drop = FALSE])
  }))

  expect_identical(dim(result$centers), c(3L, ncol(original$features)))
  expect_identical(dim(result$coord_centers), c(3L, 3L))
  expect_equal(result$centers, center_oracle, tolerance = 1e-12)
  expect_equal(result$coord_centers, coord_oracle, tolerance = 1e-12)
  expect_identical(dim(result$embedding), c(nrow(original$features), 4L))

  unified <- cluster4d_commute(
    fixture$vec, fixture$mask, n_clusters = 3L,
    ncomp = 4L, connectivity = 3L
  )
  validation <- validate_cluster4d(
    unified, vec = fixture$vec, mask = fixture$mask
  )
  expect_true(validation$valid, info = paste(validation$errors, collapse = "; "))
})
