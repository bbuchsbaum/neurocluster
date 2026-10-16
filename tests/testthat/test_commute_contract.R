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

test_that("sparse and dense eigensolvers give the same embedding and labels", {
  fixture <- make_commute_contract_data(dims = c(9L, 8L, 6L), n_time = 7L)
  for (mode in c("heat", "binary")) {
    dense <- commute_cluster(
      fixture$vec, fixture$mask, K = 12L, ncomp = 7L,
      connectivity = 6L, weight_mode = mode, eigen_solver = "dense"
    )
    sparse <- commute_cluster(
      fixture$vec, fixture$mask, K = 12L, ncomp = 7L,
      connectivity = 6L, weight_mode = mode, eigen_solver = "sparse"
    )
    expect_identical(dense$parameters$eigen_solver, "dense")
    expect_identical(sparse$parameters$eigen_solver, "sparse")
    expect_identical(sparse$metadata$spectral$solver, "sparse")
    expect_equal(sparse$eigenvalues, dense$eigenvalues, tolerance = 1e-10)
    # Eigenvector signs are arbitrary; commute-time distances are not.
    expect_equal(
      as.matrix(stats::dist(sparse$embedding)),
      as.matrix(stats::dist(dense$embedding)),
      tolerance = 1e-8
    )
    expect_identical(sparse$labels, dense$labels)
  }
})


test_that("sparse spectral embedding matches the effective-resistance oracle", {
  set.seed(11)
  n <- 30L
  edges <- rbind(cbind(seq_len(n - 1L), 2:n), cbind(1:10, 21:30))
  graph <- Matrix::sparseMatrix(
    i = c(edges[, 1L], edges[, 2L]), j = c(edges[, 2L], edges[, 1L]),
    x = rep(stats::runif(nrow(edges), 0.5, 2), 2L), dims = c(n, n)
  )
  dense <- neurocluster:::.commute_spectral_embedding(
    graph, ncomp = 5L, eigen_solver = "dense"
  )
  sparse <- neurocluster:::.commute_spectral_embedding(
    graph, ncomp = 5L, eigen_solver = "sparse"
  )
  expect_identical(sparse$solver, "sparse")
  expect_equal(sparse$eigenvalues, dense$eigenvalues, tolerance = 1e-10)
  expect_equal(
    sparse$eigenvectors %*% t(sparse$eigenvectors),
    dense$eigenvectors %*% t(dense$eigenvectors),
    tolerance = 1e-8
  )
  expect_identical(sparse$zero_eigenvalues, 1L)
  expect_s4_class(sparse$laplacian, "sparseMatrix")
  residual <- sparse$laplacian %*% sparse$eigenvectors -
    sweep(sparse$eigenvectors, 2L, sparse$eigenvalues, "*")
  expect_lt(max(sqrt(Matrix::colSums(residual^2))), 1e-8)

  # Asking for the whole spectrum is the dense problem.
  full <- neurocluster:::.commute_spectral_embedding(
    graph, ncomp = n - 1L, eigen_solver = "sparse"
  )
  expect_identical(full$solver, "dense")
})


test_that("auto solver keeps small problems dense and large ones sparse", {
  expect_identical(neurocluster:::.commute_resolve_solver("auto", 3000L), "dense")
  expect_identical(neurocluster:::.commute_resolve_solver("dense", 10L), "dense")
  expect_message(
    neurocluster:::.commute_resolve_solver("dense", 3001L),
    "dense eigendecomposition"
  )
  expect_identical(
    neurocluster:::.commute_resolve_solver("auto", 3001L), "sparse"
  )
  expect_error(
    neurocluster:::.commute_resolve_solver("eigs", 10L),
    "should be one of"
  )

  fixture <- make_commute_contract_data()
  small <- commute_cluster(fixture$vec, fixture$mask, K = 3L, ncomp = 4L,
                           connectivity = 3L)
  expect_identical(small$parameters$eigen_solver, "dense")
})


test_that("commute scales to a whole-brain-sized mask with the sparse solver", {
  skip_on_cran()
  dims <- c(30L, 30L, 20L)
  space <- neuroim2::NeuroSpace(dims, c(2, 2, 2))
  mask <- neuroim2::NeuroVol(array(1, dims), space)
  set.seed(3)
  vec <- neuroim2::NeuroVec(
    array(stats::rnorm(prod(dims) * 5L), c(dims, 5L)),
    neuroim2::NeuroSpace(c(dims, 5L), c(2, 2, 2))
  )
  # The former dense path needed a 18,000 x 18,000 eigendecomposition here
  # (several GB, many minutes).
  elapsed <- system.time(
    result <- cluster4d(vec, mask, n_clusters = 50L, method = "commute")
  )[["elapsed"]]
  expect_identical(result$n_clusters, 50L)
  expect_identical(result$parameters$eigen_solver, "sparse")
  expect_lt(elapsed, 120)
})


test_that("auto commute works above 3000 voxels with required eigencore", {
  fixture <- make_commute_contract_data(c(16L, 14L, 14L))
  result <- cluster4d(
    fixture$vec, fixture$mask, n_clusters = 3L, method = "commute",
    ncomp = 5L
  )
  expect_identical(result$parameters$eigen_solver, "sparse")
  expect_identical(result$n_clusters, 3L)
  expect_identical(dim(result$embedding), c(3136L, 5L))
  expect_true(all(is.finite(result$embedding)))
  expect_true(all(result$eigenvalues > 0))
  expect_true(validate_cluster4d(result, fixture$vec, fixture$mask)$valid)
})


test_that("eigencore sparse solves are reproducible and preserve RNG state", {
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
  run <- function() commute_cluster(
    fixture$vec, fixture$mask, K = 1L, ncomp = 3L,
    connectivity = 3L, eigen_solver = "sparse"
  )
  if (had_seed) rm(".Random.seed", envir = .GlobalEnv)
  first <- run()
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  expect_true(first$metadata$rng$used)
  expect_identical(first$metadata$rng$purpose, "eigencore")
  expect_identical(first$metadata$rng$local_seed, 1L)

  set.seed(481)
  before <- .Random.seed
  second <- run()
  expect_identical(.Random.seed, before)
  expect_identical(second$embedding, first$embedding)
})


test_that("sparse commute independently checks original-coordinate eigenpairs", {
  graph <- Matrix::bandSparse(
    4L, k = c(-1L, 1L), diagonals = list(rep(1, 3L), rep(1, 3L))
  )
  laplacian <- diag(as.numeric(Matrix::rowSums(graph))) - as.matrix(graph)
  oracle <- eigen(laplacian, symmetric = TRUE)
  selected <- order(oracle$values)[1:2]
  candidate <- list(
    values = oracle$values[selected], vectors = oracle$vectors[, selected],
    nconv = 0L, certificate = list(passed = FALSE)
  )
  testthat::local_mocked_bindings(
    eig_partial = function(...) candidate,
    .package = "eigencore"
  )
  run <- function() neurocluster:::.commute_spectral_embedding(
    graph, ncomp = 1L, eigen_solver = "sparse"
  )
  accepted <- run()
  expect_false(accepted$solver_diagnostics$certificate$passed)
  expect_lt(max(accepted$solver_diagnostics$relative_residuals), 1e-10)

  # A package certificate alone cannot establish our numerical contract.
  candidate$certificate$passed <- TRUE
  candidate$nconv <- 2L
  candidate$values <- candidate$values + 1
  expect_error(run(), "sparse eigensolver did not converge")
  # Repeated exact eigenvectors have zero residual but are not orthogonal.
  candidate$values <- rep(oracle$values[selected[1L]], 2L)
  candidate$vectors <- oracle$vectors[, rep(selected[1L], 2L)]
  expect_error(run(), "sparse eigensolver did not converge")
  candidate$vectors <- matrix(0, 3L, 2L)
  expect_error(run(), "sparse eigensolver did not converge")
  candidate$vectors <- matrix(NA_real_, 4L, 2L)
  expect_error(run(), "sparse eigensolver did not converge")
  candidate$vectors <- oracle$vectors[, selected]
  candidate$values <- c(0, Inf)
  expect_error(run(), "sparse eigensolver did not converge")
})
