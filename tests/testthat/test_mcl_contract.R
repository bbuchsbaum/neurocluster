slow_mcl_normalize <- function(x) {
  x <- as.matrix(x)
  sums <- colSums(x)
  for (column in which(sums == 0)) x[column, column] <- 1
  sweep(x, 2L, colSums(x), "/")
}

slow_mcl_prune <- function(x, max_per_col, threshold) {
  x <- as.matrix(x)
  out <- matrix(0, nrow(x), ncol(x))
  for (column in seq_len(ncol(x))) {
    values <- x[, column]
    keep <- which(values >= threshold)
    if (!length(keep)) keep <- order(-values, seq_along(values))[1L]
    if (length(keep) > max_per_col) {
      keep <- keep[order(-values[keep], keep)[seq_len(max_per_col)]]
    }
    out[keep, column] <- values[keep]
  }
  out
}

slow_mcl_step <- function(flow, inflation, expansion, prune_k, threshold) {
  expanded <- flow
  for (power in seq_len(expansion - 1L)) expanded <- expanded %*% flow
  inflated <- expanded^inflation
  normalized <- slow_mcl_normalize(inflated)
  slow_mcl_normalize(slow_mcl_prune(normalized, prune_k, threshold))
}

make_mcl_contract_data <- function() {
  dims <- c(8L, 3L, 1L)
  space <- neuroim2::NeuroSpace(dims)
  mask <- neuroim2::NeuroVol(array(1, dims), space)
  coords <- arrayInd(seq_len(prod(dims)), dims)
  volumes <- lapply(seq_len(8L), function(time) {
    values <- sin(coords[, 1] / 2 + time / 3) +
      (coords[, 1] > 4) * cos(time / 2)
    neuroim2::NeuroVol(array(values, dims), space)
  })
  list(mask = mask, vec = do.call(neuroim2::concat, volumes))
}

expect_connected_mcl_labels <- function(labels, dims) {
  coords <- arrayInd(seq_len(prod(dims)), dims)
  offsets <- rbind(
    c(-1L, 0L, 0L), c(1L, 0L, 0L),
    c(0L, -1L, 0L), c(0L, 1L, 0L),
    c(0L, 0L, -1L), c(0L, 0L, 1L)
  )
  key <- setNames(seq_len(nrow(coords)), apply(coords, 1L, paste, collapse = ":"))
  for (label in sort(unique(labels))) {
    members <- which(labels == label)
    seen <- members[1L]
    frontier <- members[1L]
    while (length(frontier)) {
      node <- frontier[1L]
      frontier <- frontier[-1L]
      candidates <- sweep(offsets, 2L, coords[node, ], "+")
      inside <- apply(candidates >= 1L &
        sweep(candidates, 2L, dims, "<="), 1L, all)
      adjacent <- unname(key[apply(
        candidates[inside, , drop = FALSE], 1L, paste, collapse = ":"
      )])
      adjacent <- adjacent[labels[adjacent] == label]
      fresh <- setdiff(adjacent, seen)
      seen <- c(seen, fresh)
      frontier <- c(frontier, fresh)
    }
    expect_setequal(seen, members)
  }
}


test_that("sparse MCL matches an independent dense reference step by step", {
  graph <- Matrix::Matrix(matrix(c(
    0, .8, .2, 0,
    .8, 0, .4, .1,
    .2, .4, 0, .7,
    0, .1, .7, 0
  ), 4L, 4L, byrow = TRUE), sparse = TRUE)
  inflation <- 1.7
  expansion <- 2L
  prune_k <- 3L
  threshold <- 0.04
  loop_weight <- 0.6

  fit <- neurocluster:::.mcl_sparse(
    graph, inflation, expansion, max_iter = 4L, tol = 1e-14,
    prune_k = prune_k, prune_threshold = threshold,
    loop_weight = loop_weight, trace = TRUE
  )
  reference <- slow_mcl_normalize(
    as.matrix(graph) + diag(loop_weight, nrow(graph))
  )
  expect_equal(as.matrix(fit$trace[[1L]]), reference, tolerance = 1e-13)
  for (iteration in seq_len(4L)) {
    reference <- slow_mcl_step(
      reference, inflation, expansion, prune_k, threshold
    )
    expect_equal(
      as.matrix(fit$trace[[iteration + 1L]]), reference,
      tolerance = 1e-12, info = paste("iteration", iteration)
    )
  }
})


test_that("flow invariants, determinism, attractors, and permutations hold", {
  graph <- Matrix::Matrix(matrix(c(
    0, .91, .17, .03, 0,
    .91, 0, .31, .07, .02,
    .17, .31, 0, .83, .11,
    .03, .07, .83, 0, .57,
    0, .02, .11, .57, 0
  ), 5L, 5L, byrow = TRUE), sparse = TRUE)
  args <- list(
    graph = graph, inflation = 1.8, expansion = 3L,
    max_iter = 5L, tol = 1e-12, prune_k = 4L,
    prune_threshold = 1e-5, loop_weight = .7, trace = TRUE
  )
  fit1 <- do.call(neurocluster:::.mcl_sparse, args)
  fit2 <- do.call(neurocluster:::.mcl_sparse, args)

  expect_identical(fit1$labels, fit2$labels)
  expect_equal(as.matrix(fit1$flow), as.matrix(fit2$flow), tolerance = 0)
  for (flow in fit1$trace) {
    expect_true(all(is.finite(flow@x)))
    expect_true(all(flow@x >= 0))
    expect_equal(Matrix::colSums(flow), rep(1, ncol(flow)), tolerance = 1e-12)
  }
  expect_identical(fit1$attractors[fit1$attractors], fit1$attractors)

  chain <- Matrix::sparseMatrix(
    i = c(2L, 3L, 3L), j = 1:3, x = 1, dims = c(3L, 3L)
  )
  expect_identical(
    neurocluster:::.mcl_attractors_from_flow(chain), rep.int(3L, 3L)
  )

  permutation <- c(3L, 5L, 1L, 4L, 2L)
  permuted_args <- args
  permuted_args$graph <- graph[permutation, permutation]
  permuted <- do.call(neurocluster:::.mcl_sparse, permuted_args)
  restored_labels <- integer(5L)
  restored_labels[permutation] <- permuted$labels
  expect_identical(
    outer(restored_labels, restored_labels, `==`),
    outer(fit1$labels, fit1$labels, `==`)
  )
  restored_flow <- matrix(0, 5L, 5L)
  restored_flow[permutation, permutation] <- as.matrix(permuted$flow)
  expect_equal(restored_flow, as.matrix(fit1$flow), tolerance = 1e-12)
})


test_that("MCL parameters fail closed and parallel is absent", {
  graph <- Matrix::Diagonal(3L)
  bad <- list(
    list(inflation = 1), list(inflation = Inf),
    list(expansion = 1L), list(expansion = 2.5),
    list(max_iter = 0L), list(loop_weight = -1),
    list(prune_k = 0L), list(prune_threshold = -1e-3),
    list(prune_threshold = 1), list(tol = 0)
  )
  for (entry in bad) {
    expect_error(do.call(neurocluster:::.mcl_sparse, c(list(graph), entry)))
  }
  expect_false("parallel" %in% names(formals(cluster4d_mcl)))
  expect_false(neurocluster:::cluster4d_method_contract("mcl")$parallel)
})


test_that("natural target search operationally selects the closest candidate", {
  data <- make_mcl_contract_data()
  low <- cluster4d_mcl(
    data$vec, data$mask, n_clusters = 1L, max_iterations = 5L,
    inflation = 1.8, exact_k = FALSE
  )
  high <- cluster4d_mcl(
    data$vec, data$mask, n_clusters = 24L, max_iterations = 5L,
    inflation = 1.8, exact_k = FALSE
  )

  for (result in list(low, high)) {
    search <- result$metadata$target_search
    target <- result$parameters$n_clusters_requested
    expect_equal(
      abs(result$actual_k - target),
      min(abs(search$achieved_k - target))
    )
    expect_identical(result$metadata$achieved_k, result$actual_k)
    expect_identical(
      result$metadata$target_policy,
      "closest_natural_k_over_deterministic_inflation_grid"
    )
  }
  expect_false(identical(low$metadata$selected_inflation,
                         high$metadata$selected_inflation))
  expect_false(identical(low$cluster, high$cluster))
})


test_that("exact-K covers one, N, under, over, and topology", {
  data <- make_mcl_contract_data()
  targets <- c(1L, 3L, 8L, 24L)
  results <- lapply(targets, function(target) {
    suppressWarnings(cluster4d_mcl(
      data$vec, data$mask, n_clusters = target,
      max_iterations = 5L, inflation = 1.8,
      connectivity = 6L, exact_k = TRUE
    ))
  })
  natural <- vapply(results, function(x) x$metadata$natural_k, integer(1))
  expect_true(any(natural > targets))
  expect_true(any(natural < targets))
  for (index in seq_along(results)) {
    result <- results[[index]]
    expect_identical(result$actual_k, targets[index])
    expect_identical(result$metadata$achieved_k, targets[index])
    expect_true(result$metadata$exact_k_applied)
    expect_connected_mcl_labels(result$labels, c(8L, 3L, 1L))
  }
  expect_identical(results[[1L]]$actual_k, 1L)
  expect_true(all(tabulate(results[[4L]]$labels, nbins = 24L) == 1L))
})


test_that("exact-K respects disconnected mask components", {
  dims <- c(3L, 1L, 1L)
  space <- neuroim2::NeuroSpace(dims)
  mask_values <- array(c(1, 0, 1), dims)
  mask <- neuroim2::NeuroVol(mask_values, space)
  volume <- neuroim2::NeuroVec(
    array(c(2, 0, 5), dim = c(dims, 1L)),
    neuroim2::add_dim(space, 1L)
  )

  expect_error(
    cluster4d_mcl(volume, mask, n_clusters = 1L, exact_k = TRUE),
    class = "cluster4d_exact_k_infeasible"
  )
  result <- cluster4d_mcl(
    volume, mask, n_clusters = 2L, exact_k = TRUE
  )
  expect_identical(result$actual_k, 2L)
  expect_identical(result$labels, 1:2)
})
