#!/usr/bin/env Rscript

# Dense commute-time eigensolver benchmark with a retained dense eigen/subspace
# and effective-resistance oracle. Run from the package root:
#   Rscript inst/benchmarks/bench_commute_contract.R

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

repetitions <- as.integer(Sys.getenv("NEUROCLUSTER_BENCH_REPS", "1"))
if (is.na(repetitions) || repetitions < 1L) {
  stop("NEUROCLUSTER_BENCH_REPS must be a positive integer")
}

make_graph <- function(n) {
  left <- seq_len(n - 1L)
  right <- left + 1L
  weights <- 1 + (left %% 7L) / 10
  Matrix::sparseMatrix(
    i = c(left, right), j = c(right, left),
    x = c(weights, weights), dims = c(n, n)
  )
}

oracle <- function(graph, ncomp) {
  degree <- Matrix::rowSums(graph)
  laplacian <- diag(degree) - as.matrix(graph)
  decomposition <- eigen(laplacian, symmetric = TRUE)
  order_idx <- order(decomposition$values)
  values <- decomposition$values[order_idx]
  vectors <- decomposition$vectors[, order_idx, drop = FALSE]
  positive <- which(values > 1e-10)[seq_len(ncomp)]
  list(values = values[positive], vectors = vectors[, positive, drop = FALSE])
}

grid <- expand.grid(
  n = c(40L, 100L, 200L),
  ncomp = c(2L, 8L, 20L),
  KEEP.OUT.ATTRS = FALSE
)
rows <- vector("list", nrow(grid))
for (case_index in seq_len(nrow(grid))) {
  specification <- grid[case_index, ]
  graph <- make_graph(specification$n)
  expression <- function() neurocluster:::.commute_spectral_embedding(
    graph, specification$ncomp, eigen_tol = 1e-10
  )
  invisible(expression())
  if (requireNamespace("bench", quietly = TRUE)) {
    measured <- bench::mark(
      result = expression(), iterations = repetitions,
      check = FALSE, memory = TRUE
    )
    seconds <- as.numeric(measured$median)
    allocation <- as.numeric(measured$mem_alloc)
  } else {
    seconds <- median(replicate(
      repetitions, system.time(expression())[["elapsed"]]
    ))
    allocation <- NA_real_
  }
  result <- expression()
  reference <- oracle(graph, specification$ncomp)
  eigen_error <- max(abs(result$eigenvalues - reference$values))
  subspace_error <- max(abs(
    result$eigenvectors %*% t(result$eigenvectors) -
      reference$vectors %*% t(reference$vectors)
  ))
  if (eigen_error > 1e-10 || subspace_error > 1e-9) {
    stop("spectral oracle exceeded tolerance")
  }
  rows[[case_index]] <- data.frame(
    n = specification$n,
    ncomp = specification$ncomp,
    runtime_seconds = seconds,
    allocation_bytes = allocation,
    dense_laplacian_bytes = 8 * specification$n^2,
    eigen_error = eigen_error,
    subspace_error = subspace_error
  )
}

# Full-eigenspace squared distances must equal graph volume times effective
# resistance on the tractable oracle graph.
small_graph <- make_graph(40L)
full <- neurocluster:::.commute_spectral_embedding(small_graph, 39L, 1e-12)
laplacian <- diag(Matrix::rowSums(small_graph)) - as.matrix(small_graph)
decomposition <- eigen(laplacian, symmetric = TRUE)
positive <- which(decomposition$values > 1e-12)
pseudoinverse <- decomposition$vectors[, positive, drop = FALSE] %*%
  diag(1 / decomposition$values[positive], nrow = length(positive)) %*%
  t(decomposition$vectors[, positive, drop = FALSE])
resistance <- outer(diag(pseudoinverse), diag(pseudoinverse), "+") -
  2 * pseudoinverse
distance_error <- max(abs(
  unname(as.matrix(stats::dist(full$embedding))^2) -
    unname(sum(Matrix::rowSums(small_graph)) * resistance)
))
if (distance_error > 1e-8) stop("effective-resistance oracle exceeded tolerance")

results <- do.call(rbind, rows)
print(results, row.names = FALSE)
cat(sprintf("Full-space effective-resistance error (N=40): %.3e\n", distance_error))
