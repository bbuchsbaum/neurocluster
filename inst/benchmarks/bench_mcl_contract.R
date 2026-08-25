#!/usr/bin/env Rscript

# Correctness-gated MCL benchmark. It varies graph size/density, inflation,
# pruning, and exact-K mode while reporting flow invariants and dense-reference
# error. Run from the package root:
#   Rscript inst/benchmarks/bench_mcl_contract.R

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

repetitions <- as.integer(Sys.getenv("NEUROCLUSTER_BENCH_REPS", "1"))
if (is.na(repetitions) || repetitions < 1L) {
  stop("NEUROCLUSTER_BENCH_REPS must be a positive integer")
}

dense_normalize <- function(x) {
  sums <- colSums(x)
  for (column in which(sums == 0)) x[column, column] <- 1
  sweep(x, 2L, colSums(x), "/")
}

dense_step <- function(flow, inflation, expansion, prune_k, threshold) {
  expanded <- flow
  for (power in seq_len(expansion - 1L)) expanded <- expanded %*% flow
  normalized <- dense_normalize(expanded^inflation)
  pruned <- matrix(0, nrow(flow), ncol(flow))
  for (column in seq_len(ncol(flow))) {
    values <- normalized[, column]
    keep <- which(values >= threshold)
    if (!length(keep)) keep <- order(-values, seq_along(values))[1L]
    if (length(keep) > prune_k) {
      keep <- keep[order(-values[keep], keep)[seq_len(prune_k)]]
    }
    pruned[keep, column] <- values[keep]
  }
  dense_normalize(pruned)
}

make_fixture <- function(side, connectivity) {
  dims <- c(side, side, 1L)
  space <- neuroim2::NeuroSpace(dims)
  mask <- neuroim2::NeuroVol(array(1, dims), space)
  coords <- arrayInd(seq_len(prod(dims)), dims)
  features <- cbind(
    sin(coords[, 1] / 2),
    cos(coords[, 2] / 3),
    sin((coords[, 1] + coords[, 2]) / 4),
    coords[, 1] / side
  )
  graph <- neurocluster:::.mcl_build_weighted_graph(
    features, coords, mask, seq_len(nrow(coords)),
    connectivity = connectivity, spatial_weight = 0.3,
    feature_metric = "euclidean", feature_sigma = 1,
    spatial_sigma = 1, min_edge_weight = 0
  )
  list(
    n = nrow(coords),
    edges = length(graph@x) %/% 2L,
    graph = graph,
    features = features,
    mask = mask,
    connectivity = connectivity,
    target_k = max(2L, as.integer(round(nrow(coords) / 8)))
  )
}

run_case <- function(fixture, inflation, threshold, exact_k) {
  fit <- neurocluster:::.mcl_sparse(
    fixture$graph, inflation = inflation, expansion = 2L,
    max_iter = 4L, tol = 1e-8, prune_k = 32L,
    prune_threshold = threshold, loop_weight = 1
  )
  labels <- fit$labels
  if (exact_k) {
    labels <- neurocluster:::force_exact_k(
      labels, fixture$features, fixture$target_k,
      fixture$mask, fixture$connectivity
    )
  }
  list(fit = fit, labels = labels)
}

dense_reference <- function(fixture, inflation, threshold, iterations) {
  flow <- dense_normalize(as.matrix(fixture$graph) + diag(fixture$n))
  for (iteration in seq_len(iterations)) {
    flow <- dense_step(flow, inflation, 2L, 32L, threshold)
  }
  flow
}

grid <- expand.grid(
  side = c(8L, 12L),
  connectivity = c(6L, 26L),
  inflation = c(1.4, 2.2),
  prune_threshold = c(0, 1e-4),
  exact_k = c(FALSE, TRUE),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

rows <- vector("list", nrow(grid))
for (case_index in seq_len(nrow(grid))) {
  specification <- grid[case_index, ]
  fixture <- make_fixture(specification$side, specification$connectivity)
  expression <- function() run_case(
    fixture, specification$inflation,
    specification$prune_threshold, specification$exact_k
  )

  if (requireNamespace("bench", quietly = TRUE)) {
    measured <- bench::mark(
      result = expression(), iterations = repetitions,
      check = FALSE, memory = TRUE
    )
    seconds <- as.numeric(measured$median)
    allocation <- as.numeric(measured$mem_alloc)
  } else {
    times <- replicate(repetitions, system.time(expression())[["elapsed"]])
    seconds <- median(times)
    allocation <- NA_real_
  }
  result <- expression()
  fit <- result$fit
  reference_error <- NA_real_
  if (fixture$n <= 64L) {
    reference <- dense_reference(
      fixture, specification$inflation,
      specification$prune_threshold, fit$iterations
    )
    reference_error <- max(abs(reference - as.matrix(fit$flow)))
    if (reference_error > 1e-11) stop("dense-reference error exceeded tolerance")
  }
  flow_error <- max(abs(Matrix::colSums(fit$flow) - 1))
  if (flow_error > 1e-10 || any(fit$flow@x < 0) ||
      any(!is.finite(fit$flow@x))) {
    stop("flow invariant failed")
  }
  if (specification$exact_k && max(result$labels) != fixture$target_k) {
    stop("exact-K benchmark case missed its target")
  }

  rows[[case_index]] <- data.frame(
    n = fixture$n,
    edges = fixture$edges,
    connectivity = specification$connectivity,
    inflation = specification$inflation,
    prune_threshold = specification$prune_threshold,
    exact_k = specification$exact_k,
    achieved_k = max(result$labels),
    runtime_seconds = seconds,
    iterations = fit$iterations,
    nnz_final = length(fit$flow@x),
    sparsity = length(fit$flow@x) / fixture$n^2,
    allocation_bytes = allocation,
    flow_error = flow_error,
    dense_reference_error = reference_error
  )
}

results <- do.call(rbind, rows)
print(results, row.names = FALSE)
