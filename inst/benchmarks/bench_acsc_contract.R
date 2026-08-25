#!/usr/bin/env Rscript

# ACSC metric and graph-builder benchmark with independent pairwise oracles.
# Run from the package root:
#   Rscript inst/benchmarks/bench_acsc_contract.R

suppressPackageStartupMessages(devtools::load_all(quiet = TRUE))

set.seed(8252026)
grid <- as.matrix(expand.grid(x = 0:7, y = 0:3, z = 0:2))
nb <- nrow(grid)
timepoints <- 60L
latent <- matrix(rnorm(12L * timepoints), 12L, timepoints)
group <- 1L + (grid[, 1L] >= 4L) + 2L * (grid[, 2L] >= 2L) +
  4L * (grid[, 3L] >= 2L)
summaries <- latent[group, , drop = FALSE] +
  matrix(rnorm(nb * timepoints, sd = 0.65), nb, timepoints)
outlier_rows <- seq.int(5L, nb, by = 17L)
summaries[cbind(outlier_rows, rep(11L, length(outlier_rows)))] <- 25
block_summary <- list(
  summaries = summaries,
  spatial = grid + 1,
  grid = grid
)
ann_k <- 8L

direct_bicor <- function(x, y) {
  transform_one <- function(z) {
    residual <- z - stats::median(z)
    scale <- stats::median(abs(residual))
    if (scale <= sqrt(.Machine$double.eps)) return(rep(0, length(z)))
    u <- residual / (9 * scale)
    residual * ifelse(abs(u) < 1, (1 - u^2)^2, 0)
  }
  xw <- transform_one(x)
  yw <- transform_one(y)
  denominator <- sqrt(sum(xw^2) * sum(yw^2))
  if (denominator <= sqrt(.Machine$double.eps)) return(0)
  sum(xw * yw) / denominator
}

direct_similarity <- function(metric) {
  switch(
    metric,
    pearson = stats::cor(t(summaries), method = "pearson"),
    spearman = stats::cor(t(summaries), method = "spearman"),
    robust = outer(seq_len(nb), seq_len(nb), Vectorize(function(i, j) {
      direct_bicor(summaries[i, ], summaries[j, ])
    }))
  )
}

oracle_pairs <- function(similarity) {
  diag(similarity) <- -Inf
  nearest <- t(vapply(seq_len(nb), function(i) {
    order(-similarity[i, ], seq_len(nb))[seq_len(ann_k)]
  }, integer(ann_k)))
  pairs <- cbind(rep(seq_len(nb), each = ann_k), as.vector(t(nearest)))
  pairs <- cbind(pmin(pairs[, 1L], pairs[, 2L]), pmax(pairs[, 1L], pairs[, 2L]))
  unique(paste(pairs[, 1L], pairs[, 2L], sep = ":"))
}

graph_metrics <- function(graph, oracle, oracle_edge_keys) {
  edges <- igraph::as_data_frame(graph, what = "edges")
  from <- as.integer(edges$from)
  to <- as.integer(edges$to)
  edge_keys <- paste(pmin(from, to), pmax(from, to), sep = ":")
  list(
    edge_recall = mean(oracle_edge_keys %in% edge_keys),
    max_correlation_error = max(abs(edges$correlation - oracle[cbind(from, to)]))
  )
}

rows <- list()
for (metric in c("pearson", "spearman", "robust")) {
  oracle_elapsed <- system.time(oracle <- direct_similarity(metric))[["elapsed"]]
  oracle_edges <- oracle_pairs(oracle)
  for (builder in c("exact", "dct12")) {
    args <- list(
      block_summary = block_summary,
      ann_k = ann_k,
      alpha = 1,
      spatial_weighting = "binary",
      block_size = 1,
      correlation_metric = metric,
      knn_proj_dim = if (builder == "exact") 0L else 12L,
      knn_proj_method = if (builder == "exact") "none" else "dct"
    )
    timings <- numeric(3L)
    for (iteration in seq_along(timings)) {
      timings[iteration] <- system.time(
        graph <- do.call(build_acsc_graph, args)
      )[["elapsed"]]
    }
    accuracy <- graph_metrics(graph, oracle, oracle_edges)
    rows[[length(rows) + 1L]] <- data.frame(
      metric = metric,
      builder = builder,
      blocks = nb,
      timepoints = timepoints,
      median_elapsed_sec = stats::median(timings),
      direct_pairwise_sec = oracle_elapsed,
      edge_recall_vs_direct = accuracy$edge_recall,
      max_edge_correlation_error = accuracy$max_correlation_error
    )
  }
}

spatial_timings <- numeric(3L)
for (iteration in seq_along(spatial_timings)) {
  spatial_timings[iteration] <- system.time({
    spatial_graph <- build_acsc_graph(
      block_summary, ann_k, alpha = 0, spatial_weighting = "binary",
      block_size = 1, correlation_metric = "robust"
    )
  })[["elapsed"]]
}
cat("spatial-only builder median seconds:", stats::median(spatial_timings),
    "edges:", igraph::ecount(spatial_graph), "\n")
print(do.call(rbind, rows), digits = 6, row.names = FALSE)
