#!/usr/bin/env Rscript

# Linear-memory comparison-metric benchmark and correctness oracle.
# Run from the package root with:
#   Rscript inst/benchmarks/bench_metrics_contract.R

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("This benchmark requires devtools.")
}
devtools::load_all(".", quiet = TRUE)

dense_pairwise_dice <- function(a, b) {
  same_a <- outer(a, a, "==")
  same_b <- outer(b, b, "==")
  keep <- upper.tri(same_a)
  true_positive <- sum(same_a[keep] & same_b[keep])
  false_positive <- sum(same_a[keep] & !same_b[keep])
  false_negative <- sum(!same_a[keep] & same_b[keep])
  denominator <- 2 * true_positive + false_positive + false_negative
  if (denominator == 0) 1 else 2 * true_positive / denominator
}

measure_seconds <- function(operation, repetitions = 3L) {
  stopifnot(is.function(operation))
  min(replicate(repetitions, system.time(operation())[["elapsed"]]))
}

measure_allocations <- function(operation) {
  stopifnot(is.function(operation))
  profile <- tempfile("neurocluster-metric-rprofmem-")
  on.exit(unlink(profile), add = TRUE)
  utils::Rprofmem(profile)
  on.exit(utils::Rprofmem(NULL), add = TRUE)
  operation()
  utils::Rprofmem(NULL)
  lines <- readLines(profile, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  bytes <- bytes[is.finite(bytes)]
  c(total_bytes = sum(bytes), peak_allocation_bytes = max(c(0, bytes)))
}

set.seed(4301)
oracle_n <- 1200L
oracle_a <- sample(seq_len(30L), oracle_n, replace = TRUE)
oracle_b <- sample(seq_len(35L), oracle_n, replace = TRUE)
linear <- neurocluster:::.partition_metrics(oracle_a, oracle_b)$pairwise_dice
dense <- dense_pairwise_dice(oracle_a, oracle_b)
stopifnot(isTRUE(all.equal(linear, dense, tolerance = 1e-14)))

cat("Correctness oracle: linear pairwise Dice equals dense pair oracle\n")
cat(sprintf("  N=%d, pairwise Dice=%.12f\n\n", oracle_n, linear))

scales <- c(1000L, 10000L, 50000L, 100000L)
rows <- lapply(scales, function(n) {
  a <- seq_len(n)
  b <- rev(a)
  elapsed <- measure_seconds(
    function() neurocluster:::.partition_metrics(a, b),
    repetitions = if (n <= 10000L) 5L else 3L
  )
  allocations <- measure_allocations(function() {
    neurocluster:::.partition_metrics(a, b)
  })
  stats <- neurocluster:::.partition_statistics(a, b)
  stopifnot(stats$n_nonzero_cells == n)
  stopifnot(as.numeric(object.size(stats)) < 20 * n)
  data.frame(
    n = n,
    elapsed_seconds = elapsed,
    retained_bytes = as.numeric(object.size(stats)),
    total_allocated_bytes = unname(allocations[["total_bytes"]]),
    peak_allocation_bytes = unname(allocations[["peak_allocation_bytes"]]),
    hypothetical_dense_two_logical_matrices_bytes = 8 * as.double(n)^2
  )
})
scaling <- do.call(rbind, rows)
print(scaling, row.names = FALSE)

large <- scaling[scaling$n == 100000L, ]
stopifnot(large$retained_bytes < 20 * large$n)
stopifnot(large$peak_allocation_bytes < 20 * large$n)

cat("\nDense-baseline timing at the bounded oracle size:\n")
cat(sprintf(
  "  linear %.6fs; dense %.6fs\n",
  measure_seconds(function() {
    neurocluster:::.partition_metrics(oracle_a, oracle_b)
  }),
  measure_seconds(function() dense_pairwise_dice(oracle_a, oracle_b))
))
cat("\nPASS: observed-cell computation remained linear-memory at N=100,000.\n")
