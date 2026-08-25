#!/usr/bin/env Rscript

# Correctness-gated throughput, allocation, and thread-scaling benchmark for
# the supervoxels assignment kernel. Run from the package root:
#   Rscript inst/benchmarks/bench_supervoxels_parallel_contract.R

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

repetitions <- as.integer(Sys.getenv("NEUROCLUSTER_BENCH_REPS", "3"))
if (is.na(repetitions) || repetitions < 1L) {
  stop("NEUROCLUSTER_BENCH_REPS must be a positive integer")
}

make_fixture <- function(dims, K = 64L, features = 12L) {
  coords <- arrayInd(seq_len(prod(dims)), dims)
  n <- nrow(coords)
  labels <- rep(seq_len(K), length.out = n)
  labels <- labels[order(order(coords[, 1] + dims[1] * (
    coords[, 2] - 1L + dims[2] * (coords[, 3] - 1L)
  )))]
  data <- vapply(seq_len(features), function(j) {
    sin(seq_len(n) / (7 + j)) + cos(rowSums(coords) / (3 + j))
  }, numeric(n))
  data <- t(scale(data, center = TRUE, scale = TRUE))
  state <- neurocluster:::.supervoxel_center_state(
    labels, data, coords, K, alpha = 0.55, sigma1 = 2,
    sigma2 = 3, parallel = FALSE
  )

  # Two deterministic local neighbors are sufficient because the binned path
  # gathers its primary candidates from the centroid grid.
  nn_index <- cbind(c(seq_len(n - 1L), n) - 1L,
                    c(1L, seq_len(n - 1L)) - 1L)
  nn_dist <- matrix(1, nrow = n, ncol = 2L)
  list(
    n = n,
    K = K,
    coords = t(coords),
    data = data,
    labels = state$labels - 1L,
    counts = state$counts,
    feature_centers = state$centers,
    coord_centers = state$coord_centers,
    nn_index = nn_index,
    nn_dist = nn_dist
  )
}

run_assignment <- function(fixture, parallel, grain_size = 256L) {
  fn <- if (parallel) {
    fused_assignment_parallel_binned
  } else {
    fused_assignment_binned
  }
  args <- list(
    fixture$nn_index, fixture$nn_dist, fixture$labels,
    fixture$coords, fixture$feature_centers, fixture$coord_centers,
    fixture$counts, fixture$data,
    1, 2, 3, 0.55
  )
  if (parallel) {
    args <- c(args, list(
      grain_size = grain_size, window_factor = 3, bin_expand = 1L
    ))
  } else {
    args <- c(args, list(window_factor = 3, bin_expand = 1L))
  }
  do.call(fn, args)
}

elapsed <- function(fn, reps) {
  timings <- replicate(reps, system.time(fn())[["elapsed"]])
  unname(median(timings))
}

fixtures <- list(
  small = make_fixture(c(24L, 24L, 8L)),
  large = make_fixture(c(32L, 32L, 16L))
)

for (name in names(fixtures)) {
  fixture <- fixtures[[name]]
  sequential_labels <- run_assignment(fixture, parallel = FALSE)
  parallel_labels <- run_assignment(fixture, parallel = TRUE)
  stopifnot(identical(sequential_labels, parallel_labels))
}

timings <- do.call(rbind, lapply(names(fixtures), function(name) {
  fixture <- fixtures[[name]]
  data.frame(
    size = name,
    voxels = fixture$n,
    mode = c("sequential", "parallel"),
    seconds = c(
      elapsed(function() run_assignment(fixture, FALSE), repetitions),
      elapsed(function() run_assignment(fixture, TRUE), repetitions)
    )
  )
}))

allocations <- NULL
if (requireNamespace("bench", quietly = TRUE)) {
  fixture <- fixtures$large
  measured <- bench::mark(
    sequential = run_assignment(fixture, FALSE),
    parallel = run_assignment(fixture, TRUE),
    iterations = repetitions,
    check = TRUE,
    memory = TRUE
  )
  allocations <- data.frame(
    mode = as.character(measured$expression),
    median_seconds = as.numeric(measured$median),
    memory_bytes = as.numeric(measured$mem_alloc)
  )
}

available_threads <- max(1L, RcppParallel::defaultNumThreads())
thread_counts <- unique(pmin(c(1L, 2L, 4L), available_threads))
thread_scaling <- do.call(rbind, lapply(thread_counts, function(threads) {
  RcppParallel::setThreadOptions(numThreads = threads)
  data.frame(
    threads = threads,
    seconds = elapsed(
      function() run_assignment(fixtures$large, TRUE), repetitions
    )
  )
}))
RcppParallel::setThreadOptions(numThreads = available_threads)

small_n <- fixtures$small$n
large_n <- fixtures$large$n
for (mode in c("sequential", "parallel")) {
  small_time <- timings$seconds[timings$size == "small" & timings$mode == mode]
  large_time <- timings$seconds[timings$size == "large" & timings$mode == mode]
  if (!is.finite(small_time) || !is.finite(large_time) ||
      small_time <= 0 || large_time <= 0) {
    stop("non-finite benchmark timing")
  }
}

cat("Supervoxels assignment correctness: exact sequential/parallel parity\n")
print(timings, row.names = FALSE)
cat(sprintf(
  "Voxel scaling: %.2fx (%d to %d voxels)\n",
  large_n / small_n, small_n, large_n
))
if (!is.null(allocations)) {
  cat("Allocation profile on the large fixture:\n")
  print(allocations, row.names = FALSE)
} else {
  cat("Allocation profile skipped: install the optional bench package\n")
}
cat("Parallel thread scaling on the large fixture:\n")
print(thread_scaling, row.names = FALSE)
