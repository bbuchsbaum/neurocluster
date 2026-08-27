#!/usr/bin/env Rscript

# Slice-MSF cost and correctness benchmark. Run from the package root after
# installing/loading neurocluster. Set SLICE_MSF_BENCH_REPS to repeat timings.

suppressPackageStartupMessages(devtools::load_all(quiet = TRUE, recompile = FALSE))

adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  choose2 <- function(x) x * (x - 1) / 2
  n <- length(a)
  total <- choose2(n)
  if (total == 0) return(1)
  same_both <- sum(choose2(tab))
  same_a <- sum(choose2(rowSums(tab)))
  same_b <- sum(choose2(colSums(tab)))
  expected <- same_a * same_b / total
  maximum <- 0.5 * (same_a + same_b)
  if (maximum == expected) return(1)
  (same_both - expected) / (maximum - expected)
}

make_fixture <- function() {
  set.seed(701)
  dims <- c(12L, 12L, 3L)
  n_time <- 40L
  tt <- 0:(n_time - 1L)
  basis <- vapply(
    seq_len(18L),
    function(k) cos(pi * (tt + 0.5) * k / n_time),
    numeric(n_time)
  )
  values <- array(0, c(dims, n_time))
  truth <- array(0L, dims)
  for (z in seq_len(dims[3L])) {
    for (y in seq_len(dims[2L])) {
      for (x in seq_len(dims[1L])) {
        region <- 1L + (x - 1L) %/% 3L + 4L * ((y - 1L) %/% 6L)
        coefficients <- rnorm(18L, sd = 0.06)
        coefficients[region] <- 1
        coefficients[region + 8L] <- 0.4
        values[x, y, z, ] <-
          as.vector(basis %*% coefficients) + rnorm(n_time, sd = 0.02)
        truth[x, y, z] <- region
      }
    }
  }
  list(
    vec = NeuroVec(values, NeuroSpace(c(dims, n_time))),
    mask = NeuroVol(array(1, dims), NeuroSpace(dims)),
    truth = as.integer(truth),
    dims = dims,
    n_time = n_time
  )
}

reps <- as.integer(Sys.getenv("SLICE_MSF_BENCH_REPS", "3"))
if (is.na(reps) || reps < 1L) stop("SLICE_MSF_BENCH_REPS must be positive")
fixture <- make_fixture()
specs <- neurocluster:::.slice_msf_run_specs(
  fixture$n_time, r = 8L, num_runs = 3L, seed = 17L,
  ensemble_fraction = 0.6
)

run_one <- function(spec) {
  slice_msf_single(
    fixture$vec, fixture$mask,
    r = 8L, k = 2 / 11, min_size = 3L, nbhd = 8L,
    stitch_z = TRUE, gamma = 0,
    dct_frequencies = spec$frequencies,
    dct_weights = spec$weights
  )
}

timed <- function(expr) {
  values <- numeric(reps)
  result <- NULL
  for (i in seq_len(reps)) {
    timing <- system.time(result <- force(expr()))
    values[i] <- unname(timing[["elapsed"]])
  }
  list(seconds_median = stats::median(values), result = result)
}

single <- timed(function() run_one(list(
  frequencies = seq_len(8L), weights = rep(1, 8L)
)))
ensemble <- timed(function() lapply(specs, run_one))
consensus <- timed(function() slice_msf_consensus(
  ensemble$result, fixture$mask,
  nbhd = 8L, k_fuse = 0.4, min_size_fuse = 16L,
  use_features = TRUE, lambda = 0.7, stitch_z = TRUE
))

run_labels <- lapply(ensemble$result, function(run) run$labels[run$labels > 0L])
consensus_labels <- consensus$result$labels[consensus$result$labels > 0L]
raw_truth_ari <- vapply(
  run_labels, adjusted_rand, numeric(1L), b = fixture$truth
)
pairwise <- combn(seq_along(run_labels), 2L, function(pair) {
  adjusted_rand(run_labels[[pair[1L]]], run_labels[[pair[2L]]])
})
consensus_truth_ari <- adjusted_rand(consensus_labels, fixture$truth)
replay <- lapply(specs, run_one)

metrics <- data.frame(
  metric = c(
    "single_seconds_median",
    "three_run_ensemble_seconds_median",
    "consensus_only_seconds_median",
    "distinct_raw_partitions",
    "raw_pairwise_ari_mean",
    "raw_truth_ari_mean",
    "consensus_truth_ari",
    "consensus_truth_ari_gain",
    "consensus_cluster_count",
    "seeded_replay_exact"
  ),
  value = c(
    single$seconds_median,
    ensemble$seconds_median,
    consensus$seconds_median,
    length(unique(vapply(run_labels, paste, character(1L), collapse = ","))),
    mean(pairwise),
    mean(raw_truth_ari),
    consensus_truth_ari,
    consensus_truth_ari - mean(raw_truth_ari),
    length(unique(consensus_labels)),
    identical(
      lapply(ensemble$result, `[[`, "labels"),
      lapply(replay, `[[`, "labels")
    )
  )
)

print(metrics, row.names = FALSE)
if (metrics$value[metrics$metric == "distinct_raw_partitions"] < 3) {
  stop("ensemble diversity contract failed")
}
if (!isTRUE(as.logical(metrics$value[metrics$metric == "seeded_replay_exact"]))) {
  stop("seeded replay contract failed")
}
