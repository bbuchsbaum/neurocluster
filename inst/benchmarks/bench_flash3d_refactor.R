#!/usr/bin/env Rscript

# Reproducible FLASH-3D accuracy/fidelity/throughput gate for the 2026-08-26
# physical-scale, hash, and jump-flood refactor. The recorded baselines were
# measured from commit 94fa747 immediately before the refactor on the same host.

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE, recompile = FALSE)
  library(neuroim2)
})

threads <- suppressWarnings(as.integer(Sys.getenv("FLASH3D_BENCH_THREADS", "8")))
if (is.na(threads) || threads < 1L) threads <- 1L
RcppParallel::setThreadOptions(numThreads = threads)

score_fixture <- function(spec, seed) {
  data <- do.call(generate_synthetic_volume, c(spec, list(seed = seed)))
  fit <- cluster4d(
    data$vec, data$mask, data$n_clusters,
    method = "flash3d", verbose = FALSE
  )
  observed <- !is.na(data$truth)
  clustering_accuracy(
    fit$cluster[observed], data$truth[observed]
  )$ari
}

fixtures <- list(
  blobs_orthogonal = list(
    scenario = "gaussian_blobs", dims = c(20, 20, 6), n_clusters = 4,
    n_time = 40, pattern_type = "orthogonal", noise_sd = 0.1
  ),
  blobs_harmonic = list(
    scenario = "gaussian_blobs", dims = c(20, 20, 6), n_clusters = 4,
    n_time = 60, pattern_type = "harmonic", noise_sd = 0.15
  ),
  z_layers = list(
    scenario = "z_layers", dims = c(20, 20, 10), n_clusters = 5,
    n_time = 60, pattern_type = "harmonic", noise_sd = 0.2
  )
)
pre_refactor_ari <- c(
  blobs_orthogonal = 0.5489404930,
  blobs_harmonic = 0.4912240392,
  z_layers = 0.0196261887
)

accuracy <- do.call(rbind, lapply(names(fixtures), function(name) {
  scores <- vapply(1:5, function(seed) score_fixture(fixtures[[name]], seed),
                   numeric(1L))
  observed <- median(scores)
  data.frame(
    fixture = name,
    pre_refactor_median_ari = pre_refactor_ari[[name]],
    refactored_median_ari = observed,
    ari_delta = observed - pre_refactor_ari[[name]],
    accuracy_gate = observed >= pre_refactor_ari[[name]] + 0.03
  )
}))

# Pairwise correlation-distance fidelity of the 64-bit proposal hash.
data <- do.call(generate_synthetic_volume, c(fixtures$blobs_harmonic, list(seed = 99L)))
mask_idx <- which(as.array(data$mask) > 0)
ts <- series(data$vec, mask_idx)
core <- flash3d_supervoxels_cpp(
  ts, mask_idx, as.integer(data$dims), data$n_clusters,
  c(0.25, 0.75, 0), 1L, 64L, 16L, c(1, 1, 1), NULL, FALSE, TRUE
)
set.seed(100)
pairs <- cbind(sample.int(ncol(ts), 12000L, replace = TRUE),
               sample.int(ncol(ts), 12000L, replace = TRUE))
pairs <- pairs[pairs[, 1L] != pairs[, 2L], , drop = FALSE]
exact <- vapply(seq_len(nrow(pairs)), function(i) {
  1 - cor(ts[, pairs[i, 1L]], ts[, pairs[i, 2L]])
}, numeric(1L))
hash <- vapply(seq_len(nrow(pairs)), function(i) {
  mean(core$voxel_hash_bits[pairs[i, 1L], ] !=
         core$voxel_hash_bits[pairs[i, 2L], ])
}, numeric(1L))
hash_spearman <- cor(exact, hash, method = "spearman")
pre_refactor_hash_spearman <- 0.6209271

# Metamorphic physical-unit check through the native scoring boundary.
scaled_core <- function(scale) {
  flash3d_supervoxels_cpp(
    ts, mask_idx, as.integer(data$dims), data$n_clusters,
    c(0.25, 0.75, 0), 4L, 64L, 16L, rep(scale, 3L), NULL, FALSE, TRUE
  )
}
unit <- scaled_core(1)
tripled <- scaled_core(3)
spacing_invariant <- identical(unit$labels, tripled$labels) &&
  isTRUE(all.equal(tripled$S2 / unit$S2, 9, tolerance = 1e-12))

# Representative native throughput fixture. Runtime is reported relative to
# the same-host pre-refactor measurement but is not a cross-machine hard gate.
dims <- c(32L, 32L, 16L)
n <- prod(dims)
set.seed(20260826)
perf_ts <- matrix(rnorm(60L * n), nrow = 60L)
run_perf <- function() {
  flash3d_supervoxels_cpp(
    perf_ts, seq_len(n), dims, 100L, c(0.25, 0.75, 0),
    4L, 64L, 16L, c(1, 1, 1), NULL, FALSE, FALSE
  )
}
invisible(run_perf())
elapsed <- replicate(7L, system.time(invisible(run_perf()))[["elapsed"]])
median_elapsed <- median(elapsed)
pre_refactor_median_elapsed <- 0.115

cat("ACCURACY\n")
write.table(accuracy, row.names = FALSE, sep = "\t", quote = FALSE)
cat("FIDELITY_AND_INVARIANCE\n")
write.table(data.frame(
  pre_refactor_hash_spearman,
  refactored_hash_spearman = hash_spearman,
  hash_fidelity_delta = hash_spearman - pre_refactor_hash_spearman,
  hash_fidelity_gate = hash_spearman >= pre_refactor_hash_spearman + 0.1,
  spacing_invariant
), row.names = FALSE, sep = "\t", quote = FALSE)
cat("THROUGHPUT\n")
write.table(data.frame(
  threads,
  n_voxels = n,
  n_features = nrow(perf_ts),
  requested_k = 100L,
  rounds = 4L,
  pre_refactor_median_elapsed_s = pre_refactor_median_elapsed,
  refactored_median_elapsed_s = median_elapsed,
  speedup = pre_refactor_median_elapsed / median_elapsed
), row.names = FALSE, sep = "\t", quote = FALSE)

stopifnot(all(accuracy$accuracy_gate),
          hash_spearman >= pre_refactor_hash_spearman + 0.1,
          spacing_invariant)
