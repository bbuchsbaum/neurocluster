#!/usr/bin/env Rscript

# Fixed-fixture FLASH-3D throughput gate. Each 64/128-bit, barrier-off/on case
# is timed without diagnostics, then rerun with score diagnostics and checked
# against an exhaustive R implementation of the native scoring equation.

suppressPackageStartupMessages(
  devtools::load_all(".", quiet = TRUE, recompile = FALSE)
)

dims <- c(24L, 20L, 3L)
n <- prod(dims)
K <- 36L
mask_idx <- seq_len(n)
coords <- arrayInd(mask_idx, dims) - 1
set.seed(20260825)
ts <- matrix(rnorm(24L * n), nrow = 24L)
vox_scale <- c(1, 1.25, 2)
boundary <- ifelse(coords[, 1] < dims[[1]] / 2, 0, 3)

run_core <- function(bits, barrier_on, diagnostics = FALSE) {
  lambda <- c(0.35, 0.65, if (barrier_on) 1.2 else 0)
  barrier <- if (barrier_on) boundary else NULL
  flash3d_supervoxels_cpp(
    ts, mask_idx, dims, K, lambda, 1L, bits, 12L,
    vox_scale, barrier, FALSE, diagnostics
  )
}

slow_oracle_check <- function(bits, barrier_on) {
  oracle_dims <- c(12L, 10L, 2L)
  oracle_n <- prod(oracle_dims)
  oracle_k <- 12L
  oracle_idx <- seq_len(oracle_n)
  oracle_coords <- arrayInd(oracle_idx, oracle_dims) - 1
  set.seed(20260825)
  oracle_ts <- matrix(rnorm(24L * oracle_n), nrow = 24L)
  oracle_boundary <- ifelse(
    oracle_coords[, 1] < oracle_dims[[1]] / 2, 0, 3
  )
  lambda <- c(0.35, 0.65, if (barrier_on) 1.2 else 0)
  barrier_arg <- if (barrier_on) oracle_boundary else NULL
  barrier_values <- if (barrier_on) oracle_boundary else rep(0, oracle_n)
  core <- flash3d_supervoxels_cpp(
    oracle_ts, oracle_idx, oracle_dims, oracle_k, lambda, 1L, bits, 12L,
    vox_scale, barrier_arg, FALSE, TRUE
  )
  scaled <- matrix(vox_scale, nrow = oracle_n, ncol = 3L, byrow = TRUE)
  scores <- matrix(NA_real_, nrow = oracle_n, ncol = oracle_k)
  for (cluster in seq_len(oracle_k)) {
    site_bits <- matrix(
      core$assignment_site_hash_bits[cluster, ],
      nrow = oracle_n, ncol = bits, byrow = TRUE
    )
    hamming <- rowSums(abs(core$voxel_hash_bits - site_bits)) / bits
    spatial <- rowSums(
      (sweep(oracle_coords, 2L, core$assignment_site_coords[cluster, ], "-") * scaled)^2
    ) / core$S2
    barrier_cost <- abs(
      barrier_values -
        barrier_values[core$assignment_site_seed_grid[cluster]]
    )
    scores[, cluster] <- core$lambda_s_final * spatial +
      lambda[[2L]] * hamming + lambda[[3L]] * barrier_cost
  }
  oracle <- max.col(-scores, ties.method = "first")
  identical(as.integer(core$labels), as.integer(oracle))
}

visible_allocated_bytes <- function(bits, barrier_on) {
  profile <- tempfile("flash3d-rprofmem-")
  on.exit(unlink(profile), add = TRUE)
  Rprofmem(profile)
  invisible(run_core(bits, barrier_on, diagnostics = FALSE))
  Rprofmem(NULL)
  sizes <- suppressWarnings(as.numeric(sub(
    " .*", "", readLines(profile, warn = FALSE)
  )))
  sum(sizes, na.rm = TRUE)
}

benchmark_case <- function(bits, barrier_on, repetitions = 9L) {
  invisible(run_core(bits, barrier_on, diagnostics = FALSE))
  gc()
  elapsed <- numeric(repetitions)
  output <- NULL
  for (iteration in seq_len(repetitions)) {
    elapsed[[iteration]] <- system.time({
      output <- run_core(bits, barrier_on, diagnostics = FALSE)
    })[["elapsed"]]
  }
  data.frame(
    bits = bits,
    barrier = if (barrier_on) "on" else "off",
    n_voxels = n,
    n_features = nrow(ts),
    requested_k = K,
    actual_k = length(unique(output$labels)),
    reseed_count = output$reseed_count,
    oracle_fixture_n = 240L,
    slow_oracle_identical = slow_oracle_check(bits, barrier_on),
    median_elapsed_s = stats::median(elapsed),
    p95_elapsed_s = stats::quantile(elapsed, 0.95, names = FALSE),
    median_voxels_per_s = n / max(stats::median(elapsed), .Machine$double.eps),
    r_visible_alloc_bytes = visible_allocated_bytes(bits, barrier_on),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, lapply(c(64L, 128L), function(bits) {
  rbind(
    benchmark_case(bits, FALSE),
    benchmark_case(bits, TRUE)
  )
}))
write.table(results, row.names = FALSE, sep = "\t", quote = FALSE)
