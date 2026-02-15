# Tune corr_slic hyperparameters by dataset family.
#
# Usage:
#   Rscript inst/benchmarks/tune_corr_slic.R
#
# Optional env vars:
#   NEUROCLUSTER_BENCH_THREADS=4
#   NEUROCLUSTER_TUNE_REPS=1
#   NEUROCLUSTER_TUNE_FAST=1

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
  } else {
    library(neurocluster)
  }
})

set.seed(123)

bench_threads <- suppressWarnings(as.integer(Sys.getenv("NEUROCLUSTER_BENCH_THREADS", "4")))
if (is.na(bench_threads) || bench_threads < 1L) bench_threads <- 4L
bench_reps <- suppressWarnings(as.integer(Sys.getenv("NEUROCLUSTER_TUNE_REPS", "1")))
if (is.na(bench_reps) || bench_reps < 1L) bench_reps <- 1L
tune_fast <- identical(Sys.getenv("NEUROCLUSTER_TUNE_FAST", "1"), "1")
score_lambda <- suppressWarnings(as.numeric(Sys.getenv("NEUROCLUSTER_TUNE_LAMBDA", "0.15")))
if (!is.finite(score_lambda) || score_lambda < 0) score_lambda <- 0.15

Sys.setenv(OMP_NUM_THREADS = as.character(bench_threads))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = as.character(bench_threads))

dataset_specs <- list(
  blobs_small = function() {
    generate_synthetic_volume(
      scenario = "gaussian_blobs",
      dims = c(20, 20, 6),
      n_clusters = 4,
      n_time = 40,
      pattern_type = "orthogonal",
      noise_sd = 0.10,
      seed = 2
    )
  },
  blobs_noisy = function() {
    generate_synthetic_volume(
      scenario = "gaussian_blobs",
      dims = c(20, 20, 6),
      n_clusters = 4,
      n_time = 50,
      pattern_type = "orthogonal",
      noise_sd = 0.50,
      seed = 7
    )
  },
  blobs_overlap_soft = function() {
    generate_synthetic_volume(
      scenario = "gaussian_blobs",
      dims = c(24, 24, 10),
      n_clusters = 6,
      n_time = 60,
      spread = c(2.0, 2.0, 2.0),
      pattern_type = "random",
      noise_sd = 0.25,
      seed = 11
    )
  },
  z_layers_soft = function() {
    generate_synthetic_volume(
      scenario = "z_layers",
      dims = c(28, 28, 10),
      n_clusters = 5,
      n_time = 60,
      pattern_type = "harmonic",
      noise_sd = 0.20,
      seed = 12
    )
  },
  touching_blobs = function() {
    generate_synthetic_volume(
      scenario = "gaussian_blobs",
      dims = c(32, 32, 8),
      n_clusters = 6,
      n_time = 50,
      pattern_type = "orthogonal",
      noise_sd = 0.10,
      spread = c(4.0, 4.0, 2.5),
      seed = 3
    )
  }
)

family_map <- list(
  clean_blobs = c("blobs_small"),
  noisy_blobs = c("blobs_noisy", "blobs_overlap_soft"),
  layered = c("z_layers_soft"),
  topology = c("touching_blobs")
)

if (isTRUE(tune_fast)) {
  param_grid <- expand.grid(
    spatial_weight = c(0.05, 0.08),
    embedding_dim = c(24L, 32L),
    embedding_basis = c("hash", "dct"),
    embedding_whiten = c(FALSE, TRUE),
    sketch_repeats = c(1L, 2L),
    quantize_assign = c(TRUE),
    assign_stride = c(2L, 3L),
    max_iterations = c(4L, 6L),
    refine_exact_iters = c(0L, 1L, 2L),
    refine_boundary_only = c(TRUE, FALSE),
    refine_stride = c(1L, 2L),
    stringsAsFactors = FALSE
  )
} else {
  param_grid <- expand.grid(
    spatial_weight = c(0.03, 0.05, 0.08, 0.10),
    embedding_dim = c(NA_integer_, 16L, 24L, 32L, 48L, 64L),
    embedding_basis = c("hash", "dct"),
    embedding_whiten = c(FALSE, TRUE),
    sketch_repeats = c(1L, 2L),
    quantize_assign = c(FALSE, TRUE),
    assign_stride = c(1L, 2L, 3L),
    max_iterations = c(3L, 4L, 5L, 6L),
    refine_exact_iters = c(0L, 1L, 2L, 3L),
    refine_boundary_only = c(TRUE, FALSE),
    refine_stride = c(1L, 2L, 4L),
    stringsAsFactors = FALSE
  )
}

# Prune obviously-irrelevant combinations to keep tuning runtime reasonable.
param_grid <- subset(param_grid, !(embedding_basis == "dct" & sketch_repeats != 1L))
param_grid <- subset(param_grid, !(refine_exact_iters == 0L & (refine_boundary_only == FALSE | refine_stride != 1L)))

run_one <- function(syn, p) {
  gc()
  emb_dim <- if (is.na(p$embedding_dim)) "auto" else as.integer(p$embedding_dim)
  adaptive <- is.na(p$embedding_dim)
  t0 <- proc.time()[["elapsed"]]
  res <- cluster4d(
    syn$vec, syn$mask,
    n_clusters = syn$n_clusters,
    method = "corr_slic",
    spatial_weight = as.numeric(p$spatial_weight),
    embedding_dim = emb_dim,
    adaptive_embedding = adaptive,
    embedding_basis = p$embedding_basis,
    embedding_whiten = isTRUE(p$embedding_whiten),
    sketch_repeats = as.integer(p$sketch_repeats),
    max_iterations = as.integer(p$max_iterations),
    connectivity = 6L,
    quantize_assign = isTRUE(p$quantize_assign),
    assign_stride = as.integer(p$assign_stride),
    refine_exact_iters = as.integer(p$refine_exact_iters),
    refine_boundary_only = isTRUE(p$refine_boundary_only),
    refine_stride = as.integer(p$refine_stride),
    parallel = (bench_threads != 1L),
    verbose = FALSE,
    seed = 7L
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  acc <- clustering_accuracy(res$cluster, syn$truth)
  list(
    elapsed_sec = elapsed,
    ari = acc$ari,
    nmi = acc$nmi,
    n_clusters_found = res$n_clusters,
    emb_dim_used = res$parameters$embedding_dim
  )
}

rows <- list()
row_i <- 1L

for (family_name in names(family_map)) {
  for (ds_name in family_map[[family_name]]) {
    syn <- dataset_specs[[ds_name]]()
    for (pid in seq_len(nrow(param_grid))) {
      p <- param_grid[pid, , drop = FALSE]
      for (rep_i in seq_len(bench_reps)) {
        out <- tryCatch(run_one(syn, p), error = function(e) e)
        if (inherits(out, "error")) {
          rows[[row_i]] <- cbind(
            data.frame(
              family = family_name,
              dataset = ds_name,
              param_id = pid,
              rep = rep_i,
              elapsed_sec = NA_real_,
              ari = NA_real_,
              nmi = NA_real_,
              n_clusters_found = NA_integer_,
              emb_dim_used = NA_integer_,
              error = conditionMessage(out),
              stringsAsFactors = FALSE
            ),
            p
          )
          row_i <- row_i + 1L
          next
        }

        rows[[row_i]] <- cbind(
          data.frame(
            family = family_name,
            dataset = ds_name,
            param_id = pid,
            rep = rep_i,
            elapsed_sec = out$elapsed_sec,
            ari = out$ari,
            nmi = out$nmi,
            n_clusters_found = out$n_clusters_found,
            emb_dim_used = out$emb_dim_used,
            error = NA_character_,
            stringsAsFactors = FALSE
          ),
          p
        )
        row_i <- row_i + 1L
      }
    }
  }
}

raw <- do.call(rbind, rows)
ok <- subset(raw, is.na(error))

# Score trades quality against runtime.
agg <- aggregate(
  cbind(elapsed_sec, ari, nmi) ~ family + param_id +
    spatial_weight + embedding_dim + embedding_basis + embedding_whiten +
    sketch_repeats + quantize_assign + assign_stride + max_iterations +
    refine_exact_iters + refine_boundary_only + refine_stride,
  data = ok,
  FUN = mean
)
agg$score <- agg$ari - score_lambda * log1p(agg$elapsed_sec)

best <- do.call(rbind, lapply(split(agg, agg$family), function(df) {
  df[order(-df$score, -df$ari, df$elapsed_sec), ][1, , drop = FALSE]
}))

out_dir <- file.path("inst", "benchmarks")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
raw_path <- file.path(out_dir, "corr_slic_tuning_raw.csv")
agg_path <- file.path(out_dir, "corr_slic_tuning_summary.csv")
best_path <- file.path(out_dir, "corr_slic_tuning_best_by_family.csv")
write.csv(raw, raw_path, row.names = FALSE)
write.csv(agg, agg_path, row.names = FALSE)
write.csv(best, best_path, row.names = FALSE)

cat("=== corr_slic tuning: best by family ===\n")
print(best, row.names = FALSE)
cat(sprintf("\nSaved raw trials to %s\n", raw_path))
cat(sprintf("Saved summary to %s\n", agg_path))
cat(sprintf("Saved best configs to %s\n", best_path))
