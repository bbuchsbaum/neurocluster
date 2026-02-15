# Benchmark corr_slic optimizations with speed/quality gates.
#
# Usage:
#   Rscript inst/benchmarks/bench_corr_slic_ablations.R
#
# Optional env vars:
#   NEUROCLUSTER_BENCH_THREADS=1   # thread count (default: 1 for stable timing)
#   NEUROCLUSTER_BENCH_REPS=3      # reps per dataset/variant (default: 3)

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
  } else {
    library(neurocluster)
  }
})

set.seed(123)

bench_threads <- suppressWarnings(as.integer(Sys.getenv("NEUROCLUSTER_BENCH_THREADS", "1")))
if (is.na(bench_threads) || bench_threads < 1L) bench_threads <- 1L
bench_reps <- suppressWarnings(as.integer(Sys.getenv("NEUROCLUSTER_BENCH_REPS", "3")))
if (is.na(bench_reps) || bench_reps < 1L) bench_reps <- 3L

Sys.setenv(OMP_NUM_THREADS = as.character(bench_threads))
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = as.character(bench_threads))

datasets <- list(
  blobs_mid = function() {
    generate_synthetic_volume(
      scenario = "gaussian_blobs",
      dims = c(24, 24, 8),
      n_clusters = 6,
      n_time = 60,
      pattern_type = "orthogonal",
      noise_sd = 0.12,
      seed = 11
    )
  },
  z_layers_soft = function() {
    generate_synthetic_volume(
      scenario = "z_layers",
      dims = c(24, 24, 8),
      n_clusters = 5,
      n_time = 60,
      pattern_type = "harmonic",
      noise_sd = 0.18,
      seed = 12
    )
  },
  checkerboard = function() {
    generate_synthetic_volume(
      scenario = "checkerboard",
      dims = c(20, 20, 6),
      n_clusters = 4,
      n_time = 50,
      pattern_type = "orthogonal",
      noise_sd = 0.10,
      seed = 13
    )
  }
)

variants <- list(
  baseline = list(assign_stride = 1L),
  stride2 = list(assign_stride = 2L),
  stride3 = list(assign_stride = 3L),
  quant8 = list(assign_stride = 1L, quantize_assign = TRUE),
  stride2_quant8 = list(assign_stride = 2L, quantize_assign = TRUE)
)

run_once <- function(syn, variant_args) {
  gc()
  t0 <- proc.time()[["elapsed"]]
  res <- cluster4d(
    syn$vec, syn$mask,
    n_clusters = syn$n_clusters,
    method = "corr_slic",
    spatial_weight = 0.08,
    embedding_dim = 64L,
    sketch_repeats = 1L,
    max_iterations = 6L,
    connectivity = 6L,
    parallel = (bench_threads != 1L),
    verbose = FALSE,
    seed = 7L,
    assign_stride = as.integer(variant_args$assign_stride),
    quantize_assign = isTRUE(variant_args$quantize_assign)
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  acc <- clustering_accuracy(res$cluster, syn$truth)
  data.frame(
    elapsed_sec = elapsed,
    ari = acc$ari,
    nmi = acc$nmi,
    accuracy = acc$accuracy,
    n_clusters_found = res$n_clusters
  )
}

rows <- list()
row_i <- 1L
for (ds_name in names(datasets)) {
  syn <- datasets[[ds_name]]()
  for (variant_name in names(variants)) {
    for (rep_i in seq_len(bench_reps)) {
      out <- run_once(syn, variants[[variant_name]])
      out$dataset <- ds_name
      out$variant <- variant_name
      out$rep <- rep_i
      out$threads <- bench_threads
      rows[[row_i]] <- out
      row_i <- row_i + 1L
      cat(sprintf(
        "%s | %s | rep=%d | t=%.3fs | ARI=%.3f | NMI=%.3f\n",
        ds_name, variant_name, rep_i, out$elapsed_sec, out$ari, out$nmi
      ))
    }
  }
}

raw <- do.call(rbind, rows)
raw <- raw[, c("dataset", "variant", "rep", "threads", "elapsed_sec", "ari", "nmi", "accuracy", "n_clusters_found")]

agg <- aggregate(
  cbind(elapsed_sec, ari, nmi, accuracy, n_clusters_found) ~ dataset + variant,
  data = raw,
  FUN = mean
)
baseline <- agg[agg$variant == "baseline", c("dataset", "elapsed_sec", "ari", "nmi", "accuracy")]
colnames(baseline) <- c("dataset", "elapsed_base", "ari_base", "nmi_base", "acc_base")
summary <- merge(agg, baseline, by = "dataset", all.x = TRUE, sort = FALSE)
summary$speedup_vs_base <- summary$elapsed_base / summary$elapsed_sec
summary$ari_delta_vs_base <- summary$ari - summary$ari_base
summary$nmi_delta_vs_base <- summary$nmi - summary$nmi_base
summary$acc_delta_vs_base <- summary$accuracy - summary$acc_base

out_dir <- file.path("inst", "benchmarks")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
raw_path <- file.path(out_dir, "corr_slic_ablations_raw.csv")
sum_path <- file.path(out_dir, "corr_slic_ablations_summary.csv")
write.csv(raw, raw_path, row.names = FALSE)
write.csv(summary, sum_path, row.names = FALSE)

cat("\n=== corr_slic ablation summary ===\n")
print(summary[order(summary$dataset, summary$variant), ], row.names = FALSE)
cat(sprintf("\nSaved raw results to %s\n", raw_path))
cat(sprintf("Saved summary results to %s\n", sum_path))
