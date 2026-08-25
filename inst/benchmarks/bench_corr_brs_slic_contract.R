#!/usr/bin/env Rscript

# Corr/BRS-SLIC speed-error benchmark against exact Pearson score evaluation.
# Run from the package root after loading a compiled development build:
#   Rscript inst/benchmarks/bench_corr_brs_slic_contract.R

suppressPackageStartupMessages(devtools::load_all(quiet = TRUE))

set.seed(8252026)
dims <- c(16L, 8L, 8L)
n <- prod(dims)
timepoints <- 48L
k <- 16L
mask0 <- 0:(n - 1L)
coords <- cbind(
  mask0 %% dims[1L],
  (mask0 %/% dims[1L]) %% dims[2L],
  mask0 %/% (dims[1L] * dims[2L])
) + 1

# Spatial blocks receive correlated signals plus voxel noise.
block <- 1L + (coords[, 1L] > 8L) + 2L * (coords[, 2L] > 4L) +
  4L * (coords[, 3L] > 4L)
signals <- matrix(rnorm(8L * timepoints), 8L, timepoints)
features <- signals[block, , drop = FALSE] + matrix(rnorm(n * timepoints, sd = 0.7), n, timepoints)

corr_args <- list(
  feat = features, mask_lin_idx = mask0, dims = dims, K = k,
  d = 16L, sketch_repeats = 1L, alpha = 0.25, max_iter = 5L,
  seed = 2026L, assign_stride = 1L, quantize_assign = FALSE,
  embed_basis = "hash", whiten_embed = FALSE,
  refine_exact_iters = 0L, refine_boundary_only = TRUE,
  refine_stride = 1L, refine_alpha = 0.25,
  connectivity = 6L, min_size = 1L, n_threads = 1L, verbose = FALSE
)
brs_args <- list(
  feat = features, mask_lin_idx = mask0, dims = dims, K = k,
  d = 16L, sketch_repeats = 1L, alpha = 0.25, coarse_iter = 5L,
  boundary_passes = 1L, global_passes = 0L,
  refine_spatial = 0.25, refine_l2 = 0, refine_stride = 1L,
  seed = 2026L, connectivity = 6L, min_size = 1L,
  n_threads = 1L, verbose = FALSE
)

replace_args <- function(base, overrides) {
  base[names(overrides)] <- overrides
  base
}

specs <- list(
  corr_float = list(fun = corrslic_core, args = corr_args),
  corr_quantized = list(
    fun = corrslic_core,
    args = replace_args(corr_args, list(quantize_assign = TRUE))
  ),
  corr_strided = list(
    fun = corrslic_core,
    args = replace_args(corr_args, list(assign_stride = 3L))
  ),
  corr_exact_refined = list(
    fun = corrslic_core,
    args = replace_args(corr_args, list(
      refine_exact_iters = 4L, refine_boundary_only = FALSE
    ))
  ),
  brs_boundary = list(fun = brs_slic_core, args = brs_args),
  brs_global = list(
    fun = brs_slic_core,
    args = replace_args(brs_args, list(boundary_passes = 2L, global_passes = 2L))
  ),
  brs_strided = list(
    fun = brs_slic_core,
    args = replace_args(brs_args, list(
      boundary_passes = 2L, global_passes = 1L, refine_stride = 3L
    ))
  )
)

elapsed <- numeric(length(specs))
results <- vector("list", length(specs))
names(elapsed) <- names(results) <- names(specs)
for (name in names(specs)) {
  timings <- numeric(3L)
  for (i in seq_along(timings)) {
    timing <- system.time({
      result <- do.call(specs[[name]]$fun, specs[[name]]$args)
    })
    timings[i] <- unname(timing[["elapsed"]])
  }
  elapsed[name] <- stats::median(timings)
  results[[name]] <- result
}

reference <- results$corr_exact_refined
reference_labels <- as.integer(reference$labels)
spatial_scale <- max(1, (n / k)^(1 / 3))

exact_objective <- function(result, spatial_weight = 0.25) {
  scores <- corrslic_exact_scores_cpp(
    features,
    result$original_centers,
    coords,
    result$centers_xyz,
    spatial_weight = spatial_weight,
    spatial_scale = spatial_scale,
    stride = 1L,
    l2_weight = 0
  )
  mean(scores[cbind(seq_len(nrow(scores)), as.integer(result$labels))])
}

objective <- vapply(results, exact_objective, numeric(1L))
reference_objective <- unname(objective[["corr_exact_refined"]])
ari <- vapply(results, function(result) {
  clustering_accuracy(as.integer(result$labels), reference_labels)$ari
}, numeric(1L))

report <- data.frame(
  mode = names(specs),
  median_elapsed_sec = unname(elapsed),
  speedup_vs_exact_refinement = unname(elapsed[["corr_exact_refined"]] / elapsed),
  ari_vs_exact_refinement = unname(ari),
  exact_pearson_objective = unname(objective),
  objective_regret = unname(objective - reference_objective),
  actual_k = vapply(results, function(x) as.integer(x$actual_k), integer(1L)),
  row.names = NULL
)

print(report, digits = 5, row.names = FALSE)
