#!/usr/bin/env Rscript

# Focused reproduction for the near-zero Slice-MSF ARI observed in
# tests/testthat/test_supervoxel_quality.R. This is diagnostic only: it compares
# the supported path with two ablations and records the intermediate partition,
# reliability, edge-separation, and exact-K behavior.

suppressPackageStartupMessages(
  devtools::load_all(quiet = TRUE, recompile = FALSE)
)

# Load only helper definitions preceding the first test_that() call. This keeps
# the fixture exactly aligned with the quality test without running all methods.
test_expressions <- parse("tests/testthat/test_supervoxel_quality.R")
for (expression in test_expressions) {
  if (is.call(expression) && identical(expression[[1L]], as.name("test_that"))) {
    break
  }
  eval(expression, envir = .GlobalEnv)
}

data <- make_spherical_clusters(
  dims = c(15, 15, 15), n_clusters = 8, n_time = 100,
  noise_sd = 0.3, size_variation = 0.5, seed = 42
)
mask_index <- which(as.vector(as.array(data$mask)) > 0)
truth <- data$true_labels

partition_row <- function(name, labels) {
  sizes <- sort(tabulate(as.integer(factor(labels))), decreasing = TRUE)
  metrics <- clustering_accuracy(labels, truth)
  data.frame(
    fit = name,
    k = length(sizes),
    ari = metrics$ari,
    largest = sizes[[1L]],
    smallest = sizes[[length(sizes)]],
    singleton_count = sum(sizes == 1L),
    size_cv = stats::sd(sizes) / mean(sizes),
    stringsAsFactors = FALSE
  )
}

run_natural <- function(gamma, z_mult = 0.5) {
  slice_msf(
    data$vec, data$mask,
    target_k_global = -1L,
    r = 16L, compactness = 4, min_size = 2L,
    num_runs = 1L, consensus = FALSE,
    stitch_z = TRUE, nbhd = 8L, gamma = gamma, z_mult = z_mult
  )
}

run_exact <- function(gamma, z_mult = 0.5) {
  cluster4d(
    data$vec, data$mask, n_clusters = 8L, method = "slice_msf",
    spatial_weight = 0.4, min_size = 2L, r = 16L,
    num_runs = 1L, consensus = FALSE,
    gamma = gamma, z_mult = z_mult, verbose = FALSE
  )
}

supported_natural <- run_natural(gamma = 1)
supported_exact <- run_exact(gamma = 1)
unit_reliability_natural <- run_natural(gamma = 0)
unit_reliability_exact <- run_exact(gamma = 0)
unit_reliability_no_z_exact <- run_exact(gamma = 0, z_mult = 0)

partitions <- do.call(rbind, list(
  partition_row("supported_natural_gamma_1", supported_natural$cluster),
  partition_row("supported_exact_k_gamma_1", supported_exact$cluster),
  partition_row("unit_reliability_natural_gamma_0", unit_reliability_natural$cluster),
  partition_row("unit_reliability_exact_k_gamma_0", unit_reliability_exact$cluster),
  partition_row("unit_reliability_no_z_exact_k_gamma_0", unit_reliability_no_z_exact$cluster)
))

# Any strictly positive gamma preserves exact zeros from non-positive
# split-half correlations. Probe the discontinuity at gamma = 0 explicitly.
gamma_sweep <- do.call(rbind, lapply(c(0, 1e-6, 0.1, 0.5, 1, 2), function(value) {
  natural <- if (value == 0) unit_reliability_natural else if (value == 1) {
    supported_natural
  } else run_natural(value)
  exact <- if (value == 0) unit_reliability_exact else if (value == 1) {
    supported_exact
  } else run_exact(value)
  natural_sizes <- sort(tabulate(natural$cluster), decreasing = TRUE)
  exact_sizes <- sort(tabulate(exact$cluster), decreasing = TRUE)
  data.frame(
    gamma = value,
    natural_k = length(natural_sizes),
    natural_largest = natural_sizes[[1L]],
    exact_ari = clustering_accuracy(exact$cluster, truth)$ari,
    exact_largest = exact_sizes[[1L]],
    exact_singletons = sum(exact_sizes == 1L)
  )
}))

native <- slice_msf_single(
  data$vec, data$mask,
  r = 16L, k = 0.4, min_size = 2L,
  nbhd = 8L, stitch_z = TRUE, gamma = 1, z_mult = 0.5
)
weights <- native$weights[mask_index]
reliability <- data.frame(
  statistic = c(
    "minimum", "q25", "median", "mean", "q75", "maximum",
    "fraction_exact_zero", "fraction_below_0.01", "fraction_below_0.10"
  ),
  value = c(
    min(weights), stats::quantile(weights, 0.25), stats::median(weights),
    mean(weights), stats::quantile(weights, 0.75), max(weights),
    mean(weights == 0), mean(weights < 0.01), mean(weights < 0.10)
  )
)

graph <- neurocluster:::.exact_k_graph(data$mask, 26L)
edges <- graph$edges
sketch <- native$sketch[, mask_index, drop = FALSE]
similarity <- colSums(
  sketch[, edges[, 1L], drop = FALSE] *
    sketch[, edges[, 2L], drop = FALSE]
)
unweighted_distance <- 1 - pmax(-1, pmin(1, similarity))
edge_reliability <- sqrt(weights[edges[, 1L]] * weights[edges[, 2L]])
weighted_distance <- edge_reliability * unweighted_distance
boundary <- truth[edges[, 1L]] != truth[edges[, 2L]]

edge_summary <- do.call(rbind, lapply(
  list(within_truth = !boundary, between_truth = boundary),
  function(select) {
    data.frame(
      n = sum(select),
      unweighted_q25 = stats::quantile(unweighted_distance[select], 0.25),
      unweighted_median = stats::median(unweighted_distance[select]),
      unweighted_q75 = stats::quantile(unweighted_distance[select], 0.75),
      weighted_q25 = stats::quantile(weighted_distance[select], 0.25),
      weighted_median = stats::median(weighted_distance[select]),
      weighted_q75 = stats::quantile(weighted_distance[select], 0.75),
      fraction_weighted_zero = mean(weighted_distance[select] == 0)
    )
  }
))
edge_summary$edge_type <- rownames(edge_summary)
rownames(edge_summary) <- NULL
edge_summary <- edge_summary[, c("edge_type", setdiff(names(edge_summary), "edge_type"))]

original <- neurocluster:::.cluster4d_original_data(
  data$vec, data$mask, "slice_msf_diagnosis"
)
one_component_exact <- neurocluster:::force_exact_k(
  rep(1L, length(truth)), original$features, 8L,
  mask = data$mask, connectivity = 26L, graph_info = graph
)
fallback_comparison <- data.frame(
  final_equals_one_component_fallback = identical(
    as.integer(supported_exact$cluster), as.integer(one_component_exact)
  ),
  partition_ari = neurocluster:::.partition_metrics(
    supported_exact$cluster, one_component_exact
  )$ari,
  fallback_sizes = paste(
    sort(tabulate(one_component_exact), decreasing = TRUE), collapse = ","
  )
)

report <- capture.output({
  cat("\nPARTITIONS\n")
  print(partitions, row.names = FALSE, digits = 6)
  cat("\nRELIABILITY\n")
  print(reliability, row.names = FALSE, digits = 6)
  cat("\nGAMMA SWEEP\n")
  print(gamma_sweep, row.names = FALSE, digits = 6)
  cat("\nEDGE SEPARATION\n")
  print(edge_summary, row.names = FALSE, digits = 6)
  cat("\nEXACT-K FALLBACK\n")
  print(fallback_comparison, row.names = FALSE)
  cat("\nSESSION\n")
  print(sessionInfo())
})
cat(report, sep = "\n")
cat("\n")

artifact_dir <- Sys.getenv(
  "SLICE_MSF_DIAG_DIR", "/tmp/neurocluster-slice-msf-diagnosis"
)
dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(report, file.path(artifact_dir, "diagnosis.txt"))
utils::write.csv(partitions, file.path(artifact_dir, "partitions.csv"), row.names = FALSE)
utils::write.csv(reliability, file.path(artifact_dir, "reliability.csv"), row.names = FALSE)
utils::write.csv(gamma_sweep, file.path(artifact_dir, "gamma-sweep.csv"), row.names = FALSE)
utils::write.csv(edge_summary, file.path(artifact_dir, "edge-summary.csv"), row.names = FALSE)
utils::write.csv(
  fallback_comparison, file.path(artifact_dir, "exact-k-fallback.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    git_sha = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    partitions = partitions,
    reliability = reliability,
    gamma_sweep = gamma_sweep,
    edge_summary = edge_summary,
    exact_k_fallback = fallback_comparison,
    session_info = sessionInfo()
  ),
  file.path(artifact_dir, "diagnosis.rds")
)
cat("Artifacts:", normalizePath(artifact_dir), "\n")
