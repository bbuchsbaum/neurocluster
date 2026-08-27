#!/usr/bin/env Rscript

# Proof-bearing Slice-MSF characterization. Run from the package root. The
# script records public defaults, natural versus repaired K, the signed
# temporal-smoothness diagnostic, unweighted edge separation, and fail-closed
# gamma behavior. It deliberately uses the same fixed spherical fixture as the
# regression test.

suppressPackageStartupMessages(
  devtools::load_all(quiet = TRUE, recompile = FALSE)
)
source("tests/testthat/helper-slice-msf-oracle.R", local = TRUE)

fixture <- slice_msf_make_spherical_fixture()
truth <- fixture$truth
mask_index <- which(as.vector(as.array(fixture$mask)) > 0)

partition_row <- function(name, labels) {
  sizes <- sort(tabulate(as.integer(factor(labels))), decreasing = TRUE)
  data.frame(
    fit = name,
    k = length(sizes),
    ari = clustering_accuracy(labels, truth)$ari,
    largest = sizes[[1L]],
    smallest = sizes[[length(sizes)]],
    singleton_count = sum(sizes == 1L),
    sizes = paste(sizes, collapse = ","),
    stringsAsFactors = FALSE
  )
}

run_natural <- function() {
  slice_msf(
    fixture$vec, fixture$mask,
    target_k_global = -1L, r = 16L, compactness = 4, min_size = 2L,
    num_runs = 1L, consensus = FALSE, stitch_z = TRUE, nbhd = 8L
  )
}

run_exact <- function(...) {
  cluster4d(
    fixture$vec, fixture$mask, n_clusters = 8L, method = "slice_msf",
    spatial_weight = 0.4, min_size = 2L, r = 16L,
    num_runs = 1L, consensus = FALSE, verbose = FALSE, ...
  )
}

natural <- run_natural()
public_exact <- run_exact()
explicit_safe_exact <- run_exact(gamma = 0, z_mult = 0)
partitions <- do.call(rbind, list(
  partition_row("public_natural", natural$cluster),
  partition_row("public_exact_k", public_exact$cluster),
  partition_row("explicit_gamma_0_z_0_exact_k", explicit_safe_exact$cluster)
))

repair <- data.frame(
  natural_k = length(unique(natural$cluster)),
  requested_k = 8L,
  repair_direction = if (length(unique(natural$cluster)) > 8L) "merge" else "split",
  repair_ran = TRUE,
  pre_sizes = paste(sort(tabulate(natural$cluster), decreasing = TRUE), collapse = ","),
  post_sizes = paste(sort(tabulate(public_exact$cluster), decreasing = TRUE), collapse = ","),
  public_default_matches_explicit_safe = identical(
    as.integer(public_exact$cluster), as.integer(explicit_safe_exact$cluster)
  ),
  stringsAsFactors = FALSE
)

capture_gamma_error <- function(path, expression) {
  condition <- tryCatch(
    {
      force(expression)
      NULL
    },
    error = identity
  )
  data.frame(
    path = path,
    rejected = inherits(condition, "slice_msf_unsupported_gamma"),
    condition_class = if (is.null(condition)) "none" else class(condition)[[1L]],
    message = if (is.null(condition)) "" else conditionMessage(condition),
    stringsAsFactors = FALSE
  )
}

gamma_contract <- rbind(
  capture_gamma_error("direct", slice_msf(
    fixture$vec, fixture$mask, min_size = 2L,
    num_runs = 1L, consensus = FALSE, gamma = 1e-6
  )),
  capture_gamma_error("cluster4d", run_exact(gamma = 1e-6))
)

native <- slice_msf_single(
  fixture$vec, fixture$mask,
  r = 16L, k = 0.4, min_size = 2L,
  nbhd = 8L, stitch_z = TRUE, gamma = 0, z_mult = 0
)
smoothness <- native$temporal_smoothness[mask_index]
temporal_smoothness <- data.frame(
  statistic = c(
    "minimum", "q25", "median", "mean", "q75", "maximum",
    "fraction_non_positive", "fraction_negative"
  ),
  value = c(
    min(smoothness), stats::quantile(smoothness, 0.25),
    stats::median(smoothness), mean(smoothness),
    stats::quantile(smoothness, 0.75), max(smoothness),
    mean(smoothness <= 0), mean(smoothness < 0)
  )
)

graph <- neurocluster:::.exact_k_graph(fixture$mask, 26L)
edges <- graph$edges
sketch <- native$sketch[, mask_index, drop = FALSE]
similarity <- colSums(
  sketch[, edges[, 1L], drop = FALSE] *
    sketch[, edges[, 2L], drop = FALSE]
)
distance <- 1 - pmax(-1, pmin(1, similarity))
boundary <- truth[edges[, 1L]] != truth[edges[, 2L]]
edge_summary <- do.call(rbind, lapply(
  list(within_truth = !boundary, between_truth = boundary),
  function(select) {
    data.frame(
      n = sum(select),
      q25 = stats::quantile(distance[select], 0.25),
      median = stats::median(distance[select]),
      q75 = stats::quantile(distance[select], 0.75),
      fraction_zero = mean(distance[select] == 0)
    )
  }
))
edge_summary$edge_type <- rownames(edge_summary)
rownames(edge_summary) <- NULL
edge_summary <- edge_summary[, c("edge_type", setdiff(names(edge_summary), "edge_type"))]

report <- capture.output({
  cat("\nPARTITIONS\n")
  print(partitions, row.names = FALSE, digits = 6)
  cat("\nREPAIR\n")
  print(repair, row.names = FALSE)
  cat("\nGAMMA CONTRACT\n")
  print(gamma_contract, row.names = FALSE)
  cat("\nTEMPORAL SMOOTHNESS (DIAGNOSTIC ONLY)\n")
  print(temporal_smoothness, row.names = FALSE, digits = 6)
  cat("\nUNWEIGHTED EDGE SEPARATION\n")
  print(edge_summary, row.names = FALSE, digits = 6)
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
utils::write.csv(repair, file.path(artifact_dir, "repair.csv"), row.names = FALSE)
utils::write.csv(gamma_contract, file.path(artifact_dir, "gamma-contract.csv"), row.names = FALSE)
utils::write.csv(
  temporal_smoothness, file.path(artifact_dir, "temporal-smoothness.csv"),
  row.names = FALSE
)
utils::write.csv(edge_summary, file.path(artifact_dir, "edge-summary.csv"), row.names = FALSE)
saveRDS(
  list(
    git_sha = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    dirty = length(system2("git", c("status", "--porcelain"), stdout = TRUE)) > 0L,
    partitions = partitions,
    repair = repair,
    gamma_contract = gamma_contract,
    temporal_smoothness = temporal_smoothness,
    edge_summary = edge_summary,
    session_info = sessionInfo()
  ),
  file.path(artifact_dir, "diagnosis.rds")
)
cat("Artifacts:", normalizePath(artifact_dir), "\n")
