#!/usr/bin/env Rscript

# Fixed-fixture SLIC benchmark. Native diagnostics time connectivity and final
# center recomputation separately; every case retains an independent rowsum
# center oracle and a flood-fill connectivity check.

suppressPackageStartupMessages(
  devtools::load_all(".", quiet = TRUE, recompile = FALSE)
)

make_case <- function(dims, n_features = 8L) {
  n <- prod(dims)
  index <- seq_len(n)
  coords <- arrayInd(index, dims)
  set.seed(20260825 + n)
  features <- matrix(rnorm(n * n_features), nrow = n)
  list(
    dims = as.integer(dims),
    features = features,
    coords = coords,
    mask_idx = as.integer(index - 1L),
    K = as.integer(max(2L, n %/% 16L))
  )
}

connected_labels <- function(labels, dims) {
  space <- neuroim2::NeuroSpace(dims)
  mask <- neuroim2::NeuroVol(array(1, dims), space)
  graph <- neurocluster:::.exact_k_graph(mask, 6L)
  normalized <- neurocluster:::.exact_k_connected_labels(
    labels, graph$graph, graph$edges
  )
  length(unique(normalized)) == length(unique(labels))
}

benchmark_case <- function(dims, enforce, repetitions = 3L) {
  fixture <- make_case(dims)
  output <- NULL
  elapsed <- numeric(repetitions)
  connectivity_ms <- center_ms <- numeric(repetitions)
  for (iteration in seq_len(repetitions)) {
    elapsed[[iteration]] <- system.time({
      output <- neurocluster:::slic4d_core(
        fixture$features, fixture$coords, fixture$mask_idx, fixture$dims,
        c(1, 1, 1), fixture$K, compactness = 4, max_iter = 2L,
        step_mm = 0, n_threads = 1L, seed_method = "mask_grid",
        enforce_connectivity = enforce, min_size = 0L, connectivity = 6L,
        strict_connectivity = enforce, preserve_k = FALSE,
        topup_iters = 0L, grad_masked = numeric(),
        seed_relocate_radius = 0L, verbose = FALSE
      )
    })[["elapsed"]]
    connectivity_ms[[iteration]] <- output$connectivity_elapsed_ms
    center_ms[[iteration]] <- output$center_recompute_elapsed_ms
  }

  counts <- tabulate(output$labels, nbins = output$actual_k)
  feature_oracle <- rowsum(
    fixture$features, output$labels, reorder = TRUE
  ) / counts
  coord_oracle <- rowsum(
    fixture$coords, output$labels, reorder = TRUE
  ) / counts

  data.frame(
    dims = paste(dims, collapse = "x"),
    n_voxels = nrow(fixture$features),
    n_features = ncol(fixture$features),
    requested_k = fixture$K,
    actual_k = output$actual_k,
    connectivity = if (enforce) "strict" else "off",
    repetitions = repetitions,
    median_elapsed_s = stats::median(elapsed),
    p95_elapsed_s = stats::quantile(elapsed, 0.95, names = FALSE),
    median_connectivity_ms = stats::median(connectivity_ms),
    median_center_recompute_ms = stats::median(center_ms),
    centers_match_oracle = isTRUE(all.equal(
      unname(output$center_feats), unname(feature_oracle), tolerance = 1e-12
    )) && isTRUE(all.equal(
      unname(output$center_coords), unname(coord_oracle), tolerance = 1e-12
    )),
    strict_labels_connected = !enforce || connected_labels(
      output$labels, fixture$dims
    ),
    stringsAsFactors = FALSE
  )
}

cases <- list(c(16L, 16L, 2L), c(24L, 24L, 2L), c(32L, 32L, 2L))
results <- do.call(rbind, lapply(cases, function(dims) {
  rbind(benchmark_case(dims, FALSE), benchmark_case(dims, TRUE))
}))
write.table(results, row.names = FALSE, sep = "\t", quote = FALSE)

parallel_fixture <- make_case(c(16L, 16L, 8L))
run_parallel_case <- function(n_threads, repetitions = 3L) {
  output <- NULL
  elapsed <- numeric(repetitions)
  for (iteration in seq_len(repetitions)) {
    elapsed[[iteration]] <- system.time({
      output <- neurocluster:::slic4d_core(
        parallel_fixture$features, parallel_fixture$coords,
        parallel_fixture$mask_idx, parallel_fixture$dims, c(1, 1, 1),
        parallel_fixture$K, compactness = 4, max_iter = 2L,
        step_mm = 0, n_threads = n_threads, seed_method = "mask_grid",
        enforce_connectivity = FALSE, min_size = 0L, connectivity = 6L,
        strict_connectivity = FALSE, preserve_k = FALSE,
        topup_iters = 0L, grad_masked = numeric(),
        seed_relocate_radius = 0L, verbose = FALSE
      )
    })[["elapsed"]]
  }
  list(output = output, median = stats::median(elapsed))
}

sequential <- run_parallel_case(1L)
automatic <- run_parallel_case(0L)
stopifnot(
  !isTRUE(sequential$output$assignment_parallel_used),
  isTRUE(automatic$output$assignment_parallel_used),
  identical(sequential$output$labels, automatic$output$labels),
  identical(sequential$output$center_feats, automatic$output$center_feats),
  identical(sequential$output$center_coords, automatic$output$center_coords)
)
write.table(
  data.frame(
    path = c("serial", "auto_parallel"),
    n_voxels = nrow(parallel_fixture$features),
    median_elapsed_s = c(sequential$median, automatic$median),
    parallel_used = c(
      sequential$output$assignment_parallel_used,
      automatic$output$assignment_parallel_used
    ),
    exact_parity = TRUE
  ),
  row.names = FALSE, sep = "\t", quote = FALSE
)
