#!/usr/bin/env Rscript

# Reproducible runtime and visible-allocation measurements for the fail-closed
# cluster4d validator. Run from the package root with:
# Rscript inst/benchmarks/bench_validate_cluster4d.R

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE, recompile = FALSE))

make_case <- function(name, dims, n_time, k) {
  set.seed(20260825 + prod(dims))
  sp <- neuroim2::NeuroSpace(
    dims, spacing = c(1.5, 1.5, 2), origin = c(-20, -30, 10)
  )
  frames <- lapply(seq_len(n_time), function(i) {
    neuroim2::NeuroVol(array(stats::rnorm(prod(dims)), dims), sp)
  })
  vec <- do.call(neuroim2::concat, frames)
  mask <- neuroim2::NeuroVol(array(1, dims), sp)
  labels <- rep(seq_len(k), length.out = prod(dims))
  raw <- structure(
    list(
      cluster = as.integer(labels),
      parameters = list(n_clusters_requested = k),
      metadata = list()
    ),
    class = c("cluster4d_result", "cluster_result", "list")
  )
  result <- neurocluster:::finalize_cluster4d_result(
    raw, vec, mask, "benchmark", raw$parameters
  )
  list(name = name, vec = vec, mask = mask, result = result)
}

profile_once <- function(case) {
  profile <- tempfile("validate-cluster4d-", fileext = ".mem")
  on.exit(unlink(profile), add = TRUE)
  gc()
  Rprofmem(profile)
  elapsed <- system.time({
    validation <- validate_cluster4d(case$result, case$vec, case$mask)
  })[["elapsed"]]
  Rprofmem(NULL)
  if (!validation$valid) stop(paste(validation$errors, collapse = "; "))
  lines <- readLines(profile, warn = FALSE)
  bytes <- suppressWarnings(as.numeric(sub(" .*", "", lines)))
  c(elapsed_s = elapsed, visible_alloc_bytes = sum(bytes, na.rm = TRUE))
}

benchmark_case <- function(case, repetitions) {
  measurements <- replicate(repetitions, profile_once(case))
  data.frame(
    case = case$name,
    n_voxels = length(case$result$labels),
    n_features = ncol(case$result$centers),
    actual_k = case$result$actual_k,
    repetitions = repetitions,
    median_elapsed_s = unname(stats::median(measurements["elapsed_s", ])),
    p95_elapsed_s = unname(stats::quantile(
      measurements["elapsed_s", ], 0.95, names = FALSE
    )),
    median_visible_alloc_bytes = unname(stats::median(
      measurements["visible_alloc_bytes", ]
    )),
    stringsAsFactors = FALSE
  )
}

roi <- make_case("representative_roi", c(20L, 20L, 10L), 20L, 100L)
whole_mask <- make_case("representative_whole_mask", c(50L, 50L, 20L), 20L, 400L)

results <- rbind(
  benchmark_case(roi, repetitions = 10L),
  benchmark_case(whole_mask, repetitions = 3L)
)
write.table(results, row.names = FALSE, sep = "\t", quote = FALSE)
