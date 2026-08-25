#!/usr/bin/env Rscript

# Fixed-fixture benchmark for the native-linkage remediation. The pinned
# pre-remediation candidate aborts on this exact SLIC path (exit 134), so it
# cannot supply a valid timing distribution. This script records the repaired
# candidate's wall time and R-visible allocation footprint without pretending
# that the crashed baseline is a performance sample.

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
script <- if (length(file_arg) == 1L) {
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
} else {
  normalizePath("inst/benchmarks/bench_native_safety.R", mustWork = TRUE)
}
root <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)

option_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

iterations <- as.integer(option_value("iterations", "20"))
if (is.na(iterations) || iterations < 3L) stop("--iterations must be at least 3")
output_path <- option_value("output")

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("pkgload is required to benchmark the source tree")
}
pkgload::load_all(root, recompile = FALSE, quiet = TRUE)

dims <- c(8L, 8L, 4L)
n <- prod(dims)
coords <- cbind(
  rep(seq_len(dims[[1L]]) - 1L, dims[[2L]] * dims[[3L]]),
  rep(rep(seq_len(dims[[2L]]) - 1L, each = dims[[1L]]), dims[[3L]]),
  rep(seq_len(dims[[3L]]) - 1L, each = dims[[1L]] * dims[[2L]])
)
features <- outer(seq_len(n), seq_len(8L), function(i, j) {
  sin(i / 11 + j / 7) + cos(i * j / 37) + i / (10 * n)
})
features <- sweep(features, 1L, sqrt(rowSums(features^2)), "/")

run_slic <- function() {
  neurocluster:::slic4d_core(
    features, coords, seq_len(n) - 1L, dims, c(1, 1, 1), 16L,
    compactness = 10, max_iter = 3L, n_threads = 1L,
    seed_method = "mask_grid", enforce_connectivity = FALSE,
    strict_connectivity = FALSE, preserve_k = FALSE
  )
}

warmup <- run_slic()
stopifnot(length(warmup$labels) == n, all(is.finite(warmup$center_feats)))

elapsed <- numeric(iterations)
for (i in seq_len(iterations)) {
  invisible(gc(FALSE))
  elapsed[[i]] <- system.time(value <- run_slic(), gcFirst = FALSE)[["elapsed"]]
  stopifnot(length(value$labels) == n, all(is.finite(value$center_feats)))
}

profile <- tempfile("neurocluster-native-alloc-", fileext = ".out")
on.exit(unlink(profile), add = TRUE)
Rprofmem(profile)
allocation_value <- run_slic()
Rprofmem(NULL)
stopifnot(length(allocation_value$labels) == n)

allocation_lines <- readLines(profile, warn = FALSE)
allocation_bytes <- suppressWarnings(as.numeric(sub(" .*", "", allocation_lines)))
allocation_bytes <- sum(allocation_bytes[is.finite(allocation_bytes)])

git_sha <- tryCatch(
  system2("git", c("-C", root, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_
)
if (length(git_sha) != 1L) git_sha <- NA_character_
git_dirty_status <- suppressWarnings(system2(
  "git", c("-C", root, "diff", "--quiet"), stdout = FALSE, stderr = FALSE
))
candidate_state <- if (identical(as.integer(git_dirty_status), 0L)) "clean" else "dirty"

result <- data.frame(
  benchmark = "slic4d_core_native_linkage",
  baseline_ref = "18ee880",
  baseline_status = "abort_exit_134",
  candidate_sha = git_sha,
  candidate_state = candidate_state,
  n_voxels = n,
  n_features = ncol(features),
  k = 16L,
  iterations = iterations,
  median_elapsed_s = unname(median(elapsed)),
  p95_elapsed_s = unname(quantile(elapsed, 0.95, names = FALSE)),
  r_allocation_bytes = allocation_bytes,
  baseline_elapsed_s = NA_real_,
  baseline_r_allocation_bytes = NA_real_,
  elapsed_overhead_s = NA_real_,
  allocation_overhead_bytes = NA_real_,
  comparison_note = "baseline aborted before a valid performance sample",
  allocation_scope = "R allocator only",
  stringsAsFactors = FALSE
)

print(result, row.names = FALSE)
if (!is.null(output_path)) {
  write.csv(result, output_path, row.names = FALSE)
}
