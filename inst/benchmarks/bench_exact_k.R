#!/usr/bin/env Rscript

# Absolute runtime and peak-memory characterization of the corrected sparse
# exact-K split engine. The two fixtures pin the production-scale N/K pairs
# from the Slice-MSF remediation epic. The parent process uses /usr/bin/time -l
# so every worker reports peak resident memory in an isolated process.

args <- commandArgs(trailingOnly = TRUE)

make_case <- function(case) {
  specification <- switch(
    case,
    n12544_k150 = list(dims = c(28L, 28L, 16L), target = 150L,
                       block = c(7L, 7L, 8L)),
    n32000_k300 = list(dims = c(40L, 40L, 20L), target = 300L,
                       block = c(8L, 8L, 10L)),
    stop("unknown exact-K benchmark case: ", case)
  )
  dims <- specification$dims
  sp <- neuroim2::NeuroSpace(dims)
  mask <- neuroim2::NeuroVol(array(1, dims), sp)
  coords <- arrayInd(seq_len(prod(dims)), dims)
  features <- cbind(
    sin(coords[, 1] * 0.31) + cos(coords[, 2] * 0.17),
    cos(coords[, 2] * 0.23) + coords[, 3] * 0.1,
    sin(rowSums(coords) * 0.11),
    coords[, 1] / dims[1],
    coords[, 2] / dims[2],
    coords[, 3] / dims[3]
  )
  labels <- interaction(
    (coords[, 1L] - 1L) %/% specification$block[[1L]],
    (coords[, 2L] - 1L) %/% specification$block[[2L]],
    (coords[, 3L] - 1L) %/% specification$block[[3L]],
    drop = TRUE
  )
  list(
    mask = mask,
    features = features,
    labels = as.integer(labels),
    target = specification$target,
    min_cluster_size = 20L
  )
}

topology_ok <- function(labels, mask) {
  graph <- neurocluster:::.exact_k_graph(mask, 6L)
  connected <- neurocluster:::.exact_k_connected_labels(
    labels, graph$graph, graph$edges
  )
  length(unique(connected)) == length(unique(labels))
}

if (length(args) && identical(args[[1]], "--worker")) {
  suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE, recompile = FALSE))
  case_name <- args[[2]]
  repetitions <- as.integer(args[[3]])
  fixture <- make_case(case_name)
  graph <- neurocluster:::.exact_k_graph(fixture$mask, 6L)
  elapsed <- numeric(repetitions)
  result <- NULL
  for (i in seq_len(repetitions)) {
    elapsed[[i]] <- system.time({
      result <- neurocluster:::force_exact_k(
        fixture$labels, fixture$features, fixture$target,
        fixture$mask, 6L, graph_info = graph,
        min_cluster_size = fixture$min_cluster_size
      )
    })[["elapsed"]]
  }
  output <- data.frame(
    case = case_name,
    method = "corrected_feature_tree",
    n_voxels = length(fixture$labels),
    n_edges = nrow(graph$edges),
    initial_k = length(unique(fixture$labels)),
    natural_k = length(unique(fixture$labels)),
    target_k = fixture$target,
    repaired_k = length(unique(result)),
    min_cluster_size = fixture$min_cluster_size,
    repetitions = repetitions,
    median_elapsed_s = stats::median(elapsed),
    p95_elapsed_s = stats::quantile(elapsed, 0.95, names = FALSE),
    exact_k = length(unique(result)) == fixture$target,
    topology_ok = topology_ok(result, fixture$mask),
    minimum_size_ok = min(tabulate(result)) >= fixture$min_cluster_size,
    smallest_cluster = min(tabulate(result)),
    largest_cluster = max(tabulate(result)),
    cluster_sizes = paste(tabulate(result), collapse = ","),
    stringsAsFactors = FALSE
  )
  write.table(output, row.names = FALSE, sep = "\t", quote = FALSE)
  quit(save = "no", status = 0L)
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- normalizePath(sub("^--file=", "", script_arg[[1]]))
rscript <- file.path(R.home("bin"), "Rscript")
time_bin <- "/usr/bin/time"
if (!file.exists(time_bin)) stop("/usr/bin/time is required for peak RSS measurement")

run_worker <- function(case, repetitions = 3L) {
  stdout <- tempfile("exact-k-bench-", fileext = ".tsv")
  stderr <- tempfile("exact-k-time-", fileext = ".txt")
  on.exit(unlink(c(stdout, stderr)), add = TRUE)
  status <- system2(
    time_bin,
    c("-l", rscript, script, "--worker", case, repetitions),
    stdout = stdout,
    stderr = stderr
  )
  if (status != 0L) stop(paste(readLines(stderr, warn = FALSE), collapse = "\n"))
  result <- utils::read.delim(stdout, check.names = FALSE)
  timing <- readLines(stderr, warn = FALSE)
  rss_line <- grep("maximum resident set size", timing, value = TRUE)
  if (!length(rss_line)) stop("Could not parse maximum resident set size")
  result$peak_rss_bytes <- as.numeric(sub("^\\s*([0-9]+).*", "\\1", rss_line[[1]]))
  result
}

repetitions <- as.integer(Sys.getenv("EXACT_K_BENCH_REPS", "3"))
if (!is.finite(repetitions) || repetitions < 1L) {
  stop("EXACT_K_BENCH_REPS must be a positive integer")
}
results <- do.call(rbind, lapply(
  c("n12544_k150", "n32000_k300"),
  run_worker, repetitions = repetitions
))
write.table(results, row.names = FALSE, sep = "\t", quote = FALSE)
