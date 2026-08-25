#!/usr/bin/env Rscript

# Fixed sparse-grid scaling benchmark for the versioned adjacency-constrained
# Ward queue. Isolated workers expose rejection and queue-growth diagnostics;
# the parent uses /usr/bin/time -l to add peak resident memory.

args <- commandArgs(trailingOnly = TRUE)

make_case <- function(side, n_features = 8L) {
  side <- as.integer(side)
  n <- side^2
  index <- matrix(seq_len(n), side, side)
  left <- c(as.vector(index[-side, ]), as.vector(index[, -side]))
  right <- c(as.vector(index[-1L, ]), as.vector(index[, -1L]))
  graph <- Matrix::sparseMatrix(
    i = c(left, right), j = c(right, left), x = 1,
    dims = c(n, n)
  )
  set.seed(20260825 + side)
  list(
    X = matrix(rnorm(n * n_features), nrow = n),
    G = graph,
    sizes = sample(1:8, n, replace = TRUE),
    target = as.integer(n %/% 4L),
    n_edges = length(left)
  )
}

if (length(args) && identical(args[[1L]], "--worker")) {
  suppressPackageStartupMessages(
    devtools::load_all(".", quiet = TRUE, recompile = FALSE)
  )
  side <- as.integer(args[[2L]])
  repetitions <- as.integer(args[[3L]])
  fixture <- make_case(side)
  elapsed <- numeric(repetitions)
  output <- NULL
  for (iteration in seq_len(repetitions)) {
    elapsed[[iteration]] <- system.time({
      output <- neurocluster:::ward_on_supervoxels_cpp(
        fixture$X, fixture$G, fixture$sizes, fixture$target
      )
    })[["elapsed"]]
  }
  total_rejections <- output$stale_version_rejections +
    output$inactive_rejections + output$adjacency_rejections +
    output$recomputed_rejections
  result <- data.frame(
    side = side,
    n_supervoxels = nrow(fixture$X),
    n_undirected_edges = fixture$n_edges,
    n_features = ncol(fixture$X),
    target_k = fixture$target,
    actual_k = output$n_clusters,
    repetitions = repetitions,
    median_elapsed_s = stats::median(elapsed),
    p95_elapsed_s = stats::quantile(elapsed, 0.95, names = FALSE),
    queue_pushes = output$queue_pushes,
    queue_pops = output$queue_pops,
    max_queue_size = output$max_queue_size,
    max_queue_per_input_edge = output$max_queue_size / fixture$n_edges,
    stale_version_rejections = output$stale_version_rejections,
    inactive_rejections = output$inactive_rejections,
    adjacency_rejections = output$adjacency_rejections,
    recomputed_rejections = output$recomputed_rejections,
    rejection_fraction_of_pops = total_rejections / output$queue_pops,
    exact_k = output$n_clusters == fixture$target,
    finite_nonnegative_costs = all(is.finite(output$merge_cost)) &&
      all(output$merge_cost >= 0),
    stringsAsFactors = FALSE
  )
  write.table(result, row.names = FALSE, sep = "\t", quote = FALSE)
  quit(save = "no", status = 0L)
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- normalizePath(sub("^--file=", "", script_arg[[1L]]))
rscript <- file.path(R.home("bin"), "Rscript")
time_bin <- "/usr/bin/time"
if (!file.exists(time_bin)) {
  stop("/usr/bin/time is required for peak RSS measurement")
}

run_worker <- function(side, repetitions = 5L) {
  stdout <- tempfile("rena-plus-ward-bench-", fileext = ".tsv")
  stderr <- tempfile("rena-plus-ward-time-", fileext = ".txt")
  on.exit(unlink(c(stdout, stderr)), add = TRUE)
  status <- system2(
    time_bin,
    c("-l", rscript, script, "--worker", side, repetitions),
    stdout = stdout,
    stderr = stderr
  )
  if (status != 0L) {
    stop(paste(readLines(stderr, warn = FALSE), collapse = "\n"))
  }
  result <- utils::read.delim(stdout, check.names = FALSE)
  timing <- readLines(stderr, warn = FALSE)
  rss_line <- grep("maximum resident set size", timing, value = TRUE)
  if (!length(rss_line)) stop("Could not parse maximum resident set size")
  result$peak_rss_bytes <- as.numeric(
    sub("^\\s*([0-9]+).*", "\\1", rss_line[[1L]])
  )
  result
}

results <- do.call(rbind, lapply(c(16L, 32L, 64L), run_worker))
write.table(results, row.names = FALSE, sep = "\t", quote = FALSE)
