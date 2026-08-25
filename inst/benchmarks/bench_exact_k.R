#!/usr/bin/env Rscript

# Fixed-case comparison of the topology-preserving exact-K engine with the
# pinned pre-remediation hclust/k-means helper. The parent process uses
# /usr/bin/time -l to capture peak resident memory for each isolated worker.

args <- commandArgs(trailingOnly = TRUE)

old_force_exact_k <- function(labels, feature_mat, K_target) {
  labels <- as.integer(labels)
  K_target <- max(1L, min(K_target, nrow(feature_mat) - 1L))
  relabel <- function(lbl) as.integer(match(lbl, sort(unique(lbl))))
  labels <- relabel(labels)
  k_curr <- length(unique(labels))
  if (k_curr == K_target) return(labels)
  centroids <- function(lbls) {
    rowsum(feature_mat, lbls) / as.numeric(table(lbls))
  }
  if (k_curr > K_target) {
    cut <- stats::cutree(
      stats::hclust(stats::dist(centroids(labels)), method = "ward.D2"),
      k = K_target
    )
    return(relabel(cut[labels]))
  }
  set.seed(1)
  relabel(stats::kmeans(
    feature_mat, centers = K_target, iter.max = 50, nstart = 5
  )$cluster)
}

make_case <- function(case) {
  dims <- c(20L, 20L, 4L)
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
  if (case == "over_target") {
    labels <- interaction(
      (coords[, 1] - 1L) %/% 4L,
      (coords[, 2] - 1L) %/% 4L,
      (coords[, 3] - 1L) %/% 2L,
      drop = TRUE
    )
  } else {
    labels <- interaction(
      (coords[, 1] - 1L) %/% 10L,
      (coords[, 2] - 1L) %/% 10L,
      (coords[, 3] - 1L) %/% 2L,
      drop = TRUE
    )
  }
  list(
    mask = mask,
    features = features,
    labels = as.integer(labels),
    target = 25L
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
  method <- args[[3]]
  repetitions <- as.integer(args[[4]])
  fixture <- make_case(case_name)
  graph <- neurocluster:::.exact_k_graph(fixture$mask, 6L)
  elapsed <- numeric(repetitions)
  result <- NULL
  for (i in seq_len(repetitions)) {
    elapsed[[i]] <- system.time({
      result <- if (method == "candidate") {
        neurocluster:::force_exact_k(
          fixture$labels, fixture$features, fixture$target,
          fixture$mask, 6L
        )
      } else {
        old_force_exact_k(
          fixture$labels, fixture$features, fixture$target
        )
      }
    })[["elapsed"]]
  }
  output <- data.frame(
    case = case_name,
    method = method,
    n_voxels = length(fixture$labels),
    n_edges = nrow(graph$edges),
    initial_k = length(unique(fixture$labels)),
    target_k = fixture$target,
    repetitions = repetitions,
    median_elapsed_s = stats::median(elapsed),
    p95_elapsed_s = stats::quantile(elapsed, 0.95, names = FALSE),
    exact_k = length(unique(result)) == fixture$target,
    topology_ok = topology_ok(result, fixture$mask),
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

run_worker <- function(case, method, repetitions = 5L) {
  stdout <- tempfile("exact-k-bench-", fileext = ".tsv")
  stderr <- tempfile("exact-k-time-", fileext = ".txt")
  on.exit(unlink(c(stdout, stderr)), add = TRUE)
  status <- system2(
    time_bin,
    c("-l", rscript, script, "--worker", case, method, repetitions),
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

results <- do.call(rbind, lapply(c("over_target", "under_target"), function(case) {
  rbind(
    run_worker(case, "baseline"),
    run_worker(case, "candidate")
  )
}))
write.table(results, row.names = FALSE, sep = "\t", quote = FALSE)
