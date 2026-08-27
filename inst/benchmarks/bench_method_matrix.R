#!/usr/bin/env Rscript
# Pinned speed/accuracy matrix for the cluster4d method family.
#
# Produces the baseline table referenced by the performance epic. Every number a
# gate is written against must come from this script so the comparison is like
# for like. Emits full provenance (git SHA, sessionInfo, thread count, timing
# policy) alongside the results.
#
# Usage: Rscript inst/benchmarks/bench_method_matrix.R [outfile.csv]
# Set METHOD_MATRIX_METHODS to a comma-separated subset for a focused replay.

suppressMessages({
  library(neurocluster)
  library(neuroim2)
})

## ---- policy ---------------------------------------------------------------
# Warm-up: one discarded call per method (JIT, lazy-load, BLAS thread spin-up).
# Repetitions: 3 timed reps when the warm-up came in under 5s, otherwise 1.
# Slow methods are the ones being fixed; a single rep is enough to size them.
# Reported: minimum elapsed over timed reps (least contaminated by other load).
# Timeout: per method, via setTimeLimit(elapsed=). Recorded as NA on trip.
TIMEOUT_S <- 300
FAST_CUTOFF <- 5
REPS_FAST <- 3L
REPS_SLOW <- 1L

METHODS <- c(
  "supervoxels", "snic", "slic", "corr_slic", "brs_slic", "slice_msf",
  "flash3d", "g3s", "rena", "rena_plus", "mcl", "acsc", "commute"
)
method_filter <- Sys.getenv("METHOD_MATRIX_METHODS", "")
if (nzchar(method_filter)) {
  requested_methods <- trimws(strsplit(method_filter, ",", fixed = TRUE)[[1L]])
  unknown_methods <- setdiff(requested_methods, METHODS)
  if (length(unknown_methods)) {
    stop("Unknown METHOD_MATRIX_METHODS: ", paste(unknown_methods, collapse = ", "))
  }
  METHODS <- requested_methods
}

## ---- fixtures -------------------------------------------------------------
speed_fixture <- function() {
  set.seed(11)
  dx <- 40
  dy <- 40
  dz <- 20
  n_time <- 120
  n_voxels <- dx * dy * dz
  rank <- 10
  x <- matrix(rnorm(n_voxels * rank), n_voxels, rank) %*%
    t(matrix(rnorm(n_time * rank), n_time, rank)) +
    matrix(rnorm(n_voxels * n_time, sd = 0.5), n_voxels, n_time)
  voxel_spacing <- c(3, 3, 3)
  list(
    vec = NeuroVec(
      array(x, c(dx, dy, dz, n_time)),
      NeuroSpace(c(dx, dy, dz, n_time), spacing = voxel_spacing)
    ),
    mask = NeuroVol(
      array(1, c(dx, dy, dz)),
      NeuroSpace(c(dx, dy, dz), spacing = voxel_spacing)
    ),
    x = x,
    n_voxels = n_voxels,
    n_time = n_time,
    target_k = 300L
  )
}

# Blocky parcels on a slice-friendly volume: independent of the package's own
# cube fixture, and shaped so a slice-wise method is not disadvantaged.
accuracy_fixture <- function() {
  set.seed(9)
  dims <- c(24, 24, 6)
  n_time <- 80
  grid <- as.matrix(expand.grid(
    seq_len(dims[[1L]]), seq_len(dims[[2L]]), seq_len(dims[[3L]])
  ))
  truth <- as.integer(factor(as.integer(interaction(
    cut(grid[, 1L], 4, labels = FALSE),
    cut(grid[, 2L], 3, labels = FALSE),
    cut(grid[, 3L], 2, labels = FALSE)
  ))))
  signal <- matrix(rnorm(max(truth) * n_time), max(truth), n_time)
  x <- signal[truth, ] +
    matrix(rnorm(nrow(grid) * n_time, sd = 0.2), nrow(grid), n_time)
  voxel_spacing <- c(3, 3, 3)
  list(
    vec = NeuroVec(
      array(x, c(dims, n_time)),
      NeuroSpace(c(dims, n_time), spacing = voxel_spacing)
    ),
    mask = NeuroVol(
      array(1, dims), NeuroSpace(dims, spacing = voxel_spacing)
    ),
    truth = truth,
    target_k = max(truth)
  )
}

## ---- provenance -----------------------------------------------------------
git <- function(...) {
  tryCatch(
    system2("git", c(...), stdout = TRUE, stderr = FALSE),
    error = function(error) NA_character_
  )
}
cat("# neurocluster method matrix\n")
cat("# generated : ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\n", sep = "")
cat("# git sha   : ", paste(git("rev-parse", "HEAD"), collapse = ""), "\n", sep = "")
cat("# git dirty : ", if (length(git("status", "--porcelain"))) "YES" else "no", "\n", sep = "")
cat("# pkg ver   : ", as.character(packageVersion("neurocluster")), "\n", sep = "")
cat("# R         : ", R.version.string, "\n", sep = "")
cat("# platform  : ", R.version$platform, "\n", sep = "")
cat("# BLAS      : ", La_library(), "\n", sep = "")
cat(
  "# threads   : RcppParallel=",
  tryCatch(RcppParallel::defaultNumThreads(), error = function(error) NA),
  " OMP_NUM_THREADS=", Sys.getenv("OMP_NUM_THREADS", "<unset>"), "\n",
  sep = ""
)
cat(
  "# policy    : warmup=1 discarded; reps=", REPS_FAST,
  " if warmup<", FAST_CUTOFF, "s else ", REPS_SLOW,
  "; report=min elapsed; timeout=", TIMEOUT_S, "s\n", sep = ""
)
cat("# methods   : ", paste(METHODS, collapse = ","), "\n", sep = "")

speed <- speed_fixture()
accuracy <- accuracy_fixture()
scaled_features <- stats::scale(speed$x)
wss <- function(labels) {
  centers <- rowsum(scaled_features, labels) / as.vector(tabulate(labels))
  sum((scaled_features - centers[labels, , drop = FALSE])^2)
}

run_once <- function(method, fixture, target_k) {
  setTimeLimit(elapsed = TIMEOUT_S, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  set.seed(1)
  suppressWarnings(suppressMessages(cluster4d(
    fixture$vec, fixture$mask,
    n_clusters = target_k, method = method, verbose = FALSE
  )))
}

rows <- list()
cat(sprintf(
  "\n%-13s %9s %6s %6s %9s %12s %9s %9s\n",
  "method", "sec", "reps", "K_in", "K_out", "WSS", "min_size", "ARI"
))
for (method in METHODS) {
  warmup <- tryCatch(
    system.time(run_once(method, speed, speed$target_k))[["elapsed"]],
    error = function(error) NA_real_
  )
  if (is.na(warmup)) {
    cat(sprintf("%-13s %9s\n", method, "TIMEOUT/ERR"))
    rows[[method]] <- data.frame(
      method = method, sec = NA_real_, reps = 0L, natural_k = NA_integer_,
      repaired_k = NA_integer_, wss = NA_real_, min_size = NA_integer_,
      ari = NA_real_, cluster_sizes = NA_character_
    )
    next
  }
  repetitions <- if (warmup < FAST_CUTOFF) REPS_FAST else REPS_SLOW
  elapsed <- rep(NA_real_, repetitions)
  result <- NULL
  for (iteration in seq_len(repetitions)) {
    elapsed[[iteration]] <- tryCatch(
      system.time({
        result <- run_once(method, speed, speed$target_k)
      })[["elapsed"]],
      error = function(error) NA_real_
    )
  }
  labels <- as.integer(result$cluster)
  sizes <- tabulate(labels)
  trace <- result$metadata$exact_k_repair
  natural_k <- if (!is.null(trace$natural_k)) trace$natural_k else NA_integer_
  observed_ari <- tryCatch({
    accuracy_result <- run_once(method, accuracy, accuracy$target_k)
    clustering_accuracy(accuracy_result$cluster, accuracy$truth)$ari
  }, error = function(error) NA_real_)
  rows[[method]] <- data.frame(
    method = method,
    sec = min(elapsed, na.rm = TRUE),
    reps = repetitions,
    natural_k = natural_k,
    repaired_k = length(unique(labels)),
    wss = wss(labels),
    min_size = min(sizes[sizes > 0L]),
    ari = observed_ari,
    cluster_sizes = paste(sizes, collapse = ","),
    stringsAsFactors = FALSE
  )
  cat(sprintf(
    "%-13s %9.2f %6d %6s %9d %12.0f %9d %9.4f\n",
    method, min(elapsed, na.rm = TRUE), repetitions,
    ifelse(is.na(natural_k), "NA", as.character(natural_k)),
    length(unique(labels)), wss(labels), min(sizes[sizes > 0L]), observed_ari
  ))
}

output <- do.call(rbind, rows)
destination <- if (length(commandArgs(TRUE))) {
  commandArgs(TRUE)[[1L]]
} else {
  "method_matrix.csv"
}
utils::write.csv(output, destination, row.names = FALSE)
cat("\n# wrote ", destination, "\n", sep = "")
