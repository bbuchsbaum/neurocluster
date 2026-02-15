# A/B benchmark for int8 dot-product kernel (NEON dotprod vs portable).
#
# This runs two fresh R processes because the kernel selection can be cached on
# first use (and because we want a clean comparison).
#
# Usage (from repo root):
#   Rscript inst/benchmarks/bench_dot_kernel_ab.R
#
# Optional env vars:
#   NEUROCLUSTER_DOT_AB_REPS=5
#   NEUROCLUSTER_DOT_AB_D=128
#   NEUROCLUSTER_DOT_AB_K=600
#   NEUROCLUSTER_DOT_AB_T=60
#   NEUROCLUSTER_DOT_AB_DIMS=48,48,16
#   NEUROCLUSTER_DOT_AB_MAX_ITER=6
#   NEUROCLUSTER_DOT_AB_ASSIGN_STRIDE=1
#   NEUROCLUSTER_DOT_AB_SKETCH_REPEATS=1
#   NEUROCLUSTER_DOT_AB_ALPHA=0.1
#   NEUROCLUSTER_DOT_AB_SEED=7
#   NEUROCLUSTER_DOT_AB_THREADS=1
#   NEUROCLUSTER_DOT_AB_USE_PKGLOAD=1   # set to 0 to prefer library(neurocluster)

parse_int <- function(x, default) {
  v <- suppressWarnings(as.integer(x))
  if (length(v) != 1L || is.na(v)) default else v
}

parse_num <- function(x, default) {
  v <- suppressWarnings(as.numeric(x))
  if (length(v) != 1L || is.na(v) || !is.finite(v)) default else v
}

parse_int_vec <- function(x, default) {
  if (!nzchar(x)) return(default)
  parts <- strsplit(x, "[, ]+", perl = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  v <- suppressWarnings(as.integer(parts))
  if (any(is.na(v)) || length(v) == 0L) return(default)
  v
}

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (!length(file_arg)) return(NA_character_)
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)
}

run_child <- function(out_path, force_portable) {
  script_path <- get_script_path()
  if (is.na(script_path) || !file.exists(script_path)) {
    stop("Could not resolve script path for child invocation.")
  }
  rscript <- Sys.which("Rscript")
  if (!nzchar(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
  if (!file.exists(rscript)) stop("Could not find Rscript executable.")

  env <- c(
    sprintf("NEUROCLUSTER_FORCE_PORTABLE_DOT_I8=%s", if (force_portable) "1" else "0"),
    sprintf("NEUROCLUSTER_DOT_AB_REPS=%s", Sys.getenv("NEUROCLUSTER_DOT_AB_REPS", "5")),
    sprintf("NEUROCLUSTER_DOT_AB_D=%s", Sys.getenv("NEUROCLUSTER_DOT_AB_D", "128")),
    sprintf("NEUROCLUSTER_DOT_AB_K=%s", Sys.getenv("NEUROCLUSTER_DOT_AB_K", "600")),
    sprintf("NEUROCLUSTER_DOT_AB_T=%s", Sys.getenv("NEUROCLUSTER_DOT_AB_T", "60")),
    sprintf("NEUROCLUSTER_DOT_AB_DIMS=%s", Sys.getenv("NEUROCLUSTER_DOT_AB_DIMS", "48,48,16")),
    sprintf("NEUROCLUSTER_DOT_AB_MAX_ITER=%s", Sys.getenv("NEUROCLUSTER_DOT_AB_MAX_ITER", "6")),
    sprintf("NEUROCLUSTER_DOT_AB_ASSIGN_STRIDE=%s", Sys.getenv("NEUROCLUSTER_DOT_AB_ASSIGN_STRIDE", "1")),
    sprintf("NEUROCLUSTER_DOT_AB_SKETCH_REPEATS=%s", Sys.getenv("NEUROCLUSTER_DOT_AB_SKETCH_REPEATS", "1")),
    sprintf("NEUROCLUSTER_DOT_AB_ALPHA=%s", Sys.getenv("NEUROCLUSTER_DOT_AB_ALPHA", "0.1")),
    sprintf("NEUROCLUSTER_DOT_AB_SEED=%s", Sys.getenv("NEUROCLUSTER_DOT_AB_SEED", "7")),
    sprintf("NEUROCLUSTER_DOT_AB_THREADS=%s", Sys.getenv("NEUROCLUSTER_DOT_AB_THREADS", "1")),
    sprintf("NEUROCLUSTER_DOT_AB_USE_PKGLOAD=%s", Sys.getenv("NEUROCLUSTER_DOT_AB_USE_PKGLOAD", "1"))
  )

  args <- c(script_path, "--child", out_path)
  status <- system2(rscript, args = args, env = env, stdout = "", stderr = "")
  if (!identical(status, 0L)) {
    stop(sprintf("Child process failed (status=%s).", status))
  }
  if (!file.exists(out_path)) {
    stop("Child process did not write output file.")
  }
  readRDS(out_path)
}

child_main <- function(out_path) {
  use_pkgload <- Sys.getenv("NEUROCLUSTER_DOT_AB_USE_PKGLOAD", "1") != "0"
  suppressPackageStartupMessages({
    if (use_pkgload && file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
      pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
    } else {
      library(neurocluster)
    }
  })

  reps <- parse_int(Sys.getenv("NEUROCLUSTER_DOT_AB_REPS", "5"), 5L)
  d <- parse_int(Sys.getenv("NEUROCLUSTER_DOT_AB_D", "128"), 128L)
  K <- parse_int(Sys.getenv("NEUROCLUSTER_DOT_AB_K", "600"), 600L)
  T <- parse_int(Sys.getenv("NEUROCLUSTER_DOT_AB_T", "60"), 60L)
  dims <- parse_int_vec(Sys.getenv("NEUROCLUSTER_DOT_AB_DIMS", "48,48,16"), c(48L, 48L, 16L))
  if (length(dims) != 3L) dims <- c(48L, 48L, 16L)

  max_iter <- parse_int(Sys.getenv("NEUROCLUSTER_DOT_AB_MAX_ITER", "6"), 6L)
  assign_stride <- parse_int(Sys.getenv("NEUROCLUSTER_DOT_AB_ASSIGN_STRIDE", "1"), 1L)
  sketch_repeats <- parse_int(Sys.getenv("NEUROCLUSTER_DOT_AB_SKETCH_REPEATS", "1"), 1L)
  alpha <- parse_num(Sys.getenv("NEUROCLUSTER_DOT_AB_ALPHA", "0.1"), 0.1)
  seed <- parse_int(Sys.getenv("NEUROCLUSTER_DOT_AB_SEED", "7"), 7L)
  n_threads <- parse_int(Sys.getenv("NEUROCLUSTER_DOT_AB_THREADS", "1"), 1L)

  kernel <- neurocluster_simd_info()$dot_i8_kernel

  set.seed(1)
  N <- prod(dims)
  feat <- matrix(rnorm(N * T), nrow = N, ncol = T)
  mask_lin_idx <- 0:(N - 1L)
  corrslic_core <- getFromNamespace("corrslic_core", "neurocluster")

  # Warmup.
  invisible(corrslic_core(
    feat = feat,
    mask_lin_idx = as.integer(mask_lin_idx),
    dims = as.integer(dims),
    K = as.integer(K),
    d = as.integer(d),
    sketch_repeats = as.integer(sketch_repeats),
    alpha = as.numeric(alpha),
    max_iter = as.integer(max_iter),
    seed = as.integer(seed),
    assign_stride = as.integer(assign_stride),
    quantize_assign = TRUE,
    embed_basis = "hash",
    whiten_embed = FALSE,
    refine_exact_iters = 0L,
    refine_boundary_only = TRUE,
    refine_stride = 1L,
    refine_alpha = -1.0,
    connectivity = 6L,
    min_size = 0L,
    n_threads = as.integer(n_threads),
    verbose = FALSE
  ))

  times <- numeric(reps)
  for (i in seq_len(reps)) {
    invisible(gc())
    times[[i]] <- unname(system.time(invisible(corrslic_core(
      feat = feat,
      mask_lin_idx = as.integer(mask_lin_idx),
      dims = as.integer(dims),
      K = as.integer(K),
      d = as.integer(d),
      sketch_repeats = as.integer(sketch_repeats),
      alpha = as.numeric(alpha),
      max_iter = as.integer(max_iter),
      seed = as.integer(seed),
      assign_stride = as.integer(assign_stride),
      quantize_assign = TRUE,
      embed_basis = "hash",
      whiten_embed = FALSE,
      refine_exact_iters = 0L,
      refine_boundary_only = TRUE,
      refine_stride = 1L,
      refine_alpha = -1.0,
      connectivity = 6L,
      min_size = 0L,
      n_threads = as.integer(n_threads),
      verbose = FALSE
    )))[["elapsed"]])
  }

  res <- list(
    kernel = kernel,
    force_portable = Sys.getenv("NEUROCLUSTER_FORCE_PORTABLE_DOT_I8", ""),
    times_sec = times,
    median_sec = unname(stats::median(times)),
    mean_sec = unname(mean(times)),
    min_sec = unname(min(times)),
    max_sec = unname(max(times)),
    params = list(
      reps = reps,
      dims = dims,
      n_vox = N,
      n_time = T,
      K = K,
      d = d,
      alpha = alpha,
      max_iter = max_iter,
      assign_stride = assign_stride,
      sketch_repeats = sketch_repeats,
      seed = seed,
      n_threads = n_threads
    )
  )

  saveRDS(res, out_path)
  cat(sprintf("RESULT kernel=%s median=%.6f mean=%.6f reps=%d\n",
              res$kernel, res$median_sec, res$mean_sec, reps))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2L && identical(args[[1]], "--child")) {
    child_main(args[[2]])
    return(invisible(NULL))
  }

  cat("== neurocluster int8 dot-kernel A/B ==\n")
  cat("This runs two fresh R processes:\n")
  cat("  1) default (NEON dotprod if available)\n")
  cat("  2) forced portable (NEUROCLUSTER_FORCE_PORTABLE_DOT_I8=1)\n\n")

  out_neon <- tempfile("dot_ab_neon_", fileext = ".rds")
  out_port <- tempfile("dot_ab_portable_", fileext = ".rds")

  neon <- run_child(out_neon, force_portable = FALSE)
  port <- run_child(out_port, force_portable = TRUE)

  cat("\n== Results ==\n")
  cat(sprintf("default kernel:  %s | median %.6fs | mean %.6fs\n",
              neon$kernel, neon$median_sec, neon$mean_sec))
  cat(sprintf("forced kernel:   %s | median %.6fs | mean %.6fs\n",
              port$kernel, port$median_sec, port$mean_sec))

  if (is.finite(neon$median_sec) && is.finite(port$median_sec) &&
      neon$median_sec > 0 && port$median_sec > 0) {
    speedup <- port$median_sec / neon$median_sec
    pct <- (1 - (neon$median_sec / port$median_sec)) * 100
    cat(sprintf("\nNEON speedup vs portable (median): %.3fx (%.1f%% faster)\n", speedup, pct))
  }

  if (identical(neon$kernel, port$kernel)) {
    cat("\nNote: both runs reported the same kernel; you may not have an ISA-specific kernel on this build.\n")
  }
}

main()

