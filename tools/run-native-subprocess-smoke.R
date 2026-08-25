#!/usr/bin/env Rscript

# Run every clustering-related native entry point in its own R process. The
# parent enforces a timeout per entry point so a native abort or hang cannot
# take down the rest of the smoke matrix.

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
script <- if (length(file_arg) == 1L) {
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
} else {
  normalizePath("tools/run-native-subprocess-smoke.R", mustWork = TRUE)
}

option_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

root <- option_value("root")
if (is.null(root)) root <- file.path(dirname(script), "..")
root <- normalizePath(root, mustWork = TRUE)
case_name <- option_value("case")
timeout <- as.integer(option_value("timeout", "60"))
if (is.na(timeout) || timeout < 1L) stop("--timeout must be a positive integer")
jobs <- as.integer(option_value("jobs", "4"))
if (is.na(jobs) || jobs < 1L) stop("--jobs must be a positive integer")

native <- function(name) get(name, envir = asNamespace("neurocluster"), inherits = FALSE)

fixture <- function() {
  dims <- c(2L, 2L, 2L)
  n <- prod(dims)
  coords <- cbind(
    rep(0:1, 4),
    rep(rep(0:1, each = 2), 2),
    rep(0:1, each = 4)
  )
  features <- outer(seq_len(n), seq_len(4L), function(i, j) {
    sin(i * (j + 0.25)) + i / 10 + j / 20
  })
  features <- sweep(features, 1L, sqrt(rowSums(features^2)), "/")
  labels <- rep(1:2, each = n / 2L)
  groups0 <- labels - 1L
  nn0 <- cbind(seq_len(n) - 1L, c(seq_len(n - 1L), 0L))
  nn1 <- nn0 + 1L
  nn_dist <- matrix(c(rep(0, n), rep(0.25, n)), nrow = n)
  data_centers <- vapply(split(seq_len(n), labels), function(idx) {
    colMeans(features[idx, , drop = FALSE])
  }, numeric(ncol(features)))
  coord_centers <- vapply(split(seq_len(n), labels), function(idx) {
    colMeans(coords[idx, , drop = FALSE])
  }, numeric(ncol(coords)))

  list(
    dims = dims,
    n = n,
    coords = coords,
    coords_t = t(coords),
    features = features,
    data = t(features),
    ts = t(features),
    labels = labels,
    groups0 = groups0,
    nn0 = nn0,
    nn1 = nn1,
    nn_dist = nn_dist,
    data_centers = data_centers,
    coord_centers = coord_centers,
    mask_lin0 = seq_len(n) - 1L,
    mask_lin1 = seq_len(n),
    mask = array(1L, dims)
  )
}

adjacency_fixture <- function(fx) {
  native("build_grid_adjacency_cpp")(fx$mask_lin1, fx$dims, 6L)
}

candidate_fixture <- function(fx) {
  replicate(fx$n, 1:2, simplify = FALSE)
}

slice_fixture <- function(fx) {
  native("slice_msf_runwise")(
    fx$ts, as.integer(fx$mask), fx$dims,
    r = 2L, fh_scale = 0.32, min_size = 1L, nbhd = 4L,
    rows_are_time = TRUE, gamma = 1.5, voxel_dim = c(1, 1, 1)
  )
}

snic_fixture <- function(fx) {
  valid_coords <- matrix(as.integer(fx$coords), ncol = 3L)
  list(
    L = array(0L, fx$dims),
    mask = fx$mask,
    centroids = fx$coords[c(1L, fx$n), , drop = FALSE],
    centroid_idx = c(0L, fx$n - 1L),
    valid_coords = valid_coords,
    norm_coords = fx$coords,
    vecmat = fx$data,
    K = 2L,
    s = 1,
    compactness = 1,
    lookup = array(seq_len(fx$n) - 1L, fx$dims)
  )
}

cases <- list(
  refine_boundaries_cpp = function(fx) {
    native("refine_boundaries_cpp")(
      fx$labels, fx$features, fx$nn1, seq_len(fx$n), max_iter = 1L
    )
  },
  find_boundary_voxels_cpp = function(fx) {
    native("find_boundary_voxels_cpp")(fx$labels, fx$nn1)
  },
  compute_centroids_parallel_fast = function(fx) {
    native("compute_centroids_parallel_fast")(
      fx$groups0, fx$data, fx$coords_t, 2L
    )
  },
  corrslic_exact_scores_cpp = function(fx) {
    native("corrslic_exact_scores_cpp")(
      fx$features, t(fx$data_centers), fx$coords, t(fx$coord_centers),
      spatial_weight = 0.5, spatial_scale = 1
    )
  },
  corrslic_core = function(fx) {
    native("corrslic_core")(
      fx$features, fx$mask_lin0, fx$dims, 2L,
      d = 8L, sketch_repeats = 1L, max_iter = 1L,
      connectivity = 6L, n_threads = 1L
    )
  },
  brs_slic_core = function(fx) {
    native("brs_slic_core")(
      fx$features, fx$mask_lin0, fx$dims, 2L,
      d = 8L, sketch_repeats = 1L, coarse_iter = 1L,
      boundary_passes = 1L, connectivity = 6L, n_threads = 1L
    )
  },
  correlation_gradient_cpp = function(fx) {
    image <- array(as.numeric(fx$features), c(fx$dims, ncol(fx$features)))
    native("correlation_gradient_cpp")(image, array(as.numeric(fx$mask), fx$dims))
  },
  build_grid_adjacency_cpp = function(fx) adjacency_fixture(fx),
  rena_rnn_coarse_cpp = function(fx) {
    native("rena_rnn_coarse_cpp")(
      fx$features, adjacency_fixture(fx), numeric(), 4L, lambda = 0, max_iter = 2L
    )
  },
  ward_on_supervoxels_cpp = function(fx) {
    native("ward_on_supervoxels_cpp")(
      fx$features, adjacency_fixture(fx), rep(1L, fx$n), 2L
    )
  },
  compute_scores = function(fx) {
    native("compute_scores")(
      fx$labels, fx$coords_t, fx$data_centers, fx$coord_centers,
      fx$data, 1, 1
    )
  },
  best_candidate = function(fx) {
    native("best_candidate")(
      candidate_fixture(fx), fx$labels, fx$coords_t,
      fx$data_centers, fx$coord_centers, fx$data, 1, 1, 0.5
    )
  },
  find_candidates = function(fx) {
    native("find_candidates")(fx$nn0, fx$nn_dist, fx$labels, 1)
  },
  best_candidate_parallel = function(fx) {
    native("best_candidate_parallel")(
      candidate_fixture(fx), fx$labels, fx$coords_t,
      fx$data_centers, fx$coord_centers, fx$data, 1, 1, 0.5, grain_size = 1L
    )
  },
  best_candidate_sequential = function(fx) {
    native("best_candidate_sequential")(
      candidate_fixture(fx), fx$labels, fx$coords_t,
      fx$data_centers, fx$coord_centers, fx$data, 1, 1, 0.5
    )
  },
  heat_kernel = function(fx) native("heat_kernel")(fx$data[, 1L], fx$data[, 2L], 1),
  normalized_heat_kernel = function(fx) {
    native("normalized_heat_kernel")(fx$data[, 1L], fx$data[, 2L], 1)
  },
  flash3d_reseed_empty_cpp = function(fx) {
    native("flash3d_reseed_empty_cpp")(
      c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L), seq_len(8L), 4L
    )
  },
  flash3d_supervoxels_cpp = function(fx) {
    native("flash3d_supervoxels_cpp")(
      fx$ts, fx$mask_lin1, fx$dims, 2L, c(0.5, 0.5, 0),
      rounds = 1L, bits = 64L, dctM = 4L, vox_scale = c(1, 1, 1)
    )
  },
  fused_assignment = function(fx) {
    native("fused_assignment")(
      fx$nn0, fx$nn_dist, fx$labels, fx$coords_t,
      fx$data_centers, fx$coord_centers, fx$data, 1, 1, 1, 0.5
    )
  },
  fused_assignment_parallel = function(fx) {
    native("fused_assignment_parallel")(
      fx$nn0, fx$nn_dist, fx$labels, fx$coords_t,
      fx$data_centers, fx$coord_centers, fx$data, 1, 1, 1, 0.5,
      grain_size = 1L
    )
  },
  compute_centroids_parallel = function(fx) {
    native("compute_centroids_parallel")(
      fx$labels, fx$data, fx$coords_t, 2L, grain_size = 1L
    )
  },
  fused_assignment_parallel_binned = function(fx) {
    output <- NULL
    for (iteration in seq_len(100L)) {
      output <- native("fused_assignment_parallel_binned")(
        fx$nn0, fx$nn_dist, fx$groups0, fx$coords_t,
        fx$data_centers, fx$coord_centers,
        tabulate(fx$groups0 + 1L, nbins = 2L), fx$data, 1, 1, 1, 0.5,
        grain_size = 1L, window_factor = 2, bin_expand = 1L
      )
      stopifnot(length(output) == fx$n, all(output %in% 0:1))
    }
    output
  },
  fused_assignment_binned = function(fx) {
    native("fused_assignment_binned")(
      fx$nn0, fx$nn_dist, fx$groups0, fx$coords_t,
      fx$data_centers, fx$coord_centers,
      tabulate(fx$groups0 + 1L, nbins = 2L), fx$data, 1, 1, 1, 0.5,
      window_factor = 2, bin_expand = 1L
    )
  },
  calculate_local_gradient = function(fx) {
    native("calculate_local_gradient")(fx$data, fx$nn1)
  },
  g3s_propagate_cpp = function(fx) {
    native("g3s_propagate_cpp")(
      fx$data, c(1L, fx$n),
      fx$nn1[, 2L, drop = FALSE], fx$nn_dist[, 2L, drop = FALSE],
      0.5, 1
    )
  },
  refine_boundaries_g3s_cpp = function(fx) {
    native("refine_boundaries_g3s_cpp")(
      fx$labels, fx$data, fx$coords, fx$nn1[, 2L, drop = FALSE],
      0.5, 1, 1L
    )
  },
  mcl_prune_sparse_cpp = function(fx) {
    native("mcl_prune_sparse_cpp")(
      c(0L, 2L, 4L), c(0L, 1L, 0L, 1L), c(0.8, 0.2, 0.3, 0.7),
      2L, 2L, 0
    )
  },
  compute_masked_distances_cpp = function(fx) {
    native("compute_masked_distances_cpp")(fx$data, 1:7, 2:8)
  },
  find_1nn_subgraph_cpp = function(fx) {
    native("find_1nn_subgraph_cpp")(fx$n, 1:7, 2:8, rep(1, 7))
  },
  find_connected_components_cpp = function(fx) {
    native("find_connected_components_cpp")(fx$n, c(1:7, 6L))
  },
  aggregate_features_cpp = function(fx) {
    native("aggregate_features_cpp")(fx$data, fx$groups0, 2L)
  },
  aggregate_coords_cpp = function(fx) {
    native("aggregate_coords_cpp")(fx$coords, fx$groups0, 2L)
  },
  contract_graph_cpp = function(fx) {
    native("contract_graph_cpp")(1:7, 2:8, fx$groups0, 2L)
  },
  prune_edges_for_k_cpp = function(fx) {
    native("prune_edges_for_k_cpp")(fx$n, c(1:7, 6L), rep(1, fx$n), 2L)
  },
  slic4d_core = function(fx) {
    standard <- native("slic4d_core")(
      fx$features, fx$coords, fx$mask_lin0, fx$dims, c(1, 1, 1), 2L,
      compactness = 10, max_iter = 1L, n_threads = 1L,
      seed_method = "mask_grid", enforce_connectivity = FALSE
    )
    one_cluster <- native("slic4d_core")(
      fx$features, fx$coords, fx$mask_lin0, fx$dims, c(1, 1, 1), 1L,
      compactness = 0, max_iter = 1L, n_threads = 1L,
      seed_method = "mask_grid", enforce_connectivity = TRUE,
      connectivity = 6L, strict_connectivity = TRUE, preserve_k = TRUE
    )
    checker_dims <- c(6L, 6L, 1L)
    checker_idx <- seq_len(prod(checker_dims))
    checker_coords <- arrayInd(checker_idx, checker_dims)
    checker <- matrix((rowSums(checker_coords) %% 2L) * 10, ncol = 1L)
    fragmented <- native("slic4d_core")(
      checker, checker_coords, checker_idx - 1L, checker_dims, c(1, 1, 1), 2L,
      compactness = 0, max_iter = 3L, step_mm = 1, n_threads = 1L,
      seed_method = "mask_grid", enforce_connectivity = TRUE,
      connectivity = 6L, strict_connectivity = TRUE, preserve_k = FALSE,
      topup_iters = 0L
    )
    stopifnot(
      length(standard$labels) == fx$n,
      identical(one_cluster$actual_k, 1L),
      fragmented$connectivity_changed,
      nrow(fragmented$center_feats) == fragmented$actual_k,
      all(is.finite(fragmented$center_feats)),
      all(is.finite(fragmented$center_coords))
    )
    list(standard = standard, one_cluster = one_cluster, fragmented = fragmented)
  },
  slice_msf_runwise = function(fx) slice_fixture(fx),
  slice_fuse_consensus = function(fx) {
    run <- slice_fixture(fx)
    native("slice_fuse_consensus")(
      list(run, run), fx$dims, nbhd = 4L, fh_scale = 0.3,
      min_size = 1L, use_features = FALSE, voxel_dim = c(1, 1, 1)
    )
  },
  update_centroid_online = function(fx) {
    native("update_centroid_online")(
      list(x = c(0, 0, 0), c = fx$data[, 1L], label = 1L, n = 1L),
      fx$coords[2L, ], fx$data[, 2L]
    )
  },
  snic_main = function(fx) {
    s <- snic_fixture(fx)
    native("snic_main")(
      s$L, s$mask, s$centroids, s$centroid_idx, s$valid_coords,
      s$norm_coords, s$vecmat, s$K, s$s, s$compactness, s$lookup
    )
  },
  snic_main_optimized = function(fx) {
    s <- snic_fixture(fx)
    native("snic_main_optimized")(
      s$L, s$mask, s$centroids, s$centroid_idx, s$valid_coords,
      s$norm_coords, s$vecmat, s$K, s$s, s$compactness, s$lookup
    )
  },
  compute_boundaryscore_3d_cpp = function(fx) {
    native("compute_boundaryscore_3d_cpp")(array(fx$labels, fx$dims), fx$mask)
  },
  detect_boundaries_2d_cpp = function(fx) {
    native("detect_boundaries_2d_cpp")(array(fx$labels, fx$dims), fx$mask)
  }
)

load_current_package <- function() {
  use_installed <- identical(tolower(Sys.getenv("NEUROCLUSTER_USE_INSTALLED")), "true")
  if (use_installed) {
    suppressPackageStartupMessages(library(neurocluster))
  } else {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("pkgload is required for source-tree native smoke tests")
    }
    pkgload::load_all(root, recompile = FALSE, quiet = TRUE)
  }
}

if (!is.null(case_name)) {
  if (!case_name %in% names(cases)) stop("Unknown smoke case: ", case_name)
  load_current_package()
  value <- cases[[case_name]](fixture())
  if (is.null(value)) stop(case_name, " returned NULL")
  cat("native-smoke-ok:", case_name, "\n")
  quit(save = "no", status = 0L)
}

export_lines <- readLines(file.path(root, "R", "RcppExports.R"), warn = FALSE)
export_match <- regexec("^([A-Za-z][A-Za-z0-9_.]*) <- function\\(", export_lines)
declared <- vapply(regmatches(export_lines, export_match), function(x) {
  if (length(x) == 2L) x[[2L]] else NA_character_
}, character(1L))
declared <- declared[!is.na(declared)]

non_clustering <- c(
  "normalize_volumes_cpp", "detrend_time_cpp", "normalize_detrend_cpp",
  "make_dct_basis", "make_poly_basis", "detrend_basis_cpp",
  "detrend_poly_cpp", "detrend_dct_cpp", "simd_info_cpp"
)
required <- setdiff(declared, non_clustering)
missing <- setdiff(required, names(cases))
extra <- setdiff(names(cases), required)
if (length(missing) || length(extra)) {
  stop(
    "Native smoke coverage mismatch. Missing: ", paste(missing, collapse = ", "),
    "; extra: ", paste(extra, collapse = ", ")
  )
}

run_case <- function(name) {
  start <- proc.time()[["elapsed"]]
  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(script, paste0("--root=", root), paste0("--case=", name)),
    stdout = TRUE,
    stderr = TRUE,
    timeout = timeout
  )
  elapsed <- proc.time()[["elapsed"]] - start
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L

  ok <- identical(as.integer(status), 0L) &&
    any(grepl(paste0("native-smoke-ok: ", name), output, fixed = TRUE))
  list(name = name, ok = ok, status = status, output = output, elapsed = elapsed)
}

if (.Platform$OS.type != "windows" && jobs > 1L) {
  results <- parallel::mclapply(required, run_case, mc.cores = jobs, mc.preschedule = FALSE)
} else {
  results <- lapply(required, run_case)
}

failures <- character()
for (result in results) {
  if (result$ok) {
    cat(sprintf("ok %-40s %7.2fs\n", result$name, result$elapsed))
  } else {
    failures <- c(failures, sprintf(
      "%s (status %d): %s",
      result$name, result$status, paste(result$output, collapse = " | ")
    ))
  }
}

if (length(failures)) {
  cat(paste(failures, collapse = "\n"), "\n", file = stderr())
  quit(save = "no", status = 1L)
}

cat(sprintf(
  "native subprocess smoke: ok (%d entry points, %.2fs total)\n",
  length(required), sum(vapply(results, `[[`, numeric(1L), "elapsed"))
))
