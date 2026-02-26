# Benchmark core clustering methods on shared synthetic datasets.
# Run manually (not during R CMD check) to produce inst/benchmarks/results.csv.
#
# Usage:
#   Rscript inst/benchmarks/bench_methods.R
#
# The script is self-contained and uses only base + package functions;
# it avoids adding new dependencies to the package.

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE, export_all = FALSE, helpers = FALSE)
  } else {
    library(neurocluster)
  }
  library(neuroim2)
})

set.seed(123)

# Thread/parallel fairness control for methods that support parallel execution.
# Use NEUROCLUSTER_BENCH_THREADS=1 to force single-thread runs.
bench_threads <- suppressWarnings(as.integer(Sys.getenv("NEUROCLUSTER_BENCH_THREADS", "0")))
if (is.na(bench_threads)) bench_threads <- 0L
if (bench_threads > 0L) {
  Sys.setenv(OMP_NUM_THREADS = as.character(bench_threads))
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = as.character(bench_threads))
}
bench_parallel <- bench_threads != 1L

choose_k <- function(mask, truth_k = NULL, max_k = 60) {
  if (!is.null(truth_k) && truth_k > 0) return(truth_k)
  k <- round(sum(mask > 0) / 80)
  k <- max(3, min(max_k, k))
  k
}

datasets <- list(
  block_small = function() {
    syn <- make_block_synthetic(dims = c(12, 12, 1), ntime = 80, noise = 0.12, seed = 1)
    syn$truth_k <- 3L
    syn
  },
  blobs_small = function() {
    # Gaussian blobs with orthogonal co-fluctuation signals
    syn <- generate_synthetic_volume(
      scenario = "gaussian_blobs", dims = c(20, 20, 6),
      n_clusters = 4, n_time = 40,
      pattern_type = "orthogonal",  # Use QR-orthogonal signals
      noise_sd = 0.1,
      seed = 2
    )
    syn$truth_k <- syn$n_clusters
    syn
  },
  touching_blobs = function() {
    # Six Gaussian blobs with smaller spread to allow touching at boundaries
    syn <- generate_synthetic_volume(
      scenario = "gaussian_blobs",
      dims = c(32, 32, 8),
      n_clusters = 6,
      n_time = 50,
      pattern_type = "orthogonal",  # QR-orthogonal signals
      noise_sd = 0.1,
      spread = c(4.0, 4.0, 2.5),  # larger spread allows more overlap
      seed = 3
    )
    syn$truth_k <- syn$n_clusters
    syn
  },
  blobs_noisy = function() {
    # Higher noise version to test robustness
    # Each cluster has distinct orthogonal temporal pattern (co-fluctuation)
    syn <- generate_synthetic_volume(
      scenario = "gaussian_blobs", dims = c(20, 20, 6),
      n_clusters = 4, n_time = 50,
      pattern_type = "orthogonal",  # Distinct temporal shapes per cluster
      noise_sd = 0.5,               # Higher noise to test robustness
      seed = 7
    )
    syn$truth_k <- syn$n_clusters
    syn
  },
  blobs_overlap_soft = function() {
    # Overlapping blobs with softer signals and higher noise -> ambiguity
    syn <- generate_synthetic_volume(
      scenario = "gaussian_blobs", dims = c(24, 24, 10),
      n_clusters = 6, n_time = 60,
      spread = c(2.0, 2.0, 2.0),      # tighter blobs -> more overlap
      pattern_type = "random",        # less separable temporal patterns
      noise_sd = 0.25,
      seed = 11
    )
    syn$truth_k <- syn$n_clusters
    syn
  },
  z_layers_soft = function() {
    # Layered slabs with jittered boundaries and harmonic signals
    syn <- generate_synthetic_volume(
      scenario = "z_layers", dims = c(28, 28, 10),
      n_clusters = 5, n_time = 60,
      pattern_type = "harmonic",      # correlated / less orthogonal
      noise_sd = 0.2,
      seed = 12
    )
    syn$truth_k <- syn$n_clusters
    syn
  },
  braided_tubes = function() {
    # Three intertwined curved tubes that approach and cross but never merge
    dims <- c(32, 32, 12)
    n_clusters <- 3
    n_time <- 70
    noise_sd <- 0.18
    radius <- 2.0
    set.seed(21)

    # Grid
    coords <- as.matrix(expand.grid(x = seq_len(dims[1]),
                                    y = seq_len(dims[2]),
                                    z = seq_len(dims[3])))
    nz <- dims[3]
    z_norm <- (coords[, 3] - 1) / (nz - 1)

    # Centerlines for three tubes (phase-shifted sin/cos)
    cx <- cbind(10 + 6 * sin(2 * pi * z_norm),
                16 + 6 * sin(2 * pi * z_norm + 2 * pi / 3),
                22 + 6 * sin(2 * pi * z_norm + 4 * pi / 3))
    cy <- cbind(10 + 6 * cos(2 * pi * z_norm),
                16 + 6 * cos(2 * pi * z_norm + 2 * pi / 3),
                22 + 6 * cos(2 * pi * z_norm + 4 * pi / 3))

    # Assign voxels within radius of any tube centerline
    d2 <- (coords[,1] - cx)^2 + (coords[,2] - cy)^2
    min_d2 <- apply(d2, 1, min)
    labels <- apply(d2, 1, which.min)
    keep <- min_d2 <= radius^2

    # Build mask and truth over kept voxels
    mask_arr <- array(FALSE, dim = dims)
    mask_arr[cbind(coords[keep,1], coords[keep,2], coords[keep,3])] <- TRUE
    mask <- neuroim2::NeuroVol(mask_arr, neuroim2::NeuroSpace(dims, spacing = c(1,1,1)))

    truth <- labels[keep]
    n_vox <- length(truth)

    # Temporal patterns and signals
    patterns <- synthetic_time_patterns(n_clusters, n_time, pattern_types = "random")
    signals <- matrix(0, nrow = n_vox, ncol = n_time)
    for (k in seq_len(n_clusters)) {
      idx <- which(truth == k)
      signals[idx, ] <- matrix(rep(patterns[k, ], each = length(idx)), nrow = length(idx))
    }
    signals <- signals + matrix(rnorm(n_vox * n_time, sd = noise_sd), nrow = n_vox)

    # Build NeuroVec from per-time NeuroVols (avoid index helpers)
    vols <- vector("list", n_time)
    space3d <- neuroim2::NeuroSpace(dims, spacing = c(1,1,1))
    for (t in seq_len(n_time)) {
      vol_t <- array(0, dim = dims)
      vol_t[mask_arr] <- signals[, t]
      vols[[t]] <- neuroim2::NeuroVol(vol_t, space3d)
    }
    vec <- do.call(neuroim2::concat, vols)

    list(
      vec = vec,
      mask = mask,
      truth = truth,
      truth_k = n_clusters,
      n_clusters = n_clusters,
      dims = dims
    )
  }
)

method_grid <- list(
  snic = list(
    list(param_id = 1, compactness = 1.5),
    list(param_id = 2, compactness = 3),
    list(param_id = 3, compactness = 0.8)   # more data-weighted
  ),
  slice_msf = list(
    # Use volumetric (3-D) graph version by always enabling z connectivity.
    list(param_id = 1, r = 12, min_size = 2, compactness = 2.5, gamma = 1.0, nbhd = 8, z_mult = 0.5, stitch_z = TRUE, num_runs = 1, consensus = FALSE),
    list(param_id = 2, r = 8, min_size = 2, compactness = 1.5, gamma = 1.0, nbhd = 8, z_mult = 0.5, stitch_z = TRUE, num_runs = 1, consensus = FALSE),
    list(param_id = 3, r = 12, min_size = 2, compactness = 2.5, gamma = 1.0, nbhd = 4, z_mult = 0.0, stitch_z = TRUE, num_runs = 1, consensus = FALSE)
  ),
  supervoxels = list(
    list(param_id = 1, alpha = 0.3, iterations = 30, sigma1 = 1.0, sigma2 = 1.8, use_gradient = FALSE, connectivity = 6),
    list(param_id = 2, alpha = 0.4, iterations = 30, sigma1 = 1.0, sigma2 = 2.0, use_gradient = FALSE, connectivity = 18)
  ),
  flash3d = list(
    list(param_id = 1, dctM = 8, lambda_t = 1.0, lambda_s = 0.15, bits = 64, rounds = 2),
    list(param_id = 2, dctM = 16, lambda_t = 1.4, lambda_s = 0.35, bits = 64, rounds = 4),
    list(param_id = 3, dctM = 16, lambda_t = 1.4, lambda_s = 0.25, bits = 64, rounds = 4)
  ),
  slic = list(
    list(param_id = 1, spatial_weight = 0.1, connectivity = 26),
    list(param_id = 2, spatial_weight = 0.25, connectivity = 26),
    list(param_id = 3, spatial_weight = 0.05, connectivity = 26) # more data-weighted
  ),
  corr_slic = list(
    list(param_id = 1, tuned = TRUE)
  ),
  brs_slic = list(
    list(param_id = 1, spatial_weight = 0.05, embedding_dim = 24, sketch_repeats = 1, connectivity = 6, max_iterations = 3, boundary_passes = 1, global_passes = 0, refine_spatial_weight = 0.01, refine_l2_weight = 0.02),
    list(param_id = 2, spatial_weight = 0.05, embedding_dim = 24, sketch_repeats = 1, connectivity = 6, max_iterations = 3, boundary_passes = 1, global_passes = 0, refine_spatial_weight = 0.00, refine_l2_weight = 0.08),
    list(param_id = 3, spatial_weight = 0.05, embedding_dim = 24, sketch_repeats = 1, connectivity = 6, max_iterations = 3, boundary_passes = 2, global_passes = 0, refine_spatial_weight = 0.00, refine_l2_weight = 0.08)
  ),
  g3s = list(
    list(param_id = 1, spatial_weight = 0.15),
    list(param_id = 2, spatial_weight = 0.35),
    list(param_id = 3, spatial_weight = 0.1)  # more data-weighted
  ),
  consensus_combo = list(
    list(
      param_id = 1,
      snic_compactness = 1.5,
      slice_r = 8, slice_min_size = 6, slice_compactness = 3, slice_stitch_z = FALSE,
      super_alpha = 0.3, super_sigma1 = 1.0, super_sigma2 = 1.8, super_connectivity = 6,
      flash_dctM = 8, flash_lambda_t = 1.0, flash_lambda_s = 0.6, flash_bits = 64,
      slic_spatial_weight = 0.1, slic_connectivity = 26,
      corrslic_spatial_weight = 0.08, corrslic_embedding_dim = "auto", corrslic_sketch_repeats = 1, corrslic_connectivity = 6, corrslic_assign_stride = 1, corrslic_quantize_assign = TRUE, corrslic_embedding_basis = "dct", corrslic_embedding_whiten = FALSE, corrslic_refine_exact_iters = 1, corrslic_refine_boundary_only = TRUE,
      g3s_spatial_weight = 0.15,
      rena_conn = 6,
      rena_plus_conn = 6,
      acsc_alpha = 0.4
    )
  ),
  rena = list(
    list(param_id = 1, connectivity = 6),
    list(param_id = 2, connectivity = 26)
  ),
  rena_plus = list(
    list(param_id = 1, connectivity = 6),
    list(param_id = 2, connectivity = 26)
  ),
  mcl = list(
    # Speed/accuracy frontier from targeted MCL profiling/tuning.
    list(param_id = 1, inflation = 1.2, connectivity = 6, exact_k = TRUE, max_iterations = 6, prune_k = 32),
    list(param_id = 2, inflation = 1.4, connectivity = 6, exact_k = TRUE, max_iterations = 6, prune_k = 16),
    list(param_id = 3, inflation = 1.6, connectivity = 6, exact_k = TRUE, max_iterations = 6, prune_k = 16)
  ),
  acsc = list(
    list(param_id = 1, alpha = 0.4),
    list(param_id = 2, alpha = 0.6)
  ),
  commute = list(
    list(param_id = 1, spatial_weight = 0.2, ncomp = NA),
    list(param_id = 2, spatial_weight = 0.1, ncomp = NA)
  )
)

corr_slic_family <- function(dataset_name) {
  if (dataset_name %in% c("blobs_small", "block_small")) return("clean_blobs")
  if (dataset_name %in% c("blobs_noisy", "blobs_overlap_soft")) return("noisy_blobs")
  if (dataset_name %in% c("z_layers_soft")) return("layered")
  if (dataset_name %in% c("touching_blobs", "braided_tubes")) return("topology")
  "clean_blobs"
}

corr_slic_tuned <- list(
  clean_blobs = list(
    spatial_weight = 0.05, embedding_dim = 24L, adaptive_embedding = FALSE,
    embedding_basis = "hash", embedding_whiten = FALSE,
    sketch_repeats = 1L,
    quantize_assign = TRUE, assign_stride = 3L,
    max_iterations = 4L,
    refine_exact_iters = 0L, refine_boundary_only = TRUE,
    refine_stride = 1L,
    seed = 7L
  ),
  noisy_blobs = list(
    spatial_weight = 0.08, embedding_dim = 24L, adaptive_embedding = FALSE,
    embedding_basis = "hash", embedding_whiten = FALSE,
    sketch_repeats = 1L,
    quantize_assign = TRUE, assign_stride = 2L,
    max_iterations = 4L,
    refine_exact_iters = 2L, refine_boundary_only = TRUE,
    refine_stride = 1L,
    seed = 7L
  ),
  layered = list(
    spatial_weight = 0.05, embedding_dim = 24L, adaptive_embedding = FALSE,
    embedding_basis = "hash", embedding_whiten = FALSE,
    sketch_repeats = 2L,
    quantize_assign = TRUE, assign_stride = 2L,
    max_iterations = 6L,
    refine_exact_iters = 2L, refine_boundary_only = TRUE,
    refine_stride = 2L,
    seed = 7L
  ),
  topology = list(
    spatial_weight = 0.05, embedding_dim = 24L, adaptive_embedding = FALSE,
    embedding_basis = "hash", embedding_whiten = FALSE,
    sketch_repeats = 1L,
    quantize_assign = TRUE, assign_stride = 2L,
    max_iterations = 4L,
    refine_exact_iters = 0L, refine_boundary_only = TRUE,
    refine_stride = 1L,
    seed = 7L
  )
)

get_corr_slic_cfg <- function(dataset_name) {
  fam <- corr_slic_family(dataset_name)
  cfg <- corr_slic_tuned[[fam]]
  cfg$family <- fam
  cfg
}

run_once <- function(fun, args) {
  gc()
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(do.call(fun, args), error = function(e) e)
  elapsed <- proc.time()[["elapsed"]] - t0
  list(res = res, elapsed = elapsed)
}

# For methods that overshoot K badly on tiny datasets (flash3d), coarsen to K
post_merge_to_k <- function(res, k_target) {
  if (inherits(res, "error")) return(res)
  if (is.null(res$cluster)) return(res)
  labs <- res$cluster
  n_clusters <- length(unique(labs))
  # Only merge if we have MORE clusters than target AND enough centers for kmeans
  if (n_clusters <= k_target || is.null(res$centers)) return(res)
  if (nrow(res$centers) < k_target) return(res)  # Can't merge to more clusters than we have
  # k-means on existing cluster centroids to map to k_target
  km <- kmeans(res$centers, centers = k_target, iter.max = 50, nstart = 3)
  remap <- km$cluster
  res$cluster <- remap[labs]
  res$n_clusters <- length(unique(res$cluster))
  res
}

# Simple adjusted Rand index (no external deps)
adj_rand_index <- function(labels1, labels2) {
  if (length(labels1) != length(labels2)) return(NA_real_)
  tab <- table(labels1, labels2)
  sum_comb <- sum(choose(tab, 2))
  sum_rows <- sum(choose(rowSums(tab), 2))
  sum_cols <- sum(choose(colSums(tab), 2))
  n <- length(labels1)
  expected <- sum_rows * sum_cols / choose(n, 2)
  max_idx <- 0.5 * (sum_rows + sum_cols)
  if (max_idx == expected) return(0)
  (sum_comb - expected) / (max_idx - expected)
}

bench_rows <- list()

for (dname in names(datasets)) {
  cat("Dataset:", dname, "\n")
  ds <- datasets[[dname]]()
  vec <- ds$vec
  mask <- ds$mask
  ntime <- dim(vec)[4]
  nvox <- sum(mask > 0)
  truth_k <- if (!is.null(ds$truth_k)) ds$truth_k else NULL
  base_k <- choose_k(mask, truth_k = truth_k)

  for (m in names(method_grid)) {
    param_list <- method_grid[[m]]
    for (p in param_list) {
      args <- list()
      if (m == "snic") {
        args <- list(vec = vec, mask = mask, compactness = p$compactness, K = base_k)
        fn <- snic
      } else if (m == "slice_msf") {
        args <- list(
          vec = vec, mask = mask,
          r = p$r, min_size = p$min_size, compactness = p$compactness,
          stitch_z = p$stitch_z,
          num_runs = if (!is.null(p$num_runs)) p$num_runs else 1,
          consensus = if (!is.null(p$consensus)) p$consensus else FALSE,
          nbhd = if (!is.null(p$nbhd)) p$nbhd else 8,
          gamma = if (!is.null(p$gamma)) p$gamma else 1.5,
          z_mult = if (!is.null(p$z_mult)) p$z_mult else 0.0,
          target_k_global = base_k
        )
        fn <- slice_msf
      } else if (m == "supervoxels") {
        args <- list(
          bvec = vec, mask = mask, K = base_k,
          alpha = p$alpha,
          iterations = p$iterations,
          sigma1 = p$sigma1,
          sigma2 = p$sigma2,
          use_gradient = p$use_gradient,
          connectivity = p$connectivity
        )
        # Remove NULL arguments
        args <- args[!vapply(args, is.null, logical(1))]
        fn <- supervoxels
      } else if (m == "flash3d") {
        args <- list(
          vec = vec,
          mask = mask,
          n_clusters = base_k,
          method = "flash3d",
          dctM = p$dctM,
          bits = if (!is.null(p$bits)) p$bits else 64,
          spatial_weight = if (!is.null(p$lambda_s)) p$lambda_s else 0.6,
          lambda_t = if (!is.null(p$lambda_t)) p$lambda_t else 1.0,
          max_iterations = if (!is.null(p$rounds)) p$rounds else 2
        )
        fn <- cluster4d
      } else if (m == "slic") {
        args <- list(vec = vec, mask = mask, n_clusters = base_k,
                     method = "slic",
                     spatial_weight = p$spatial_weight,
                     connectivity = p$connectivity,
                     parallel = bench_parallel,
                     preserve_k = TRUE,
                     max_iterations = 12)
        fn <- cluster4d
      } else if (m == "corr_slic") {
        corr_cfg <- get_corr_slic_cfg(dname)
        args <- list(
          vec = vec, mask = mask, n_clusters = base_k,
          method = "corr_slic",
          spatial_weight = corr_cfg$spatial_weight,
          embedding_dim = corr_cfg$embedding_dim,
          adaptive_embedding = isTRUE(corr_cfg$adaptive_embedding),
          embedding_basis = corr_cfg$embedding_basis,
          embedding_whiten = isTRUE(corr_cfg$embedding_whiten),
          sketch_repeats = if (!is.null(corr_cfg$sketch_repeats)) as.integer(corr_cfg$sketch_repeats) else 1L,
          connectivity = 6L,
          max_iterations = if (!is.null(corr_cfg$max_iterations)) as.integer(corr_cfg$max_iterations) else 6L,
          assign_stride = as.integer(corr_cfg$assign_stride),
          quantize_assign = isTRUE(corr_cfg$quantize_assign),
          refine_exact_iters = as.integer(corr_cfg$refine_exact_iters),
          refine_boundary_only = isTRUE(corr_cfg$refine_boundary_only),
          refine_stride = if (!is.null(corr_cfg$refine_stride)) as.integer(corr_cfg$refine_stride) else 1L,
          seed = if (!is.null(corr_cfg$seed)) as.integer(corr_cfg$seed) else 1L,
          parallel = bench_parallel,
          verbose = FALSE
        )
        p <- c(p, corr_cfg)
        fn <- cluster4d
      } else if (m == "brs_slic") {
        args <- list(
          vec = vec, mask = mask, n_clusters = base_k,
          method = "brs_slic",
          spatial_weight = p$spatial_weight,
          embedding_dim = p$embedding_dim,
          sketch_repeats = if (!is.null(p$sketch_repeats)) p$sketch_repeats else 1L,
          connectivity = p$connectivity,
          max_iterations = p$max_iterations,
          boundary_passes = p$boundary_passes,
          global_passes = if (!is.null(p$global_passes)) p$global_passes else 0L,
          refine_spatial_weight = p$refine_spatial_weight,
          refine_l2_weight = if (!is.null(p$refine_l2_weight)) p$refine_l2_weight else 0,
          refine_stride = if (!is.null(p$refine_stride)) p$refine_stride else NULL,
          parallel = bench_parallel,
          verbose = FALSE
        )
        fn <- cluster4d
      } else if (m == "g3s") {
        args <- list(vec = vec, mask = mask, n_clusters = base_k,
                     method = "g3s",
                     spatial_weight = p$spatial_weight,
                     max_iterations = 10)
        fn <- cluster4d
      } else if (m == "consensus_combo") {
        # Run all individual methods then merge via consensus
        t0_local <- proc.time()[["elapsed"]]

        # SNIC
        snic_res <- tryCatch(
          snic(vec = vec, mask = mask, compactness = p$snic_compactness, K = base_k),
          error = function(e) e
        )
        # slice_msf
        slice_res <- tryCatch(
          slice_msf(
            vec = vec, mask = mask,
            r = p$slice_r, min_size = p$slice_min_size, compactness = p$slice_compactness,
            stitch_z = p$slice_stitch_z, num_runs = 1,
            target_k_global = base_k
          ),
          error = function(e) e
        )
        # supervoxels
        super_res <- tryCatch(
          supervoxels(
            bvec = vec, mask = mask, K = base_k,
            alpha = p$super_alpha, iterations = 30,
            sigma1 = p$super_sigma1, sigma2 = p$super_sigma2,
            use_gradient = FALSE, connectivity = p$super_connectivity
          ),
          error = function(e) e
        )
        # flash3d
        flash_res <- tryCatch(
          cluster4d(
            vec = vec, mask = mask, n_clusters = base_k, method = "flash3d",
            dctM = p$flash_dctM, bits = p$flash_bits,
            spatial_weight = p$flash_lambda_s, lambda_t = p$flash_lambda_t,
            max_iterations = 2
          ),
          error = function(e) e
        )
        # slic
        slic_res <- tryCatch(
          cluster4d(
            vec = vec, mask = mask, n_clusters = base_k, method = "slic",
            spatial_weight = p$slic_spatial_weight, connectivity = p$slic_connectivity,
            parallel = bench_parallel,
            preserve_k = TRUE, max_iterations = 12
          ),
          error = function(e) e
        )
        # corr_slic
        corrslic_res <- tryCatch(
          cluster4d(
            vec = vec, mask = mask, n_clusters = base_k, method = "corr_slic",
            spatial_weight = get_corr_slic_cfg(dname)$spatial_weight,
            embedding_dim = get_corr_slic_cfg(dname)$embedding_dim,
            adaptive_embedding = isTRUE(get_corr_slic_cfg(dname)$adaptive_embedding),
            embedding_basis = get_corr_slic_cfg(dname)$embedding_basis,
            embedding_whiten = isTRUE(get_corr_slic_cfg(dname)$embedding_whiten),
            sketch_repeats = if (!is.null(get_corr_slic_cfg(dname)$sketch_repeats)) as.integer(get_corr_slic_cfg(dname)$sketch_repeats) else p$corrslic_sketch_repeats,
            connectivity = p$corrslic_connectivity,
            assign_stride = as.integer(get_corr_slic_cfg(dname)$assign_stride),
            quantize_assign = isTRUE(get_corr_slic_cfg(dname)$quantize_assign),
            refine_exact_iters = as.integer(get_corr_slic_cfg(dname)$refine_exact_iters),
            refine_boundary_only = isTRUE(get_corr_slic_cfg(dname)$refine_boundary_only),
            refine_stride = if (!is.null(get_corr_slic_cfg(dname)$refine_stride)) as.integer(get_corr_slic_cfg(dname)$refine_stride) else 1L,
            seed = if (!is.null(get_corr_slic_cfg(dname)$seed)) as.integer(get_corr_slic_cfg(dname)$seed) else 1L,
            max_iterations = if (!is.null(get_corr_slic_cfg(dname)$max_iterations)) as.integer(get_corr_slic_cfg(dname)$max_iterations) else 6L,
            parallel = bench_parallel,
            verbose = FALSE
          ),
          error = function(e) e
        )
        # g3s
        g3s_res <- tryCatch(
          cluster4d(
            vec = vec, mask = mask, n_clusters = base_k, method = "g3s",
            spatial_weight = p$g3s_spatial_weight,
            max_iterations = 10
          ),
          error = function(e) e
        )
        # rena
        rena_res <- tryCatch(
          rena(bvec = vec, mask = mask, K = base_k, connectivity = p$rena_conn, verbose = FALSE),
          error = function(e) e
        )
        # rena_plus
        renap_res <- tryCatch(
          rena_plus(bvec = vec, mask = mask, K = base_k,
                    connectivity = p$rena_plus_conn,
                    r = 1, lambda = 0,
                    max_iterations = 12, verbose = FALSE),
          error = function(e) e
        )
        # acsc
        acsc_res <- tryCatch(
          acsc(
            bvec = vec,
            mask = mask,
            K = base_k,
            alpha = p$acsc_alpha,
            block_size = 4,
            ann_k = 5,
            refine = TRUE,
            max_refine_iter = 3
          ),
          error = function(e) e
        )

        # If any failed, propagate the error
        sub_errors <- vapply(
          list(snic_res, slice_res, super_res, flash_res, slic_res, corrslic_res, g3s_res, rena_res, renap_res, acsc_res),
          function(obj) inherits(obj, "error"), logical(1)
        )
        if (any(sub_errors)) {
          err_msg <- paste0("component errors: ",
                            paste(c("snic", "slice_msf", "supervoxels", "flash3d", "slic", "corr_slic", "g3s", "rena", "rena_plus", "acsc")[sub_errors],
                                  vapply(list(snic_res, slice_res, super_res, flash_res, slic_res, corrslic_res, g3s_res, rena_res, renap_res, acsc_res)[sub_errors],
                                         conditionMessage, character(1)),
                                  sep = "=", collapse = "; "))
          res <- list(res = simpleError(err_msg), elapsed = NA_real_)
        } else {
          cons_res <- tryCatch(
            merge_clus(snic_res, method = "SE", slice_res, super_res, flash_res, slic_res, corrslic_res, g3s_res, rena_res, renap_res, acsc_res),
            error = function(e) e
          )
          elapsed <- proc.time()[["elapsed"]] - t0_local
          res <- list(res = cons_res, elapsed = elapsed)
        }
        fn <- NULL
      } else if (m == "rena") {
        args <- list(bvec = vec, mask = mask, K = base_k, connectivity = p$connectivity, verbose = FALSE)
        fn <- rena
      } else if (m == "rena_plus") {
        args <- list(bvec = vec, mask = mask, K = base_k,
                     connectivity = p$connectivity,
                     r = 1,            # minimal over-clustering for speed/fairness
                     lambda = 0,       # disable gradient penalty for benchmark
                     max_iterations = 12,
                     verbose = FALSE)
        fn <- rena_plus
      } else if (m == "mcl") {
        args <- list(
          vec = vec,
          mask = mask,
          n_clusters = base_k,
          method = "mcl",
          connectivity = p$connectivity,
          inflation = p$inflation,
          exact_k = p$exact_k,
          max_iterations = p$max_iterations,
          prune_k = if (!is.null(p$prune_k)) p$prune_k else NULL,
          parallel = bench_parallel,
          verbose = FALSE
        )
        fn <- cluster4d
      } else if (m == "acsc") {
        args <- list(
          bvec = vec,
          mask = mask,
          K = base_k,
          alpha = p$alpha,
          block_size = 4,
          ann_k = 5,
          refine = TRUE,
      max_refine_iter = 3
    )
    fn <- acsc
  } else if (m == "commute") {
    args <- list(
      vec = vec,
      mask = mask,
      n_clusters = base_k,
      method = "commute",
      spatial_weight = p$spatial_weight,
      ncomp = if (!is.na(p$ncomp)) p$ncomp else NULL,
      verbose = FALSE
    )
    fn <- cluster4d
  } else {
    next
  }

      if (m == "consensus_combo") {
        # already computed as 'res' above
      } else {
        res <- run_once(fn, args)
        # Flash3d now runs with explicit tuning; skip legacy post-merge heuristic.
      }

      # Standardize consensus ClusteredNeuroVol to list with $cluster
      # FIXED: Check for ClusteredNeuroVol FIRST (before accessing $cluster on S4 object)
      if (!inherits(res$res, "error") && methods::is(res$res, "ClusteredNeuroVol")) {
        # S4 objects don't support $ operator, convert to list using @clusters slot
        res$res <- list(
          cluster = res$res@clusters,
          clusvol = res$res,
          n_clusters = length(unique(res$res@clusters)),
          method = paste0(m, "_consensus")
        )
      }
      # If cluster exists but length mismatches, try to pull from clusvol
      if (!inherits(res$res, "error") && !is.null(res$res$cluster) && length(res$res$cluster) != nvox) {
        if (!is.null(res$res$clusvol)) {
          lab_from_vol <- res$res$clusvol@clusters
          if (length(lab_from_vol) == nvox) {
            res$res$cluster <- lab_from_vol
            res$res$n_clusters <- length(unique(lab_from_vol))
          } else {
            res$res <- simpleError("consensus length mismatch")
          }
        } else {
          res$res <- simpleError("consensus missing clusvol/cluster of correct length")
        }
      }

      n_clusters <- NA_integer_
      ari <- NA_real_
      err <- NA_character_
      if (inherits(res$res, "error")) {
        err <- conditionMessage(res$res)
      } else if (!is.null(res$res$cluster)) {
        n_clusters <- length(unique(res$res$cluster))
        if (n_clusters <= 1) {
          err <- "returned single partition"
        } else if (!is.null(ds$truth) && length(ds$truth) == length(res$res$cluster)) {
          ari <- adj_rand_index(res$res$cluster, ds$truth)
        }
      }

      bench_rows[[length(bench_rows) + 1]] <- data.frame(
        dataset = dname,
        n_vox = nvox,
        n_time = ntime,
        method = m,
        param_id = ifelse("param_id" %in% names(p), p$param_id, NA_integer_),
        params = paste(names(p), p, collapse = ";"),
        elapsed_sec = res$elapsed,
        n_clusters = n_clusters,
        ari = ari,
        error = err,
        stringsAsFactors = FALSE
      )
      cat(sprintf("  %s (%s): %.3fs, clusters=%s\n",
                  m, paste(names(p), p, collapse = ";"),
                  res$elapsed, ifelse(is.na(n_clusters), "NA", n_clusters)))
    }
  }
}

results <- do.call(rbind, bench_rows)
out_dir <- file.path("inst", "benchmarks")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(results, file.path(out_dir, "results.csv"), row.names = FALSE)
cat("\nSaved benchmark results to inst/benchmarks/results.csv\n")
