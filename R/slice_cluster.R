# Seeded DCT subspace specifications for genuine ensemble runs.

.slice_msf_run_specs <- function(n_time, r, num_runs, seed,
                                 ensemble_fraction) {
  r_available <- min(as.integer(r), as.integer(n_time - 1L))
  if (r_available < 1L) stop("slice_msf: at least two timepoints are required")
  if (num_runs == 1L) {
    return(list(list(
      frequencies = seq_len(r_available),
      weights = rep(1, r_available)
    )))
  }

  pool_size <- min(
    as.integer(n_time - 1L),
    max(r_available, as.integer(ceiling(r_available / ensemble_fraction)))
  )
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  lapply(seq_len(num_runs), function(run_id) {
    frequencies <- sort(sample.int(pool_size, r_available, replace = FALSE))
    weights <- exp(stats::rnorm(r_available, sd = 0.35))
    weights <- weights / sqrt(mean(weights^2))
    list(frequencies = frequencies, weights = weights)
  })
}

#' Volumetric DCT minimum-spanning-forest clustering
#'
#' Compresses each voxel time series into a non-DC DCT sketch, weights feature
#' distances by split-half reliability, and applies Felzenszwalb-Huttenlocher
#' segmentation on a masked 3-D neighborhood graph. Multi-run fits use seeded
#' DCT subspace and mode-weight resampling followed by consensus segmentation.
#'
#' @param vec A `NeuroVec` or `SparseNeuroVec` containing the time series.
#' @param mask A `NeuroVol`; exactly finite values greater than zero are included.
#' @param target_k_global `-1` for the natural segmentation or a positive exact
#'   cluster count. Exact K uses the shared adjacency-preserving split/merge engine.
#' @param r Requested number of non-DC DCT modes. It is capped at `T - 1`.
#' @param compactness Finite value in `[0, 10]` controlling FH scale. Larger
#'   values use a smaller scale and generally produce finer components.
#' @param min_size Positive minimum component size for each run.
#' @param num_runs Positive number of runs. Values above one require consensus.
#' @param consensus Whether to fuse a multi-run ensemble.
#' @param stitch_z Whether graph edges may cross axial slices. When `FALSE`,
#'   slices are independent and a global exact-K target is unavailable.
#' @param nbhd `4`, `6`, or `8`. Values `4` and `6` select axial neighbors;
#'   `8` selects the full diagonal neighborhood.
#' @param gamma Non-negative exponent applied to positive split-half reliability;
#'   it directly scales feature edge distances.
#' @param k_fuse Positive FH scale for consensus; defaults to the run scale.
#' @param min_size_fuse Positive minimum consensus component size; defaults to
#'   `min_size`.
#' @param use_features Include sketch similarity as well as label agreement in
#'   multi-run consensus.
#' @param lambda Consensus mixture in `[0, 1]`; one uses label agreement only,
#'   zero uses feature similarity only. Active only when `use_features = TRUE`.
#' @param z_mult DCT-sketch smoothing fraction in `[0, 1]` across axial slices.
#' @param seed Integer seed used to construct multi-run DCT subspaces.
#' @param ensemble_fraction Fraction in `(0, 1]` controlling the candidate DCT
#'   frequency pool for a multi-run ensemble.
#'
#' @return A `slice_msf_cluster_result` and `cluster4d_result`. `cluster` is in
#'   mask order, `centers` is K by T, and `coord_centers` is K by 3. Multi-run
#'   results also contain native `runs` and complete ensemble provenance in
#'   `metadata`.
#'
#' @references Felzenszwalb, P. F., & Huttenlocher, D. P. (2004). Efficient
#'   graph-based image segmentation. IJCV, 59(2), 167-181.
#' @importFrom neuroim2 series spacing
#' @export
slice_msf <- function(vec, mask, 
                      target_k_global = -1,
                      r = 12,
                      compactness = 5,
                      min_size = 80,
                      num_runs = 3,
                      consensus = TRUE,
                      stitch_z = TRUE,
                      nbhd = 8,
                      gamma = 1.5,
                      k_fuse = NULL,
                      min_size_fuse = NULL,
                      use_features = FALSE,
                      lambda = 0.7,
                      z_mult = 0.0,
                      seed = 1L,
                      ensemble_fraction = 0.8) {
  method <- "slice_msf"
  r <- .cluster4d_scalar_number(r, "r", method, lower = 1, integer = TRUE)
  compactness <- .cluster4d_scalar_number(
    compactness, "compactness", method, lower = 0, upper = 10
  )
  min_size <- .cluster4d_scalar_number(
    min_size, "min_size", method, lower = 1, integer = TRUE
  )
  num_runs <- .cluster4d_scalar_number(
    num_runs, "num_runs", method, lower = 1, integer = TRUE
  )
  consensus <- .cluster4d_scalar_logical(consensus, "consensus", method)
  stitch_z <- .cluster4d_scalar_logical(stitch_z, "stitch_z", method)
  nbhd <- .cluster4d_scalar_number(nbhd, "nbhd", method, integer = TRUE)
  if (!nbhd %in% c(4L, 6L, 8L)) stop("slice_msf: nbhd must be 4, 6, or 8")
  if (nbhd == 6L) nbhd <- 4L
  gamma <- .cluster4d_scalar_number(gamma, "gamma", method, lower = 0)
  lambda <- .cluster4d_scalar_number(lambda, "lambda", method, lower = 0, upper = 1)
  use_features <- .cluster4d_scalar_logical(
    use_features, "use_features", method
  )
  z_mult <- .cluster4d_scalar_number(z_mult, "z_mult", method, lower = 0, upper = 1)
  seed <- .cluster4d_scalar_number(
    seed, "seed", method,
    lower = -.Machine$integer.max, upper = .Machine$integer.max, integer = TRUE
  )
  ensemble_fraction <- .cluster4d_scalar_number(
    ensemble_fraction, "ensemble_fraction", method,
    lower = .Machine$double.eps, upper = 1
  )
  if (num_runs > 1L && !consensus) {
    stop("slice_msf: num_runs > 1 requires consensus = TRUE")
  }
  if (num_runs == 1L && use_features) {
    stop("slice_msf: use_features is only active for multi-run consensus")
  }
  target_k_global <- .cluster4d_scalar_number(
    target_k_global, "target_k_global", method, integer = TRUE
  )
  if (target_k_global == 0L || target_k_global < -1L) {
    stop("slice_msf: target_k_global must be -1 or a positive integer")
  }
  if (target_k_global > 0L && !stitch_z) {
    stop("slice_msf: target_k_global requires stitch_z = TRUE")
  }
  
  # Check spatial dimension compatibility
  vec_dims <- dim(vec)[1:3]  # First 3 dimensions are spatial (x, y, z)
  mask_dims <- dim(mask)
  if (!identical(vec_dims, mask_dims)) {
    stop("NeuroVec and mask must have identical spatial dimensions. ",
         "NeuroVec dimensions: ", paste(vec_dims, collapse="x"), 
         ", mask dimensions: ", paste(mask_dims, collapse="x"))
  }
  
  effective_k <- if (target_k_global > 0L) target_k_global else 1L
  input <- validate_cluster4d_inputs(vec, mask, effective_k, method)
  mask.idx <- input$mask_idx

  # For very small test problems, avoid spinning up many RcppParallel threads.
  # This reduces overhead and prevents rare hangs on some platforms/CI.
  if (length(mask.idx) < 2000L && requireNamespace("RcppParallel", quietly = TRUE)) {
    old_opts <- RcppParallel::setThreadOptions(numThreads = 1L)
    on.exit({
      tryCatch({
        if (is.list(old_opts) && !is.null(old_opts$numThreads)) {
          if (!is.null(old_opts$stackSize)) {
            RcppParallel::setThreadOptions(numThreads = old_opts$numThreads, stackSize = old_opts$stackSize)
          } else {
            RcppParallel::setThreadOptions(numThreads = old_opts$numThreads)
          }
        } else if (is.null(old_opts)) {
          tryCatch(RcppParallel::setThreadOptions(numThreads = "auto"), error = function(e) NULL)
        } else if (is.character(old_opts)) {
          tryCatch(RcppParallel::setThreadOptions(numThreads = old_opts), error = function(e) NULL)
        } else {
          RcppParallel::setThreadOptions(numThreads = as.integer(old_opts)[1])
        }
      }, error = function(e) NULL)
    }, add = TRUE)
  }
  # Get time series matrix for all voxels
  # The C++ function expects data for the full volume
  all_idx <- seq_len(prod(dim(mask)))
  # series() returns T x N when vec is a NeuroVec
  TS <- series(vec, all_idx)
  if (!is.matrix(TS)) {
    TS <- matrix(TS, nrow = dim(vec)[4])  # Ensure 1 x N for single-volume inputs
  }
  if (nrow(TS) < 2L) stop("slice_msf: at least two timepoints are required")
  
  # Prepare mask and dimensions
  included_mask <- cluster4d_mask_array(mask, method)
  mask_flat <- as.integer(as.vector(included_mask))
  vol_dim <- dim(mask)
  
  # Map compactness to fh_scale parameter (inverse relationship)
  # Higher compactness -> smaller fh_scale -> more compact clusters
  fh_scale <- 2.0 / (compactness + 1.0)  # Maps compactness 5 -> fh_scale ≈ 0.33
  
  # Get voxel dimensions
  voxel_dim <- spacing(mask)
  
  # Set fusion parameters
  if (is.null(k_fuse)) k_fuse <- fh_scale
  if (is.null(min_size_fuse)) min_size_fuse <- min_size
  k_fuse <- .cluster4d_scalar_number(
    k_fuse, "k_fuse", method, lower = .Machine$double.eps
  )
  min_size_fuse <- .cluster4d_scalar_number(
    min_size_fuse, "min_size_fuse", method, lower = 1, integer = TRUE
  )

  run_specs <- .slice_msf_run_specs(
    nrow(TS), r, num_runs, seed, ensemble_fraction
  )
  run_one <- function(run_id) {
    spec <- run_specs[[run_id]]
    slice_msf_runwise(
      TS = TS,
      mask = mask_flat,
      vol_dim = vol_dim,
      r = r,
      fh_scale = fh_scale,
      min_size = min_size,
      nbhd = nbhd,
      stitch_z = stitch_z,
      rows_are_time = TRUE,
      gamma = gamma,
      voxel_dim = voxel_dim,
      spatial_beta = 0.0,
      target_k_global = -1L,
      target_k_per_slice = -1L,
      z_mult = z_mult,
      w_threshold = 0,
      dct_frequencies = as.integer(spec$frequencies),
      dct_weights = as.numeric(spec$weights)
    )
  }

  if (num_runs == 1L) {
    result <- run_one(1L)
    runs <- NULL
  } else {
    runs <- lapply(seq_len(num_runs), run_one)
    result <- slice_fuse_consensus(
      run_results = runs,
      vol_dim = vol_dim,
      nbhd = nbhd,
      fh_scale = k_fuse,
      min_size = min_size_fuse,
      use_features = use_features,
      lambda = lambda,
      voxel_dim = voxel_dim,
      spatial_beta = 0.0,
      target_k_global = -1L,
      target_k_per_slice = -1L,
      stitch_z = stitch_z
    )
  }

  labels <- as.integer(result$labels)
  # Extract labels for mask voxels only
  cluster_ids <- labels[mask.idx]
  if (any(cluster_ids <= 0L)) {
    stop("slice_msf: native result omitted one or more included mask voxels")
  }
  original <- .cluster4d_original_data(vec, mask, method)
  graph_info <- .exact_k_graph(mask, if (nbhd == 4L) 6L else 26L)
  cluster_ids <- .exact_k_connected_labels(
    cluster_ids, graph_info$graph, graph_info$edges
  )
  if (target_k_global > 0L) {
    cluster_ids <- force_exact_k(
      cluster_ids, original$features, target_k_global,
      mask = mask, connectivity = if (nbhd == 4L) 6L else 26L
    )
  }
  labels[] <- 0L
  labels[mask.idx] <- cluster_ids

  run_metadata <- lapply(seq_len(num_runs), function(i) {
    list(
      run_id = as.integer(i),
      dct_frequencies = as.integer(run_specs[[i]]$frequencies),
      dct_weights = as.numeric(run_specs[[i]]$weights),
      native_params = if (is.null(runs)) result$params else runs[[i]]$params,
      n_clusters = as.integer(length(unique(
        (if (is.null(runs)) result$labels else runs[[i]]$labels)[mask.idx]
      )))
    )
  })
  parameters <- list(
    target_k_global = if (target_k_global > 0L) target_k_global else NULL,
    r = r,
    compactness = compactness,
    fh_scale = fh_scale,
    min_size = min_size,
    num_runs = num_runs,
    consensus = consensus,
    stitch_z = stitch_z,
    nbhd = nbhd,
    gamma = gamma,
    k_fuse = k_fuse,
    min_size_fuse = min_size_fuse,
    use_features = use_features,
    lambda = lambda,
    z_mult = z_mult,
    seed = seed,
    ensemble_fraction = ensemble_fraction
  )
  ret <- structure(
    list(
      cluster = cluster_ids,
      cluster_map = array(labels, dim = vol_dim),
      method = "slice_msf",
      parameters = parameters,
      metadata = list(
        native = result$params,
        ensemble = list(
          seeded = TRUE,
          seed = seed,
          subspace_fraction = ensemble_fraction,
          runs = run_metadata
        ),
        consensus = if (num_runs > 1L) result$params else NULL,
        topology = list(
          volumetric = stitch_z,
          connectivity = if (nbhd == 4L) 6L else 26L,
          exact_k_engine = if (target_k_global > 0L) "shared_adjacency_preserving" else NULL
        )
      )
    ),
    class = c(
      "slice_msf_cluster_result", "cluster4d_result", "cluster_result", "list"
    )
  )
  if (!is.null(runs)) {
    ret$runs <- runs
  }
  finalize_cluster4d_result(ret, vec, mask, "slice_msf", parameters)
}

#' Run one Slice-MSF segmentation
#'
#' Low-level diagnostic interface returning native full-volume labels, split-half
#' reliability weights, and an `r_used` by voxel DCT sketch.
#'
#' @inheritParams slice_msf
#' @param k Positive FH segmentation scale.
#' @param dct_frequencies Optional unique non-DC DCT frequencies in `[1, T - 1]`.
#' @param dct_weights Optional positive weights paired with `dct_frequencies`.
#'
#' @return A list with full-volume `labels`, `weights`, `sketch`, and `params`.
#' @export
slice_msf_single <- function(vec, mask, 
                            r = 12,
                            k = 0.32,
                            min_size = 80,
                            nbhd = 8,
                            stitch_z = TRUE,
                            gamma = 1.5,
                            z_mult = 0.0,
                            dct_frequencies = NULL,
                            dct_weights = NULL) {
  method <- "slice_msf_single"
  validate_cluster4d_inputs(vec, mask, 1L, method)
  r <- .cluster4d_scalar_number(r, "r", method, lower = 1, integer = TRUE)
  k <- .cluster4d_scalar_number(
    k, "k", method, lower = .Machine$double.eps
  )
  min_size <- .cluster4d_scalar_number(
    min_size, "min_size", method, lower = 1, integer = TRUE
  )
  nbhd <- .cluster4d_scalar_number(nbhd, "nbhd", method, integer = TRUE)
  if (!nbhd %in% c(4L, 6L, 8L)) {
    stop("slice_msf_single: nbhd must be 4, 6, or 8", call. = FALSE)
  }
  if (nbhd == 6L) nbhd <- 4L
  stitch_z <- .cluster4d_scalar_logical(stitch_z, "stitch_z", method)
  gamma <- .cluster4d_scalar_number(gamma, "gamma", method, lower = 0)
  z_mult <- .cluster4d_scalar_number(
    z_mult, "z_mult", method, lower = 0, upper = 1
  )

  all_idx <- seq_len(prod(dim(mask)))
  TS <- series(vec, all_idx)
  if (!is.matrix(TS)) {
    TS <- matrix(TS, nrow = dim(vec)[4])
  }
  if (nrow(TS) < 2L) {
    stop("slice_msf_single: at least two timepoints are required", call. = FALSE)
  }

  if (is.null(dct_frequencies)) {
    dct_frequencies <- seq_len(min(r, nrow(TS) - 1L))
  }
  if (!is.numeric(dct_frequencies) || length(dct_frequencies) < 1L ||
      anyNA(dct_frequencies) || any(!is.finite(dct_frequencies)) ||
      any(dct_frequencies != as.integer(dct_frequencies)) ||
      any(dct_frequencies < 1L | dct_frequencies >= nrow(TS)) ||
      anyDuplicated(dct_frequencies)) {
    stop(
      "slice_msf_single: dct_frequencies must be unique integers in [1, T - 1]",
      call. = FALSE
    )
  }
  dct_frequencies <- as.integer(dct_frequencies)
  if (is.null(dct_weights)) dct_weights <- rep(1, length(dct_frequencies))
  if (!is.numeric(dct_weights) || length(dct_weights) != length(dct_frequencies) ||
      anyNA(dct_weights) || any(!is.finite(dct_weights)) || any(dct_weights <= 0)) {
    stop(
      "slice_msf_single: dct_weights must be positive finite values matching dct_frequencies",
      call. = FALSE
    )
  }

  mask_flat <- as.integer(as.vector(cluster4d_mask_array(mask, method)))
  vol_dim <- dim(mask)
  voxel_dim <- spacing(mask)

  slice_msf_runwise(
    TS = TS,
    mask = mask_flat,
    vol_dim = vol_dim,
    r = r,
    fh_scale = k,
    min_size = min_size,
    nbhd = nbhd,
    stitch_z = stitch_z,
    rows_are_time = TRUE,
    gamma = gamma,
    voxel_dim = voxel_dim,
    spatial_beta = 0.0,
    target_k_global = -1L,
    target_k_per_slice = -1L,
    z_mult = z_mult,
    w_threshold = 0,
    dct_frequencies = dct_frequencies,
    dct_weights = as.numeric(dct_weights)
  )
}

#' Fuse Slice-MSF runs on a masked graph
#'
#' Builds neighbor-edge co-association scores from at least two run results and
#' optionally mixes those scores with DCT-sketch similarity.
#'
#' @param run_results At least two results from [slice_msf_single()].
#' @param mask The original `NeuroVol` mask. Outside-mask labels are discarded.
#' @param nbhd `4`, `6`, or `8`; see [slice_msf()].
#' @param k_fuse Positive FH scale for the consensus graph.
#' @param min_size_fuse Positive minimum consensus component size.
#' @param use_features Whether DCT-sketch similarity contributes to edges.
#' @param lambda Mixture in `[0, 1]`; active when `use_features = TRUE`.
#' @param stitch_z Whether consensus edges may cross axial slices.
#'
#' @return A list with full-volume consensus `labels` and complete `params`.
#' @export
slice_msf_consensus <- function(run_results, mask,
                               nbhd = 8,
                               k_fuse = 0.30,
                               min_size_fuse = 80,
                               use_features = FALSE,
                               lambda = 0.7,
                               stitch_z = TRUE) {
  method <- "slice_msf_consensus"
  if (!inherits(mask, "NeuroVol")) {
    stop("slice_msf_consensus: mask must be a NeuroVol", call. = FALSE)
  }
  if (!is.list(run_results) || length(run_results) < 2L) {
    stop("slice_msf_consensus: run_results must contain at least two runs", call. = FALSE)
  }
  nbhd <- .cluster4d_scalar_number(nbhd, "nbhd", method, integer = TRUE)
  if (!nbhd %in% c(4L, 6L, 8L)) {
    stop("slice_msf_consensus: nbhd must be 4, 6, or 8", call. = FALSE)
  }
  if (nbhd == 6L) nbhd <- 4L
  k_fuse <- .cluster4d_scalar_number(
    k_fuse, "k_fuse", method, lower = .Machine$double.eps
  )
  min_size_fuse <- .cluster4d_scalar_number(
    min_size_fuse, "min_size_fuse", method, lower = 1, integer = TRUE
  )
  use_features <- .cluster4d_scalar_logical(
    use_features, "use_features", method
  )
  lambda <- .cluster4d_scalar_number(
    lambda, "lambda", method, lower = 0, upper = 1
  )
  stitch_z <- .cluster4d_scalar_logical(stitch_z, "stitch_z", method)

  vol_dim <- dim(mask)
  voxel_dim <- spacing(mask)
  included <- as.vector(cluster4d_mask_array(mask, method))
  n_voxels <- prod(vol_dim)
  run_results <- lapply(run_results, function(run) {
    if (!is.list(run) || is.null(run$labels) || length(run$labels) != n_voxels) {
      stop(
        "slice_msf_consensus: every run must have one label per mask voxel",
        call. = FALSE
      )
    }
    run$labels <- as.integer(run$labels)
    if (any(run$labels[included] <= 0L)) {
      stop(
        "slice_msf_consensus: every run must label every included mask voxel",
        call. = FALSE
      )
    }
    run$labels[!included] <- 0L
    if (!is.null(run$weights)) run$weights[!included] <- 0
    if (!is.null(run$sketch)) run$sketch[, !included] <- 0
    run
  })

  slice_fuse_consensus(
    run_results = run_results,
    vol_dim = vol_dim,
    nbhd = nbhd,
    fh_scale = k_fuse,
    min_size = min_size_fuse,
    use_features = use_features,
    lambda = lambda,
    voxel_dim = voxel_dim,
    spatial_beta = 0.0,
    target_k_global = -1L,
    target_k_per_slice = -1L,
    stitch_z = stitch_z
  )
}

# Helper function to compute cluster centers
compute_slice_cluster_centers <- function(vec, mask, cluster_ids, mask.idx) {
  unique_clusters <- sort(unique(cluster_ids[cluster_ids > 0]))
  n_clusters <- length(unique_clusters)
  
  if (n_clusters == 0) {
    return(list(centers = NULL, coord_centers = NULL))
  }
  
  # Get time series data
  vecmat <- series(vec, mask.idx)
  if (!is.matrix(vecmat)) {
    vecmat <- matrix(vecmat, nrow = dim(vec)[4])  # Keep time in rows for single-volume inputs
  }
  
  # Compute data centers (mean time series per cluster)
  centers <- matrix(0, nrow = n_clusters, ncol = nrow(vecmat))
  coord_centers <- matrix(0, nrow = n_clusters, ncol = 3)
  
  coords <- index_to_coord(mask, mask.idx)
  
  for (i in seq_along(unique_clusters)) {
    clust_id <- unique_clusters[i]
    idx <- which(cluster_ids == clust_id)
    
    if (length(idx) > 0) {
      # Mean time series
      centers[i, ] <- rowMeans(vecmat[, idx, drop = FALSE])
      
      # Spatial centroid
      coord_centers[i, ] <- colMeans(coords[idx, , drop = FALSE])
    }
  }
  
  list(centers = centers, coord_centers = coord_centers)
}
