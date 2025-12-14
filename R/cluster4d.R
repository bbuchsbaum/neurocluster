#' Unified 4D Clustering for Neuroimaging Data
#'
#' Performs spatially-constrained clustering on 4D neuroimaging data using
#' various algorithms. This is the main entry point for all clustering methods
#' in the neurocluster package.
#'
#' @param vec A \code{NeuroVec} instance supplying the 4D data (x, y, z, time) to cluster
#' @param mask A \code{NeuroVol} mask defining the voxels to include in clustering.
#'   If numeric, nonzero values define included voxels. If logical, TRUE values
#'   define included voxels.
#' @param n_clusters Target number of clusters (default 100). Note that some methods
#'   may produce slightly different numbers of clusters due to algorithmic constraints.
#' @param method Clustering algorithm to use. Options:
#'   \itemize{
#'     \item \code{"supervoxels"}: Iterative heat kernel-based clustering (default)
#'     \item \code{"snic"}: Simple Non-Iterative Clustering
#'     \item \code{"slic"}: SLIC superpixels extended to 4D
#'     \item \code{"slice_msf"}: Slice-wise Minimum Spanning Forest (fast but may show z-artifacts)
#'     \item \code{"flash3d"}: Fast Low-rank Approximate Superclusters
#'     \item \code{"g3s"}: Gradient-Guided Geodesic Supervoxels (NEW - recommended for best quality/speed)
#'     \item \code{"rena"}: Recursive Nearest Agglomeration (fast, balanced, topology-aware)
#'     \item \code{"acsc"}: Adaptive Correlation Superclustering (graph-based with boundary refinement)
#'   }
#' @param spatial_weight Balance between spatial and feature similarity (0-1).
#'   Higher values emphasize spatial compactness. Default 0.5.
#'   Maps to method-specific parameters:
#'   \itemize{
#'     \item supervoxels: \code{alpha = 1 - spatial_weight} (0 = all spatial, 1 = all feature)
#'     \item snic/slic: \code{compactness = spatial_weight * 20} (typical range 1-20)
#'     \item slice_msf: \code{compactness = spatial_weight * 10} (typical range 1-10)
#'     \item flash3d: \code{lambda_s = spatial_weight} (direct mapping)
#'   }
#' @param max_iterations Maximum iterations for iterative methods. Default 10.
#'   Maps to: \code{iterations} (supervoxels), \code{max_iter} (snic/slic),
#'   \code{rounds} (flash3d).
#' @param connectivity Neighborhood connectivity (6 or 26). Default 26.
#'   6 = face neighbors only, 26 = face + edge + corner neighbors.
#' @param parallel Enable parallel processing where supported. Default TRUE.
#' @param verbose Print progress messages. Default FALSE.
#' @param ... Additional method-specific parameters. See method documentation for details.
#'
#' @return A \code{cluster4d_result} object (also inherits from \code{cluster_result})
#'   containing:
#'   \item{clusvol}{A \code{ClusteredNeuroVol} with cluster assignments}
#'   \item{cluster}{Integer vector of cluster assignments for each masked voxel}
#'   \item{centers}{Matrix of cluster centers in feature space (n_clusters x timepoints)}
#'   \item{coord_centers}{Matrix of cluster spatial centers (n_clusters x 3)}
#'   \item{n_clusters}{Actual number of clusters produced}
#'   \item{method}{Clustering method used}
#'   \item{parameters}{List of all parameters used}
#'   \item{metadata}{Method-specific additional information}
#'
#' @section Algorithm Comparison:
#'
#' \tabular{lllllll}{
#'   \strong{Method} \tab \strong{Speed} \tab \strong{3D Continuity} \tab \strong{Memory} \tab \strong{Parallel} \tab \strong{Best For} \cr
#'   supervoxels \tab Slow \tab Excellent \tab High \tab Yes \tab Small-medium data, smooth parcels \cr
#'   snic \tab Fast \tab Good \tab Low \tab No \tab Large data, non-iterative \cr
#'   slic \tab Fast \tab Good \tab Medium \tab Yes \tab Balanced speed/quality \cr
#'   slice_msf \tab Very Fast \tab Moderate \tab Low \tab Yes \tab High-res data, accept z-artifacts \cr
#'   flash3d \tab Fast \tab Good \tab Medium \tab Partial \tab Large data, hash-based \cr
#'   rena \tab Fast \tab Excellent \tab Low \tab No \tab Balanced clusters, topology-aware \cr
#' }
#'
#' @section Parameter Guidelines:
#' 
#' \strong{For whole-brain parcellation:}
#' \itemize{
#'   \item n_clusters: 100-1000 depending on desired granularity
#'   \item spatial_weight: 0.4-0.6 for balanced clustering
#'   \item connectivity: 26 for smoother boundaries
#' }
#' 
#' \strong{For ROI analysis:}
#' \itemize{
#'   \item n_clusters: 10-100 depending on ROI size
#'   \item spatial_weight: 0.2-0.4 to emphasize functional similarity
#'   \item connectivity: 6 for more discrete parcels
#' }
#' 
#' \strong{For high-resolution data (< 2mm):}
#' \itemize{
#'   \item method: "slice_msf" or "flash3d" for speed
#'   \item n_clusters: Scale with voxel count (roughly n_voxels/200)
#' }
#'
#' @examples
#' \dontrun{
#' # Simple synthetic example (runs quickly for testing)
#' library(neuroim2)
#' mask <- NeuroVol(array(1, c(4,4,4)), NeuroSpace(c(4,4,4)))
#' vec <- NeuroVec(array(rnorm(4*4*4*10), c(4,4,4,10)), 
#'                 NeuroSpace(c(4,4,4,10)))
#' result <- cluster4d(vec, mask, n_clusters = 3, method = "g3s", 
#'                    max_iterations = 1)
#' print(result$n_clusters)
#' }
#'
#' \dontrun{
#' # More realistic examples with larger data
#' mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
#' vec <- replicate(50, NeuroVol(array(runif(20*20*20), c(20,20,20)),
#'                               NeuroSpace(c(20,20,20))), simplify=FALSE)
#' vec <- do.call(concat, vec)
#' 
#' # Basic usage with default supervoxels method
#' result <- cluster4d(vec, mask, n_clusters = 100)
#' 
#' # Fast clustering with FLASH-3D (hash-based)
#' result <- cluster4d(vec, mask, n_clusters = 100, method = "flash3d")
#' 
#' # Emphasize spatial compactness
#' result <- cluster4d(vec, mask, n_clusters = 100, spatial_weight = 0.8)
#' 
#' # Use specific method with custom parameters
#' result <- cluster4d(vec, mask, n_clusters = 100, 
#'                    method = "slice_msf",
#'                    num_runs = 3,  # slice_msf-specific parameter
#'                    consensus = TRUE)
#' 
#' # Get parameter suggestions for your data
#' n_vox <- sum(mask > 0)
#' n_time <- dim(vec)[4]
#' params <- suggest_cluster4d_params(n_vox, n_time, priority = "quality")
#' result <- cluster4d(vec, mask, 
#'                    n_clusters = params$n_clusters,
#'                    method = params$recommended_method)
#' }
#' @seealso 
#' Method-specific functions: \code{\link{cluster4d_supervoxels}}, 
#' \code{\link{cluster4d_snic}}, \code{\link{cluster4d_slic}},
#' \code{\link{cluster4d_slice_msf}}, \code{\link{cluster4d_flash3d}},
#' \code{\link{cluster4d_commute}}
#' 
#' Legacy functions (deprecated): \code{\link{supervoxels}}, \code{\link{snic}},
#' \code{\link{slic4d_supervoxels}}, \code{\link{slice_msf}}, \code{\link{supervoxels_flash3d}}
#'
#' @export
#' @importFrom neuroim2 NeuroVec NeuroVol ClusteredNeuroVol series index_to_coord spacing
cluster4d <- function(vec, mask,
                     n_clusters = 100,
                     method = c("supervoxels", "snic", "slic", "slice_msf", "flash3d", "g3s", "rena", "rena_plus", "acsc", "commute"),
                     spatial_weight = 0.5,
                     max_iterations = NULL,
                     connectivity = NULL,
                     parallel = TRUE,
                     verbose = FALSE,
                     ...) {

  # Allow users to pass a single-volume NeuroVol; wrap to NeuroVec for downstream code
  vec <- ensure_neurovec(vec)

  method <- match.arg(method)

  # ------------------------------
  # Method-specific sane defaults
  # ------------------------------
  if (is.null(max_iterations)) {
    max_iterations <- switch(method,
      supervoxels = 50,   # algorithm expects 30–50 iters for stability
    flash3d    = 2,     # native FLASH rounds; keeps tests fair
    snic       = 100,
    slic       = 10,
    slice_msf  = 10,
    g3s        = 5,     # refinement iterations
    rena       = 12,
    rena_plus  = 12,
    acsc       = 5,
    commute    = 1,     # single pass (no iterative refinement)
    10
  )
  }

  if (is.null(connectivity)) {
    connectivity <- switch(method,
      supervoxels = 27,
      snic        = 26,
      slic        = 26,
      slice_msf   = 26,
      rena        = 26,
      rena_plus   = 26,
      acsc        = 26,
      g3s         = 26,
      flash3d     = NA_integer_,
      26
    )
  }

  # Validate common inputs
  validate_cluster4d_inputs(vec, mask, n_clusters, paste0("cluster4d:", method))

  # Validate spatial_weight
  if (!is.numeric(spatial_weight) || spatial_weight < 0 || spatial_weight > 1) {
    stop("spatial_weight must be between 0 and 1")
  }

  # Validate connectivity (only when method uses it)
  if (!is.na(connectivity)) {
    if (!connectivity %in% c(6, 18, 26, 27)) {
      stop("connectivity must be 6, 18, 26, or 27")
    }
  }

  # Store original parameters
  orig_params <- list(
    n_clusters = n_clusters,
    method = method,
    spatial_weight = spatial_weight,
    max_iterations = max_iterations,
    connectivity = connectivity,
    parallel = parallel,
    verbose = verbose
  )

  # Dispatch to method-specific implementation
  result <- switch(method,
    supervoxels = cluster4d_supervoxels(vec, mask, n_clusters, spatial_weight,
                                        max_iterations, connectivity, parallel,
                                        verbose, ...),
    snic = cluster4d_snic(vec, mask, n_clusters, spatial_weight,
                         max_iterations, verbose, ...),
    slic = cluster4d_slic(vec, mask, n_clusters, spatial_weight,
                         max_iterations, connectivity, parallel, verbose, ...),
    slice_msf = cluster4d_slice_msf(vec, mask, n_clusters, spatial_weight,
                                   connectivity, parallel, verbose, ...),
    flash3d = cluster4d_flash3d(vec, mask, n_clusters, spatial_weight,
                               max_iterations, verbose, ...),
    g3s = cluster4d_g3s(vec, mask, n_clusters,
                       alpha = 1 - spatial_weight,  # Convert to feature weight
                       max_refinement_iter = max_iterations,
                       verbose = verbose, ...),
    rena = cluster4d_rena(vec, mask, n_clusters, spatial_weight,
                         max_iterations, connectivity, verbose, ...),
    rena_plus = cluster4d_rena_plus(vec, mask, n_clusters, spatial_weight,
                                    connectivity = connectivity,
                                    max_iterations = max_iterations,
                                    verbose = verbose, ...),
    acsc = cluster4d_acsc(vec, mask, n_clusters, spatial_weight,
                          max_iterations, verbose, ...),
    commute = cluster4d_commute(vec, mask, n_clusters, spatial_weight,
                                verbose, ...)
  )

  # Ensure result has cluster4d_result class
  if (!"cluster4d_result" %in% class(result)) {
    class(result) <- c("cluster4d_result", class(result))
  }

  # Add original parameters if not present
  if (is.null(result$parameters)) {
    result$parameters <- orig_params
  }

  result
}

#' Cluster4d using commute-time spectral method
#'
#' Wrapper for commute_cluster with standardized interface.
#'
#' @inheritParams cluster4d
#' @param ncomp Number of embedding components (defaults to sqrt(2K)).
#' @return A cluster4d_result object
#' @export
#' @note Commute-time embedding is O(N^3) eigen; use only on small ROIs.
cluster4d_commute <- function(vec, mask, n_clusters = 100,
                             spatial_weight = 0.5,
                             verbose = FALSE,
                             ncomp = NULL,
                             ...) {

  # Map spatial_weight (0=all spatial, 1=all feature) to commute alpha (feature weight)
  alpha <- 1 - spatial_weight
  if (is.null(ncomp)) {
    base <- ceiling(sqrt(n_clusters * 2))
    n_vox <- sum(mask > 0)
    # For small problems (<=2000 vox), allow richer embedding to capture structure
    ncomp <- if (n_vox <= 2000) max(base, min(n_clusters * 2, n_vox - 1)) else base
  }

  result <- commute_cluster(
    bvec = vec,
    mask = mask,
    K = n_clusters,
    ncomp = ncomp,
    alpha = alpha,
    verbose = verbose,
    ...
  )

  if (!"cluster4d_result" %in% class(result)) {
    class(result) <- c("cluster4d_result", class(result))
  }

  result$method <- "commute"
  result$n_clusters <- length(unique(result$cluster[!is.na(result$cluster)]))
  dots <- list(...)
  result$parameters <- c(
    list(
      n_clusters_requested = n_clusters,
      spatial_weight = spatial_weight,
      alpha = alpha,
      ncomp = ncomp
    ),
    dots
  )

  result
}

#' Cluster4d using supervoxels method
#'
#' Wrapper for supervoxels algorithm with standardized interface.
#'
#' @inheritParams cluster4d
#' @param sigma1 Bandwidth of heat kernel for features (supervoxels-specific)
#' @param sigma2 Bandwidth of heat kernel for coordinates (supervoxels-specific)
#' @param use_gradient Use gradient-based initialization
#' @param converge_thresh Convergence threshold
#' @param ... Additional parameters passed to supervoxels
#'
#' @return A cluster4d_result object
#' @export
cluster4d_supervoxels <- function(vec, mask, n_clusters = 100,
                                 spatial_weight = 0.5,
                                 max_iterations = 50,
                                 connectivity = 27,
                                 parallel = TRUE,
                                 verbose = FALSE,
                                 sigma1 = 1, sigma2 = 2.5,
                                 use_gradient = TRUE,
                                 converge_thresh = 0.001,
                                 ...) {
  
  # Convert spatial_weight to alpha (feature weight)
  # spatial_weight = 0.5 -> alpha = 0.5 (balanced)
  # spatial_weight = 1 -> alpha = 0 (all spatial)
  # spatial_weight = 0 -> alpha = 1 (all feature)
  alpha <- 1 - spatial_weight
  
  # Call original supervoxels with mapped parameters
  result <- supervoxels(
    bvec = vec,
    mask = mask,
    K = n_clusters,
    alpha = alpha,
    iterations = max_iterations,
    connectivity = connectivity,
    parallel = parallel,
    verbose = verbose,
    sigma1 = sigma1,
    sigma2 = sigma2,
    use_gradient = use_gradient,
    converge_thresh = converge_thresh,
    ...
  )
  
  # Standardize result structure
  if (!"cluster4d_result" %in% class(result)) {
    class(result) <- c("cluster4d_result", class(result))
  }
  
  # Add method info
  result$method <- "supervoxels"
  result$n_clusters <- length(unique(result$cluster[!is.na(result$cluster)]))
  
  # Store all parameters including those passed through ...
  dots <- list(...)
  result$parameters <- c(
    list(
      n_clusters_requested = n_clusters,
      spatial_weight = spatial_weight,
      alpha = alpha,
      max_iterations = max_iterations,
      connectivity = connectivity,
      parallel = parallel,
      sigma1 = sigma1,
      sigma2 = sigma2,
      use_gradient = use_gradient,
      converge_thresh = converge_thresh
    ),
    dots  # Include any additional parameters passed through
  )
  
  result
}

#' Cluster4d using SNIC method
#'
#' Wrapper for SNIC algorithm with standardized interface.
#'
#' @inheritParams cluster4d
#' @param ... Additional parameters passed to snic
#'
#' @return A cluster4d_result object
#' @export
cluster4d_snic <- function(vec, mask, n_clusters = 100,
                          spatial_weight = 0.5,
                          max_iterations = 100,
                          verbose = FALSE,
                          ...) {
  
  # Convert spatial_weight to compactness
  # SNIC uses smaller compactness values (1-10 typical)
  compactness <- spatial_weight * 10
  
  # Call original snic
  result <- snic(
    vec = vec,
    mask = mask,
    K = n_clusters,
    compactness = compactness,
    max_iter = max_iterations,
    ...
  )
  
  # Standardize result structure
  if (!"cluster4d_result" %in% class(result)) {
    class(result) <- c("cluster4d_result", class(result))
  }
  
  # Add method info
  result$method <- "snic"
  result$n_clusters <- length(unique(result$cluster[!is.na(result$cluster)]))
  
  # Compute centers if missing
  if (is.null(result$centers) || is.null(result$coord_centers)) {
    data_prep <- prepare_cluster4d_data(vec, mask)
    center_info <- compute_cluster_centers(result$cluster, data_prep$features, data_prep$coords)
    result$centers <- center_info$centers
    result$coord_centers <- center_info$coord_centers
  }
  
  # Store all parameters
  result$parameters <- list(
    n_clusters_requested = n_clusters,
    spatial_weight = spatial_weight,
    compactness = compactness,
    max_iterations = max_iterations
  )
  
  result
}

#' Cluster4d using SLIC method
#'
#' Wrapper for SLIC algorithm with standardized interface.
#'
#' @inheritParams cluster4d
#' @param preserve_k Ensure exactly K clusters
#' @param seed_relocate Seed relocation method
#' @param ... Additional parameters passed to slic4d_supervoxels
#'
#' @return A cluster4d_result object
#' @export
cluster4d_slic <- function(vec, mask, n_clusters = 100,
                          spatial_weight = 0.5,
                          max_iterations = 10,
                          connectivity = 26,
                          parallel = TRUE,
                          verbose = FALSE,
                          preserve_k = FALSE,
                          seed_relocate = "none",
                          ...) {
  
  # Convert spatial_weight to compactness
  # SLIC uses larger compactness values (1-20 typical)
  compactness <- spatial_weight * 20
  
  # Determine thread count
  n_threads <- if (parallel) 0 else 1
  
  # Call slic4d_supervoxels
  result <- slic4d_supervoxels(
    bvec = vec,
    mask = mask,
    K = n_clusters,
    compactness = compactness,
    max_iter = max_iterations,
    connectivity = connectivity,
    n_threads = n_threads,
    verbose = verbose,
    preserve_k = preserve_k,
    seed_relocate = seed_relocate,
    ...
  )
  
  # Standardize result structure
  if (!"cluster4d_result" %in% class(result)) {
    class(result) <- c("cluster4d_result", class(result))
  }
  
  # Add method info
  result$method <- "slic"
  result$n_clusters <- length(unique(result$cluster[!is.na(result$cluster)]))
  
  # Store all parameters including those passed through ...
  dots <- list(...)
  result$parameters <- c(
    list(
      n_clusters_requested = n_clusters,
      spatial_weight = spatial_weight,
      compactness = compactness,
      max_iterations = max_iterations,
      connectivity = connectivity,
      parallel = parallel,
      preserve_k = preserve_k,
      seed_relocate = seed_relocate
    ),
    dots  # Include any additional parameters passed through
  )
  
  result
}

#' Cluster4d using slice_msf method
#'
#' Wrapper for slice_msf algorithm with standardized interface.
#'
#' @inheritParams cluster4d
#' @param num_runs Number of independent runs
#' @param consensus Use consensus fusion
#' @param stitch_z Stitch clusters across z-slices
#' @param ... Additional parameters passed to slice_msf
#'
#' @return A cluster4d_result object
#' @export
cluster4d_slice_msf <- function(vec, mask, n_clusters = 100,
                               spatial_weight = 0.5,
                               connectivity = 8,
                               parallel = TRUE,
                               verbose = FALSE,
                               num_runs = 1,
                               consensus = FALSE,
                               stitch_z = TRUE,
                               theta_link = 0.85,
                               min_contact = 1,
                               r = 12,
                               gamma = 1.5,
                               min_size = NULL,
                               ...) {

  # Convert spatial_weight to compactness
  # slice_msf uses medium compactness values (1-10 typical)
  compactness <- spatial_weight * 10

  # Map connectivity (slice_msf uses 2D connectivity)
  nbhd <- if (connectivity == 26) 8 else 4

  # Compute appropriate min_size based on expected cluster size if not specified
  # Allow clusters as small as half the expected size to avoid forced merging
  if (is.null(min_size)) {
    n_voxels <- sum(mask > 0)
    expected_size <- n_voxels / n_clusters
    min_size <- max(1, floor(expected_size * 0.5))
  }

  # Helper to run slice_msf with overrides (used for retries)
  run_once <- function(compactness_val, min_size_val) {
    slice_msf(
      vec = vec,
      mask = mask,
      target_k_global = n_clusters,
      compactness = compactness_val,
      min_size = min_size_val,
      nbhd = nbhd,
      num_runs = num_runs,
      consensus = consensus,
      use_features = consensus,  # Required for exact-K consensus
      stitch_z = stitch_z,
      theta_link = theta_link,
      min_contact = min_contact,
      r = r,
      gamma = gamma,
      ...
    )
  }

  result <- run_once(compactness, min_size)

  # Adaptive retries when we fail to reach or meaningfully overshoot the requested K.
  # We want the initial FH segmentation to produce >K components so RAG merges have
  # room to respect boundaries. If we come in under, tighten fh_scale (higher compactness)
  # and relax min_size.
  if (!inherits(result, "error") && n_clusters > 0) {
    get_stats <- function(res) {
      list(
        found_k = length(unique(res$cluster[!is.na(res$cluster)])),
        n_components_fh = if (!is.null(res$metadata$n_components_fh))
          res$metadata$n_components_fh else NA_integer_
      )
    }

    stats <- get_stats(result)
    attempts <- 0
    while (attempts < 2 &&
           !is.na(stats$found_k) &&
           (stats$found_k < n_clusters ||
            (!is.na(stats$n_components_fh) &&
             stats$n_components_fh < ceiling(1.4 * n_clusters)))) {

      compactness <- min(10, compactness * 2.0)
      min_size <- max(1, floor(min_size * 0.5))
      attempts <- attempts + 1
      if (verbose) message(
        sprintf("slice_msf retry %d: fh tighter (compactness=%.2f) min_size=%d (found %d, fh comps %s)",
                attempts, compactness, min_size,
                stats$found_k,
                ifelse(is.na(stats$n_components_fh), "NA", stats$n_components_fh))
      )
      result <- run_once(compactness, min_size)
      stats <- get_stats(result)
    }
  }
  
  # Standardize result structure
  if (!"cluster4d_result" %in% class(result)) {
    class(result) <- c("cluster4d_result", class(result))
  }
  
  # Add method info
  result$method <- "slice_msf"
  result$n_clusters <- length(unique(result$cluster[!is.na(result$cluster)]))
  if (!is.null(result$centers)) {
    if (nrow(result$centers) != result$n_clusters &&
        ncol(result$centers) == result$n_clusters) {
      result$centers <- t(result$centers)
    }
  }
  
  # Store all parameters including those passed through ...
  dots <- list(...)
  result$parameters <- c(
    list(
      n_clusters_requested = n_clusters,
      spatial_weight = spatial_weight,
      compactness = compactness,
      connectivity = connectivity,
      min_size = min_size,
      num_runs = num_runs,
      consensus = consensus,
      stitch_z = stitch_z,
      theta_link = theta_link,
      min_contact = min_contact,
      r = r,
      gamma = gamma
    ),
    dots  # Include any additional parameters passed through
  )
  
  result
}

#' Cluster4d using FLASH-3D method
#'
#' Wrapper for FLASH-3D algorithm with standardized interface.
#'
#' @inheritParams cluster4d
#' @param lambda_t Temporal weight for Hamming distance
#' @param bits Hash length (64 or 128)
#' @param dctM Number of DCT coefficients
#' @param ... Additional parameters passed to supervoxels_flash3d
#'
#' @return A cluster4d_result object
#' @export
cluster4d_flash3d <- function(vec, mask, n_clusters = 100,
                             spatial_weight = 0.6,
                             max_iterations = 2,
                             verbose = FALSE,
                             lambda_t = 1.0,
                             bits = 64,
                             dctM = 12,
                             ...) {
  
  # Call supervoxels_flash3d
  result <- supervoxels_flash3d(
    vec = vec,
    mask = mask,
    K = n_clusters,
    lambda_s = spatial_weight,
    lambda_t = lambda_t,
    rounds = max_iterations,
    bits = bits,
    dctM = dctM,
    verbose = verbose,
    ...
  )
  
  # Standardize result structure
  if (!"cluster4d_result" %in% class(result)) {
    class(result) <- c("cluster4d_result", class(result))
  }
  
  # Ensure method info
  if (is.null(result$method)) {
    result$method <- "flash3d"
  }
  
  if (is.null(result$n_clusters)) {
    result$n_clusters <- length(unique(result$cluster[!is.na(result$cluster)]))
  }
  
  # Store all parameters including those passed through ...
  dots <- list(...)
  result$parameters <- c(
    list(
      n_clusters_requested = n_clusters,
      spatial_weight = spatial_weight,
      lambda_s = spatial_weight,
      lambda_t = lambda_t,
      max_iterations = max_iterations,
      rounds = max_iterations,
      bits = bits,
      dctM = dctM
    ),
    dots  # Include any additional parameters passed through
  )

  result
}

#' Cluster4d using ACSC method
#'
#' Wrapper for Adaptive Correlation Superclustering algorithm with standardized interface.
#' ACSC uses graph-based clustering with Louvain community detection and optional
#' boundary refinement.
#'
#' @inheritParams cluster4d
#' @param block_size Approximate side length of initial blocks. Default 2.
#' @param refine Logical; whether to refine boundaries. Default TRUE.
#' @param max_refine_iter Maximum iterations for boundary refinement. Default 5.
#' @param ... Additional parameters passed to acsc
#'
#' @return A cluster4d_result object
#' @export
cluster4d_acsc <- function(vec, mask, n_clusters = 100,
                           spatial_weight = 0.5,
                           max_iterations = 5,
                           verbose = FALSE,
                           block_size = 2,
                           refine = TRUE,
                           max_refine_iter = NULL,
                           ...) {

  # Map spatial_weight to alpha (acsc uses alpha for correlation vs spatial balance)
  # alpha = 0.5 means equal weight; higher alpha = more correlation-based
  alpha <- 1 - spatial_weight

  # Use max_iterations for refinement if max_refine_iter not specified

  if (is.null(max_refine_iter)) {
    max_refine_iter <- max_iterations
  }

  # Call acsc
  result <- acsc(
    bvec = vec,
    mask = mask,
    K = n_clusters,
    block_size = block_size,
    alpha = alpha,
    refine = refine,
    max_refine_iter = max_refine_iter,
    ...
  )

  # Standardize result structure - acsc already returns cluster and clusvol
  if (!"cluster4d_result" %in% class(result)) {
    class(result) <- c("cluster4d_result", class(result))
  }

  # Ensure method is set
  result$method <- "acsc"

  # Compute centers if not present
  if (is.null(result$centers)) {
    mask.idx <- which(mask > 0)
    features <- t(series(vec, mask.idx))  # N x T
    coords <- index_to_coord(mask, mask.idx)
    centers_info <- compute_cluster_centers(result$cluster, features, coords)
    result$centers <- centers_info$centers
    result$coord_centers <- centers_info$coord_centers
  }

  # Store parameters
  dots <- list(...)
  result$parameters <- c(
    list(
      n_clusters_requested = n_clusters,
      spatial_weight = spatial_weight,
      alpha = alpha,
      block_size = block_size,
      refine = refine,
      max_refine_iter = max_refine_iter
    ),
    dots
  )

  result
}
