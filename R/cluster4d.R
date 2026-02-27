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
#'     \item \code{"brs_slic"}: Boundary-Refined Sketch SLIC (coarse sketch + boundary exact-correlation refinement)
#'     \item \code{"slice_msf"}: Slice-wise Minimum Spanning Forest (fast but may show z-artifacts)
#'     \item \code{"flash3d"}: Fast Low-rank Approximate Superclusters
#'     \item \code{"g3s"}: Gradient-Guided Geodesic Supervoxels (NEW - recommended for best quality/speed)
#'     \item \code{"rena"}: Recursive Nearest Agglomeration (fast, balanced, topology-aware)
#'     \item \code{"mcl"}: Sparse Markov Cluster Algorithm on a weighted voxel graph
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
#'     \item mcl: direct blend between feature and spatial edge similarities
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
#'   mcl \tab Fast \tab Good \tab Medium \tab Partial \tab Sparse graph clustering with tunable granularity \cr
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
#' \code{\link{cluster4d_mcl}},
#' \code{\link{cluster4d_commute}}
#' 
#' Legacy functions (deprecated): \code{\link{supervoxels}}, \code{\link{snic}},
#' \code{\link{slic4d_supervoxels}}, \code{\link{slice_msf}}, \code{\link{supervoxels_flash3d}}
#'
#' @importFrom methods new
#' @export
#' @importFrom neuroim2 NeuroVec NeuroVol ClusteredNeuroVol series index_to_coord spacing
cluster4d <- function(vec, mask,
                     n_clusters = 100,
                     method = c("supervoxels", "snic", "slic", "corr_slic", "brs_slic", "slice_msf", "flash3d", "g3s", "rena", "rena_plus", "mcl", "acsc", "commute"),
                     spatial_weight = 0.5,
                     max_iterations = NULL,
                     connectivity = NULL,
                     parallel = TRUE,
                     verbose = FALSE,
                     ...) {

  # Allow users to pass a single-volume NeuroVol; wrap to NeuroVec for downstream code
  vec <- ensure_neurovec(vec)

  spatial_weight_missing <- missing(spatial_weight)
  method <- match.arg(method)

  # corr_slic is generally over-regularized with global default 0.5.
  if (spatial_weight_missing && method == "corr_slic") {
    spatial_weight <- 0.1
  }
  if (spatial_weight_missing && method == "brs_slic") {
    spatial_weight <- 0.05
  }
  if (spatial_weight_missing && method == "slice_msf") {
    spatial_weight <- 0.2
  }
  if (spatial_weight_missing && method == "flash3d") {
    spatial_weight <- 0.25
  }
  if (spatial_weight_missing && method == "mcl") {
    spatial_weight <- 0.2
  }

  # ------------------------------
  # Method-specific sane defaults
  # ------------------------------
  if (is.null(max_iterations)) {
    max_iterations <- switch(method,
      supervoxels = 50,   # algorithm expects 30–50 iters for stability
    flash3d    = 4,     # better default quality with modest runtime increase
    snic       = 100,
    slic       = 10,
    corr_slic  = 6,
    brs_slic   = 2,
    slice_msf  = 10,
    g3s        = 5,     # refinement iterations
    rena       = 12,
    rena_plus  = 12,
    mcl        = 8,
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
      corr_slic   = 6,
      brs_slic    = 6,
      slice_msf   = 26,
      rena        = 26,
      rena_plus   = 26,
      mcl         = 6,
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
    corr_slic = cluster4d_corrslic(vec, mask, n_clusters, spatial_weight,
                                   max_iterations, connectivity, parallel, verbose, ...),
    brs_slic = cluster4d_brsslic(vec, mask, n_clusters, spatial_weight,
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
    mcl = cluster4d_mcl(vec, mask, n_clusters, spatial_weight,
                        max_iterations, connectivity, parallel, verbose, ...),
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

#' Cluster4d using correlation-embedded SLIC method
#'
#' Wrapper for a correlation-first SLIC variant that builds a compact random
#' projection embedding of voxel time series, then runs local 3D SLIC updates.
#'
#' @inheritParams cluster4d
#' @param embedding_dim Embedding dimension used to approximate correlations.
#'   Use `"auto"` (or set `adaptive_embedding = TRUE`) for data-adaptive selection.
#' @param adaptive_embedding Logical; if TRUE, choose embedding dimension from
#'   data size and time length.
#' @param embedding_basis Embedding basis, either `"hash"` (CountSketch-like) or
#'   `"dct"` (demeaned DCT basis).
#' @param embedding_whiten Logical; if TRUE, whiten embedding dimensions across voxels.
#' @param sketch_repeats Number of independent hash sketches combined per voxel.
#' @param alpha Compactness weight in the SLIC distance.
#' @param assign_stride Assignment subsampling stride for coarse updates.
#'   1 disables subsampling; values > 1 process one z-slice phase per iteration
#'   followed by a final full assignment pass.
#' @param quantize_assign Logical; if TRUE, quantize embeddings and centroid
#'   features to int8 during assignment for faster dot-product distance checks.
#' @param refine_exact_iters Number of exact-correlation refinement passes after
#'   coarse embedding-SLIC iterations.
#' @param refine_boundary_only Logical; if TRUE, exact refinement is restricted
#'   to boundary voxels each pass.
#' @param refine_stride Time-axis subsampling stride for exact refinement.
#'   Values > 1 speed up refinement by computing correlations on every
#'   `refine_stride`-th time point (approximate).
#' @param refine_alpha Optional compactness for exact refinement. If `NULL`,
#'   uses `alpha`.
#' @param seed Random seed for embedding hash and seed initialization.
#' @param min_size Minimum component size for connectivity enforcement.
#' @param ... Additional arguments (currently unused; reserved for forward compatibility).
#'
#' @return A cluster4d_result object
cluster4d_corrslic <- function(vec, mask, n_clusters = 100,
                               spatial_weight = 0.5,
                               max_iterations = 6,
                               connectivity = 6,
                               parallel = TRUE,
                               verbose = FALSE,
                               embedding_dim = 64L,
                               adaptive_embedding = FALSE,
                               embedding_basis = c("hash", "dct"),
                               embedding_whiten = FALSE,
                               sketch_repeats = 1L,
                               alpha = NULL,
                               assign_stride = 1L,
                               quantize_assign = FALSE,
                               refine_exact_iters = 0L,
                               refine_boundary_only = TRUE,
                               refine_stride = 1L,
                               refine_alpha = NULL,
                               seed = 1L,
                               min_size = NULL,
                               ...) {

  validate_cluster4d_inputs(vec, mask, n_clusters, "cluster4d_corrslic")

  if (!connectivity %in% c(6, 26)) {
    stop("cluster4d_corrslic: connectivity must be 6 or 26")
  }

  # Direct mapping: higher spatial_weight means stronger compactness.
  if (is.null(alpha)) {
    alpha <- spatial_weight
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0) {
    stop("cluster4d_corrslic: alpha must be a positive scalar")
  }

  embedding_basis <- match.arg(embedding_basis)
  if (!is.logical(adaptive_embedding) || length(adaptive_embedding) != 1 || is.na(adaptive_embedding)) {
    stop("cluster4d_corrslic: adaptive_embedding must be TRUE/FALSE")
  }
  if (!is.logical(embedding_whiten) || length(embedding_whiten) != 1 || is.na(embedding_whiten)) {
    stop("cluster4d_corrslic: embedding_whiten must be TRUE/FALSE")
  }
  if (!is.numeric(refine_exact_iters) || length(refine_exact_iters) != 1 || refine_exact_iters < 0) {
    stop("cluster4d_corrslic: refine_exact_iters must be a non-negative integer")
  }
  if (!is.logical(refine_boundary_only) || length(refine_boundary_only) != 1 || is.na(refine_boundary_only)) {
    stop("cluster4d_corrslic: refine_boundary_only must be TRUE/FALSE")
  }
  if (!is.numeric(refine_stride) || length(refine_stride) != 1 || refine_stride < 1) {
    stop("cluster4d_corrslic: refine_stride must be a positive integer")
  }
  if (!is.null(refine_alpha) &&
      (!is.numeric(refine_alpha) || length(refine_alpha) != 1 || !is.finite(refine_alpha) || refine_alpha <= 0)) {
    stop("cluster4d_corrslic: refine_alpha must be NULL or a positive finite scalar")
  }

  if (is.character(embedding_dim)) {
    if (length(embedding_dim) != 1 || !identical(tolower(embedding_dim), "auto")) {
      stop("cluster4d_corrslic: embedding_dim character value must be 'auto'")
    }
    adaptive_embedding <- TRUE
  }
  if (!is.numeric(sketch_repeats) || length(sketch_repeats) != 1 || sketch_repeats < 1) {
    stop("cluster4d_corrslic: sketch_repeats must be a positive integer")
  }
  if (!is.numeric(assign_stride) || length(assign_stride) != 1 || assign_stride < 1) {
    stop("cluster4d_corrslic: assign_stride must be a positive integer")
  }
  if (!is.logical(quantize_assign) || length(quantize_assign) != 1 || is.na(quantize_assign)) {
    stop("cluster4d_corrslic: quantize_assign must be TRUE/FALSE")
  }

  n_threads <- if (parallel) 0L else 1L
  min_size_in <- if (is.null(min_size)) 0L else as.integer(min_size)

  data_prep <- prepare_cluster4d_data(
    vec = vec,
    mask = mask,
    scale_features = FALSE,
    scale_coords = FALSE
  )
  n_time <- ncol(data_prep$features)
  n_vox <- nrow(data_prep$features)

  choose_embedding_dim <- function(n_time, n_vox) {
    max_d <- max(8L, min(128L, as.integer(n_time - 1L)))
    target <- as.integer(round(8 * log2(max(4, n_vox))))
    target <- max(16L, min(max_d, target))
    if (n_time <= 40L) target <- min(target, 32L)
    if (n_time >= 120L) target <- max(target, 64L)
    target <- as.integer(8L * round(target / 8L))
    target <- max(8L, min(max_d, target))
    target
  }

  d_eff <- NULL
  if (isTRUE(adaptive_embedding)) {
    d_eff <- choose_embedding_dim(n_time, n_vox)
  } else {
    if (!is.numeric(embedding_dim) || length(embedding_dim) != 1 || embedding_dim < 8) {
      stop("cluster4d_corrslic: embedding_dim must be >= 8 (or use 'auto')")
    }
    d_eff <- as.integer(embedding_dim)
  }

  core <- corrslic_core(
    feat = data_prep$features,
    mask_lin_idx = as.integer(data_prep$mask_idx) - 1L,
    dims = as.integer(data_prep$dims),
    K = as.integer(n_clusters),
    d = as.integer(d_eff),
    sketch_repeats = as.integer(sketch_repeats),
    alpha = as.numeric(alpha),
    assign_stride = as.integer(assign_stride),
    quantize_assign = isTRUE(quantize_assign),
    embed_basis = embedding_basis,
    whiten_embed = isTRUE(embedding_whiten),
    refine_exact_iters = as.integer(refine_exact_iters),
    refine_boundary_only = isTRUE(refine_boundary_only),
    refine_stride = as.integer(refine_stride),
    refine_alpha = if (is.null(refine_alpha)) -1 else as.numeric(refine_alpha),
    max_iter = as.integer(max_iterations),
    seed = as.integer(seed),
    connectivity = as.integer(connectivity),
    min_size = min_size_in,
    n_threads = as.integer(n_threads),
    verbose = verbose
  )

  labels <- as.integer(core$labels)
  valid <- !is.na(labels) & labels > 0L
  label_values <- sort(unique(labels[valid]))
  if (length(label_values) == 0L) {
    stop("cluster4d_corrslic: no valid labels returned from C++ core")
  }

  grp <- factor(labels, levels = label_values)
  counts <- as.numeric(tabulate(as.integer(grp), nbins = length(label_values)))
  counts[counts == 0] <- NA_real_
  coord_sum <- rowsum(data_prep$coords, grp, reorder = TRUE)
  coord_centers <- coord_sum / counts
  dimnames(coord_centers) <- NULL

  centers <- NULL
  if (!is.null(core$centers) && is.matrix(core$centers)) {
    centers <- core$centers[label_values, , drop = FALSE]
    dimnames(centers) <- NULL
  } else {
    # Fallback for older cores that do not return centers.
    center_info <- compute_cluster_centers(labels, data_prep$features, data_prep$coords, method = "mean")
    centers <- center_info$centers
    coord_centers <- center_info$coord_centers
  }

  centers_xyz_voxel <- NULL
  if (!is.null(core$centers_xyz) && is.matrix(core$centers_xyz)) {
    centers_xyz_voxel <- core$centers_xyz[label_values, , drop = FALSE]
    dimnames(centers_xyz_voxel) <- NULL
  }

  result <- structure(
    list(
      clusvol = ClusteredNeuroVol(mask > 0, clusters = labels),
      cluster = labels,
      centers = centers,
      coord_centers = coord_centers,
      n_clusters = length(label_values),
      method = "corr_slic",
      parameters = list(
        n_clusters_requested = n_clusters,
        spatial_weight = spatial_weight,
        embedding_dim = as.integer(d_eff),
        adaptive_embedding = isTRUE(adaptive_embedding),
        embedding_basis = embedding_basis,
        embedding_whiten = isTRUE(embedding_whiten),
        sketch_repeats = as.integer(sketch_repeats),
        alpha = as.numeric(alpha),
        assign_stride = as.integer(assign_stride),
        quantize_assign = isTRUE(quantize_assign),
        refine_exact_iters = as.integer(refine_exact_iters),
        refine_boundary_only = isTRUE(refine_boundary_only),
        refine_stride = as.integer(refine_stride),
        refine_alpha = if (is.null(refine_alpha)) NULL else as.numeric(refine_alpha),
        max_iterations = as.integer(max_iterations),
        connectivity = as.integer(connectivity),
        parallel = isTRUE(parallel),
        seed = as.integer(seed),
        min_size = if (is.null(min_size)) NULL else as.integer(min_size)
      ),
      metadata = list(
        cpp_params = core$params,
        label_values = label_values,
        centers_xyz_voxel = centers_xyz_voxel
      )
    ),
    class = c("cluster4d_result", "cluster_result", "list")
  )

  result
}

#' Cluster4d using boundary-refined sketch SLIC
#'
#' Hybrid method: fast sketch-SLIC coarse assignment, then exact-correlation
#' refinement only on boundary voxels.
#'
#' @inheritParams cluster4d
#' @param embedding_dim Embedding dimension for coarse sketch stage.
#' @param sketch_repeats Number of independent sketches to average.
#' @param coarse_alpha Compactness weight in the coarse sketch-SLIC stage.
#' @param boundary_passes Number of boundary-only refinement passes.
#' @param global_passes Number of local-global refinement passes after boundary passes.
#' @param refine_spatial_weight Spatial regularization during boundary refinement.
#' @param refine_l2_weight Optional L2 prototype penalty during boundary refinement.
#' @param refine_stride Optional temporal stride for refinement (NULL = auto).
#' @param seed Random seed.
#' @param min_size Minimum component size for connectivity enforcement.
#' @param ... Reserved for future extensions.
#'
#' @return A cluster4d_result object
cluster4d_brsslic <- function(vec, mask, n_clusters = 100,
                              spatial_weight = 0.05,
                              max_iterations = 2,
                              connectivity = 6,
                              parallel = TRUE,
                              verbose = FALSE,
                              embedding_dim = 24L,
                              sketch_repeats = 1L,
                              coarse_alpha = NULL,
                              boundary_passes = 1L,
                              global_passes = 0L,
                              refine_spatial_weight = 0,
                              refine_l2_weight = 0,
                              refine_stride = NULL,
                              seed = 1L,
                              min_size = NULL,
                              ...) {

  validate_cluster4d_inputs(vec, mask, n_clusters, "cluster4d_brsslic")

  if (!connectivity %in% c(6, 26)) {
    stop("cluster4d_brsslic: connectivity must be 6 or 26")
  }

  if (is.null(coarse_alpha)) {
    coarse_alpha <- spatial_weight
  }

  if (!is.numeric(coarse_alpha) || length(coarse_alpha) != 1 || coarse_alpha <= 0) {
    stop("cluster4d_brsslic: coarse_alpha must be a positive scalar")
  }
  if (!is.numeric(refine_spatial_weight) || length(refine_spatial_weight) != 1 || refine_spatial_weight < 0) {
    stop("cluster4d_brsslic: refine_spatial_weight must be >= 0")
  }
  if (!is.numeric(refine_l2_weight) || length(refine_l2_weight) != 1 || refine_l2_weight < 0) {
    stop("cluster4d_brsslic: refine_l2_weight must be >= 0")
  }
  if (!is.null(refine_stride) && (!is.numeric(refine_stride) || length(refine_stride) != 1 || refine_stride < 1)) {
    stop("cluster4d_brsslic: refine_stride must be NULL or >= 1")
  }
  if (!is.numeric(embedding_dim) || embedding_dim < 8) {
    stop("cluster4d_brsslic: embedding_dim must be >= 8")
  }
  if (!is.numeric(sketch_repeats) || length(sketch_repeats) != 1 || sketch_repeats < 1) {
    stop("cluster4d_brsslic: sketch_repeats must be a positive integer")
  }
  if (!is.numeric(boundary_passes) || length(boundary_passes) != 1 || boundary_passes < 0) {
    stop("cluster4d_brsslic: boundary_passes must be >= 0")
  }
  if (!is.numeric(global_passes) || length(global_passes) != 1 || global_passes < 0) {
    stop("cluster4d_brsslic: global_passes must be >= 0")
  }

  n_threads <- if (parallel) 0L else 1L
  min_size_in <- if (is.null(min_size)) 0L else as.integer(min_size)

  data_prep <- prepare_cluster4d_data(
    vec = vec,
    mask = mask,
    scale_features = FALSE,
    scale_coords = FALSE
  )

  core <- brs_slic_core(
    feat = data_prep$features,
    mask_lin_idx = as.integer(data_prep$mask_idx) - 1L,
    dims = as.integer(data_prep$dims),
    K = as.integer(n_clusters),
    d = as.integer(embedding_dim),
    sketch_repeats = as.integer(sketch_repeats),
    alpha = as.numeric(coarse_alpha),
    coarse_iter = as.integer(max_iterations),
    boundary_passes = as.integer(boundary_passes),
    global_passes = as.integer(global_passes),
    refine_spatial = as.numeric(refine_spatial_weight),
    refine_l2 = as.numeric(refine_l2_weight),
    refine_stride = if (is.null(refine_stride)) 0L else as.integer(refine_stride),
    seed = as.integer(seed),
    connectivity = as.integer(connectivity),
    min_size = min_size_in,
    n_threads = as.integer(n_threads),
    verbose = verbose
  )

  labels <- as.integer(core$labels)
  valid <- !is.na(labels) & labels > 0L
  label_values <- sort(unique(labels[valid]))
  if (length(label_values) == 0L) {
    stop("cluster4d_brsslic: no valid labels returned from C++ core")
  }

  grp <- factor(labels, levels = label_values)
  counts <- as.numeric(tabulate(as.integer(grp), nbins = length(label_values)))
  counts[counts == 0] <- NA_real_
  coord_sum <- rowsum(data_prep$coords, grp, reorder = TRUE)
  coord_centers <- coord_sum / counts
  dimnames(coord_centers) <- NULL

  centers <- NULL
  if (!is.null(core$centers) && is.matrix(core$centers)) {
    centers <- core$centers[label_values, , drop = FALSE]
    dimnames(centers) <- NULL
  } else {
    center_info <- compute_cluster_centers(labels, data_prep$features, data_prep$coords, method = "mean")
    centers <- center_info$centers
    coord_centers <- center_info$coord_centers
  }

  centers_xyz_voxel <- NULL
  if (!is.null(core$centers_xyz) && is.matrix(core$centers_xyz)) {
    centers_xyz_voxel <- core$centers_xyz[label_values, , drop = FALSE]
    dimnames(centers_xyz_voxel) <- NULL
  }

  result <- structure(
    list(
      clusvol = ClusteredNeuroVol(mask > 0, clusters = labels),
      cluster = labels,
      centers = centers,
      coord_centers = coord_centers,
      n_clusters = length(label_values),
      method = "brs_slic",
      parameters = list(
        n_clusters_requested = n_clusters,
        spatial_weight = spatial_weight,
        embedding_dim = as.integer(embedding_dim),
        sketch_repeats = as.integer(sketch_repeats),
        coarse_alpha = as.numeric(coarse_alpha),
        coarse_iter = as.integer(max_iterations),
        boundary_passes = as.integer(boundary_passes),
        global_passes = as.integer(global_passes),
        refine_spatial_weight = as.numeric(refine_spatial_weight),
        refine_l2_weight = as.numeric(refine_l2_weight),
        refine_stride = if (is.null(refine_stride)) NULL else as.integer(refine_stride),
        connectivity = as.integer(connectivity),
        parallel = isTRUE(parallel),
        seed = as.integer(seed),
        min_size = if (is.null(min_size)) NULL else as.integer(min_size)
      ),
      metadata = list(
        cpp_params = core$params,
        label_values = label_values,
        centers_xyz_voxel = centers_xyz_voxel
      )
    ),
    class = c("cluster4d_result", "cluster_result", "list")
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
#' @param theta_link Minimum cosine similarity for stitching across z-slices.
#' @param min_contact Minimum number of contacting voxels for z-stitch edges.
#' @param r Neighborhood radius for local feature extraction.
#' @param gamma Gamma parameter for edge weighting.
#' @param z_mult Optional z-axis edge multiplier passed to slice_msf (0-1).
#' @param min_size Minimum component size after connectivity enforcement. Default NULL (auto).
#' @param ... Additional parameters passed to slice_msf
#'
#' @return A cluster4d_result object
#' @export
cluster4d_slice_msf <- function(vec, mask, n_clusters = 100,
                               spatial_weight = 0.2,
                               connectivity = 26,
                               parallel = TRUE,
                               verbose = FALSE,
                               num_runs = 1,
                               consensus = FALSE,
                               stitch_z = TRUE,
                               theta_link = 0.85,
                               min_contact = 1,
                               r = 12,
                               gamma = 1.0,
                               z_mult = 0.5,
                               min_size = NULL,
                               ...) {

  # Convert spatial_weight to compactness
  # slice_msf uses medium compactness values (1-10 typical)
  compactness <- spatial_weight * 10

  # Map connectivity (slice_msf uses 2D connectivity)
  nbhd <- if (connectivity == 26) 8 else 4

  # Compute a conservative default min_size from expected cluster size.
  # Keep this low so exact-K merging can preserve boundaries on heterogeneous data.
  if (is.null(min_size)) {
    n_voxels <- sum(mask > 0)
    expected_size <- n_voxels / n_clusters
    min_size <- max(2L, min(20L, floor(expected_size * 0.1)))
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
      z_mult = z_mult,
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
      gamma = gamma,
      z_mult = z_mult
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
                             spatial_weight = 0.25,
                             max_iterations = 4,
                             verbose = FALSE,
                             lambda_t = 1.4,
                             bits = 64,
                             dctM = 16,
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
