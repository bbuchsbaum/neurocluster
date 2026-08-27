#' Unified 4D Clustering for Neuroimaging Data
#'
#' Performs spatially-constrained clustering on 4D neuroimaging data using
#' various algorithms. This is the main entry point for all clustering methods
#' in the neurocluster package.
#'
#' @param vec A \code{NeuroVec} instance supplying the 4D data (x, y, z, time) to cluster
#' @param mask A \code{NeuroVol} mask defining the voxels to include in clustering.
#'   Values must be finite. Voxels are included exactly when the mask value is
#'   strictly greater than zero; zero and negative values are excluded.
#' @param n_clusters Finite integer target number of clusters (default 100), from
#'   one through the number of included mask voxels. Note that some methods
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
#'   Both endpoints are supported: 0 disables the spatial term and 1 disables
#'   the feature term. Higher values emphasize spatial compactness. Default 0.5.
#'   This parameter is inactive for \code{"rena"} and \code{"rena_plus"};
#'   supplying it explicitly for either method is an error.
#'   Maps to method-specific parameters:
#'   \itemize{
#'     \item supervoxels: \code{alpha = 1 - spatial_weight} (0 = all spatial, 1 = all feature)
#'     \item snic: \code{compactness = spatial_weight * 10} (range 0-10)
#'     \item slic: \code{compactness = spatial_weight * 20} (typical range 1-20)
#'     \item corr_slic/brs_slic: direct convex blend between correlation and
#'       scaled spatial distance
#'     \item slice_msf: \code{compactness = spatial_weight * 10} (typical range 1-10)
#'     \item flash3d: \code{lambda_s = spatial_weight} (direct mapping)
#'     \item mcl: direct blend between feature and spatial edge similarities
#'   }
#' @param max_iterations Positive finite integer iteration limit for iterative
#'   methods. It is inactive for \code{"snic"}, \code{"slice_msf"}, and
#'   \code{"commute"}; supplying it explicitly for those methods is an error.
#' @param connectivity Integer neighborhood connectivity. Supported values are
#'   method-specific: supervoxels accepts 6, 18, 26, or 27; slic, corr_slic,
#'   brs_slic, and slice_msf accept 6 or 26; G3S, ReNA variants, and MCL accept
#'   6, 18, or 26. It is inactive for snic, flash3d, acsc, and commute.
#' @param parallel Enable parallel processing for supervoxels and slic. For all
#'   other methods, an explicit \code{TRUE} is rejected; omitted or
#'   \code{FALSE} selects the serial implementation. Default TRUE for supported
#'   methods.
#' @param verbose Print progress messages. Default FALSE.
#' @param ... Additional method-specific parameters. See method documentation for details.
#'
#' @return A \code{cluster4d_result} object (also inherits from \code{cluster_result})
#'   containing:
#'   \item{labels}{Contiguous positive integer final labels in mask voxel order}
#'   \item{cluster}{Backward-compatible alias that is identical to \code{labels}}
#'   \item{clusvol}{A \code{ClusteredNeuroVol} with cluster assignments}
#'   \item{centers}{Actual final-label means in the original feature space,
#'     with shape actual K by timepoints}
#'   \item{coord_centers}{Actual final-label coordinate means in physical mm,
#'     with shape actual K by 3}
#'   \item{actual_k}{Actual number of clusters produced}
#'   \item{n_clusters}{Backward-compatible alias of \code{actual_k}}
#'   \item{label_ids}{Explicit mapping from center rows to contiguous label IDs}
#'   \item{method}{Clustering method used}
#'   \item{parameters}{List of all parameters used}
#'   \item{provenance}{Typed label, original-feature, physical-coordinate,
#'     geometry, and mask provenance}
#'   \item{metadata}{Method-specific additional information}
#'
#' @section Algorithm Comparison:
#'
#' \tabular{lllllll}{
#'   \strong{Method} \tab \strong{Speed} \tab \strong{3D Continuity} \tab \strong{Memory} \tab \strong{Parallel} \tab \strong{Best For} \cr
#'   supervoxels \tab Slow \tab Excellent \tab High \tab Yes \tab Small-medium data, smooth parcels \cr
#'   snic \tab Fast \tab Good \tab Low \tab No \tab Large data, non-iterative \cr
#'   slic \tab Fast \tab Good \tab Medium \tab Yes \tab Balanced speed/quality \cr
#'   slice_msf \tab Very Fast \tab Moderate \tab Low \tab No \tab High-res data, accept z-artifacts \cr
#'   flash3d \tab Fast \tab Good \tab Medium \tab No \tab Large data, hash-based \cr
#'   rena \tab Fast \tab Excellent \tab Low \tab No \tab Balanced clusters, topology-aware \cr
#'   mcl \tab Fast \tab Good \tab Medium \tab No \tab Sparse graph clustering with tunable granularity \cr
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
#' @importFrom methods new as
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

  spatial_weight_missing <- missing(spatial_weight)
  max_iterations_missing <- missing(max_iterations) || is.null(max_iterations)
  connectivity_missing <- missing(connectivity) || is.null(connectivity)
  parallel_missing <- missing(parallel)

  # Allow users to pass a single-volume NeuroVol; wrap to NeuroVec for downstream code
  vec <- ensure_neurovec(vec)

  method <- match.arg(method)
  method_contract <- cluster4d_method_contract(method)
  method_parameters <- list(...)

  # corr_slic is generally over-regularized with global default 0.5.
  if (spatial_weight_missing && method == "corr_slic") {
    spatial_weight <- if (!is.null(method_parameters$alpha)) {
      method_parameters$alpha
    } else 0.1
  }
  if (spatial_weight_missing && method == "brs_slic") {
    spatial_weight <- if (!is.null(method_parameters$coarse_alpha)) {
      method_parameters$coarse_alpha
    } else 0.05
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
      snic        = NA_integer_,
      slic        = 26,
      corr_slic   = 6,
      brs_slic    = 6,
      slice_msf   = 26,
      rena        = 26,
      rena_plus   = 26,
      mcl         = 6,
      g3s         = 26,
      acsc        = NA_integer_,
      flash3d     = NA_integer_,
      NA_integer_
    )
  }

  # Validate common inputs
  contract <- validate_cluster4d_inputs(
    vec, mask, n_clusters, paste0("cluster4d:", method)
  )
  n_clusters <- contract$n_clusters

  if (!method_contract$spatial_weight && !spatial_weight_missing) {
    stop("cluster4d:", method, ": spatial_weight is not supported", call. = FALSE)
  }
  if (method_contract$spatial_weight) {
    spatial_weight <- .cluster4d_scalar_number(
      spatial_weight, "spatial_weight", paste0("cluster4d:", method)
    )
    if (spatial_weight < 0 || spatial_weight > 1) {
      stop("spatial_weight must be between 0 and 1", call. = FALSE)
    }
  } else {
    spatial_weight <- 0
  }

  if (!method_contract$iterations && !max_iterations_missing) {
    stop("cluster4d:", method, ": max_iterations is not supported", call. = FALSE)
  }
  max_iterations <- .cluster4d_scalar_number(
    max_iterations, "max_iterations", paste0("cluster4d:", method),
    lower = 1, integer = TRUE
  )

  if (is.null(method_contract$connectivity)) {
    if (!connectivity_missing) {
      stop("cluster4d:", method, ": connectivity is not supported", call. = FALSE)
    }
    connectivity <- NA_integer_
  } else {
    connectivity <- .cluster4d_scalar_number(
      connectivity, "connectivity", paste0("cluster4d:", method),
      lower = 1, integer = TRUE
    )
    if (!connectivity %in% method_contract$connectivity) {
      if (identical(method_contract$connectivity, c(6L, 18L, 26L, 27L))) {
        stop("connectivity must be 6, 18, 26, or 27", call. = FALSE)
      } else {
        stop(
          "cluster4d:", method, ": connectivity must be one of ",
          paste(method_contract$connectivity, collapse = ", "),
          call. = FALSE
        )
      }
    }
  }

  if (!method_contract$parallel && !parallel_missing) {
    requested_parallel <- .cluster4d_scalar_logical(
      parallel, "parallel", paste0("cluster4d:", method)
    )
    if (requested_parallel) {
      stop("cluster4d:", method, ": parallel is not supported", call. = FALSE)
    }
  }
  parallel <- if (method_contract$parallel) {
    .cluster4d_scalar_logical(parallel, "parallel", paste0("cluster4d:", method))
  } else FALSE
  verbose <- .cluster4d_scalar_logical(
    verbose, "verbose", paste0("cluster4d:", method)
  )

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
  if (!method_contract$spatial_weight) orig_params$spatial_weight <- NULL
  if (!method_contract$iterations) orig_params$max_iterations <- NULL
  if (is.null(method_contract$connectivity)) orig_params$connectivity <- NULL
  if (!method_contract$parallel) orig_params$parallel <- NULL

  # Dispatch to method-specific implementation
  result <- switch(method,
    supervoxels = cluster4d_supervoxels(vec, mask, n_clusters, spatial_weight,
                                        max_iterations, connectivity, parallel,
                                        verbose, ...),
    snic = cluster4d_snic(
      vec, mask, n_clusters, spatial_weight,
      verbose = verbose, ...
    ),
    slic = cluster4d_slic(
      vec, mask, n_clusters, spatial_weight,
      max_iterations, connectivity, parallel, verbose, ...,
      .input_contract = contract
    ),
    corr_slic = cluster4d_corrslic(
      vec, mask, n_clusters, spatial_weight,
      max_iterations, connectivity, parallel, verbose, ...,
      .input_contract = contract
    ),
    brs_slic = cluster4d_brsslic(vec, mask, n_clusters, spatial_weight,
                                 max_iterations, connectivity, parallel, verbose, ...),
    slice_msf = cluster4d_slice_msf(
      vec, mask, n_clusters, spatial_weight, connectivity,
      verbose = verbose, ...
    ),
    flash3d = cluster4d_flash3d(vec, mask, n_clusters, spatial_weight,
                               max_iterations, verbose, ...),
    g3s = cluster4d_g3s(vec, mask, n_clusters,
                       alpha = 1 - spatial_weight,  # Convert to feature weight
                       max_refinement_iter = max_iterations,
                       connectivity = connectivity,
                       verbose = verbose, ...),
    rena = cluster4d_rena(
      vec, mask, n_clusters,
      max_iterations = max_iterations,
      connectivity = connectivity,
      verbose = verbose, ...
    ),
    rena_plus = cluster4d_rena_plus(
      vec, mask, n_clusters,
      connectivity = connectivity,
      max_iterations = max_iterations,
      verbose = verbose, ...
    ),
    mcl = cluster4d_mcl(
      vec, mask, n_clusters, spatial_weight,
      max_iterations, connectivity,
      verbose = verbose, ...
    ),
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
  if (!method_contract$spatial_weight) result$parameters$spatial_weight <- NULL
  if (!method_contract$iterations) result$parameters$max_iterations <- NULL
  if (is.null(method_contract$connectivity)) result$parameters$connectivity <- NULL
  if (!method_contract$parallel) result$parameters$parallel <- NULL

  if (identical(result$provenance$schema_version, "1.0.0")) {
    result$parameters <- .cluster4d_merge_parameters(
      result$parameters, orig_params
    )
    result
  } else {
    finalize_cluster4d_result(
      result = result,
      vec = vec,
      mask = mask,
      method = method,
      parameters = orig_params,
      input_contract = contract
    )
  }
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

  validate_cluster4d_inputs(vec, mask, n_clusters, "cluster4d_commute")

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

  finalize_cluster4d_result(result, vec, mask, "commute", result$parameters)
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

  validate_cluster4d_inputs(vec, mask, n_clusters, "cluster4d_supervoxels")
  
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
  
  finalize_cluster4d_result(result, vec, mask, "supervoxels", result$parameters)
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

  if (!missing(max_iterations)) {
    stop("cluster4d_snic: max_iterations is not supported", call. = FALSE)
  }
  validate_cluster4d_inputs(vec, mask, n_clusters, "cluster4d_snic")
  
  # Direct endpoint-preserving mapping: 0 is feature-only and 10 spatial-only.
  compactness <- spatial_weight * 10
  
  # Call original snic
  result <- snic(
    vec = vec,
    mask = mask,
    K = n_clusters,
    compactness = compactness,
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
  result$parameters <- .cluster4d_merge_parameters(
    result$parameters,
    list(
      n_clusters_requested = n_clusters,
      spatial_weight = spatial_weight,
      compactness = compactness
    )
  )
  
  finalize_cluster4d_result(result, vec, mask, "snic", result$parameters)
}

#' Cluster4d using SLIC method
#'
#' Wrapper for SLIC algorithm with standardized interface.
#'
#' @inheritParams cluster4d
#' @param preserve_k Request exactly K connected clusters. Exact K is feasible
#'   only when K is at least the number of connected mask components; below
#'   that topological minimum, connectivity takes precedence and the minimum
#'   feasible cluster count is returned with a warning.
#' @param seed_relocate Seed relocation method
#' @param ... Additional parameters passed to slic4d_supervoxels
#' @param .input_contract Internal prevalidated input receipt used by
#'   `cluster4d()`; direct callers should leave this as `NULL`.
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
                          ...,
                          .input_contract = NULL) {

  input_contract <- .cluster4d_resolve_input_contract(
    .input_contract, vec, mask, n_clusters, "cluster4d_slic"
  )
  n_clusters <- input_contract$n_clusters
  
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
    ...,
    .input_contract = input_contract
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
#' @param alpha Optional backward-compatible alias for the spatial blend weight.
#'   It must be in `[0, 1]`: 0 is correlation-only and 1 is spatial-only.
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
#' @param refine_alpha Optional spatial blend weight for exact refinement. If
#'   `NULL`, uses `alpha`; both endpoints are supported.
#' @param seed Random seed for embedding hash and seed initialization.
#' @param min_size Minimum component size for connectivity enforcement.
#' @param ... Additional arguments (currently unused; reserved for forward compatibility).
#' @param .input_contract Internal prevalidated input receipt used by
#'   `cluster4d()`; direct callers should leave this as `NULL`.
#'
#' @details `spatial_weight` (or its alias `alpha`) is a convex blend:
#' `(1 - w) * correlation_distance + w * scaled_spatial_distance`.
#' Coarse correlation is approximated by a hash or DCT sketch. Assignment may
#' additionally use int8 quantization or z-plane subsampling; exact refinement
#' may use a temporal stride. These choices and their stage are recorded under
#' `metadata$approximation`. The public `centers` are always recomputed from the
#' final labels in the original time-series space; sketch centers are separately
#' typed metadata and are never substituted for them.
#'
#' @return A cluster4d_result object
cluster4d_corrslic <- function(vec, mask, n_clusters = 100,
                               spatial_weight = 0.5,
                               max_iterations = 6,
                               connectivity = 6,
                               parallel = FALSE,
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
                               ...,
                               .input_contract = NULL) {

  spatial_weight_missing <- missing(spatial_weight)
  input_contract <- .cluster4d_resolve_input_contract(
    .input_contract, vec, mask, n_clusters, "cluster4d_corrslic"
  )
  n_clusters <- input_contract$n_clusters
  spatial_weight <- .cluster4d_scalar_number(
    spatial_weight, "spatial_weight", "cluster4d_corrslic", lower = 0, upper = 1
  )
  max_iterations <- .cluster4d_scalar_number(
    max_iterations, "max_iterations", "cluster4d_corrslic",
    lower = 1, upper = 1000, integer = TRUE
  )
  connectivity <- .cluster4d_scalar_number(
    connectivity, "connectivity", "cluster4d_corrslic", integer = TRUE
  )
  if (!connectivity %in% c(6L, 26L)) {
    stop("cluster4d_corrslic: connectivity must be 6 or 26", call. = FALSE)
  }
  parallel <- .cluster4d_scalar_logical(
    parallel, "parallel", "cluster4d_corrslic"
  )
  if (parallel) {
    stop(
      "cluster4d_corrslic: parallel is not supported by the package build",
      call. = FALSE
    )
  }
  verbose <- .cluster4d_scalar_logical(
    verbose, "verbose", "cluster4d_corrslic"
  )

  # Direct mapping: higher spatial_weight means stronger compactness.
  if (is.null(alpha)) {
    alpha <- spatial_weight
  }

  alpha <- .cluster4d_scalar_number(
    alpha, "alpha", "cluster4d_corrslic", lower = 0, upper = 1
  )
  if (!spatial_weight_missing && !isTRUE(all.equal(alpha, spatial_weight))) {
    stop(
      "cluster4d_corrslic: alpha and spatial_weight must agree when both are supplied",
      call. = FALSE
    )
  }
  spatial_weight <- alpha

  embedding_basis <- match.arg(embedding_basis)
  adaptive_embedding <- .cluster4d_scalar_logical(
    adaptive_embedding, "adaptive_embedding", "cluster4d_corrslic"
  )
  embedding_whiten <- .cluster4d_scalar_logical(
    embedding_whiten, "embedding_whiten", "cluster4d_corrslic"
  )
  refine_exact_iters <- .cluster4d_scalar_number(
    refine_exact_iters, "refine_exact_iters", "cluster4d_corrslic",
    lower = 0, upper = 16, integer = TRUE
  )
  refine_boundary_only <- .cluster4d_scalar_logical(
    refine_boundary_only, "refine_boundary_only", "cluster4d_corrslic"
  )
  refine_stride <- .cluster4d_scalar_number(
    refine_stride, "refine_stride", "cluster4d_corrslic",
    lower = 1, upper = 64, integer = TRUE
  )
  if (!is.null(refine_alpha)) {
    refine_alpha <- .cluster4d_scalar_number(
      refine_alpha, "refine_alpha", "cluster4d_corrslic", lower = 0, upper = 1
    )
  }

  if (is.character(embedding_dim)) {
    if (length(embedding_dim) != 1 || !identical(tolower(embedding_dim), "auto")) {
      stop("cluster4d_corrslic: embedding_dim character value must be 'auto'")
    }
    adaptive_embedding <- TRUE
  }
  sketch_repeats <- .cluster4d_scalar_number(
    sketch_repeats, "sketch_repeats", "cluster4d_corrslic",
    lower = 1, upper = 8, integer = TRUE
  )
  assign_stride <- .cluster4d_scalar_number(
    assign_stride, "assign_stride", "cluster4d_corrslic",
    lower = 1, upper = 16, integer = TRUE
  )
  quantize_assign <- .cluster4d_scalar_logical(
    quantize_assign, "quantize_assign", "cluster4d_corrslic"
  )
  seed <- .cluster4d_scalar_number(
    seed, "seed", "cluster4d_corrslic",
    lower = -.Machine$integer.max, upper = .Machine$integer.max, integer = TRUE
  )
  if (!is.null(min_size)) {
    min_size <- .cluster4d_scalar_number(
      min_size, "min_size", "cluster4d_corrslic", lower = 0, integer = TRUE
    )
  }

  n_threads <- 1L
  min_size_in <- if (is.null(min_size)) 0L else as.integer(min_size)

  data_prep <- prepare_cluster4d_data(
    vec = vec,
    mask = mask,
    scale_features = FALSE,
    scale_coords = FALSE,
    input_contract = input_contract
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
    d_eff <- .cluster4d_scalar_number(
      embedding_dim, "embedding_dim", "cluster4d_corrslic",
      lower = 8, integer = TRUE
    )
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
  label_values <- as.integer(core$label_ids)
  if (!identical(label_values, seq_len(as.integer(core$actual_k))) ||
      !identical(sort(unique(labels)), label_values)) {
    stop("cluster4d_corrslic: C++ core returned a non-contiguous final label map")
  }

  centers <- unname(as.matrix(core$original_centers))
  centers_xyz_voxel <- unname(as.matrix(core$centers_xyz))

  result <- structure(
    list(
      labels = labels,
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
        seed = as.integer(seed),
        min_size = if (is.null(min_size)) NULL else as.integer(min_size)
      ),
      metadata = list(
        cpp_params = core$params,
        label_values = label_values,
        centers_xyz_voxel = centers_xyz_voxel,
        approximation = list(
          schema_version = "1.0.0",
          family = "correlation_sketch_slic",
          embedding = list(
            basis = embedding_basis,
            dimension = as.integer(d_eff),
            repeats = as.integer(sketch_repeats),
            whitened = isTRUE(embedding_whiten),
            centers = unname(as.matrix(core$embedding_centers)),
            row_to_label = as.integer(core$label_ids),
            statistic = "final_label_mean_of_unit_sketch_vectors"
          ),
          assignment = list(
            stride = as.integer(assign_stride),
            subsampled_iterations = as.integer(core$params$assign_subsample_iters),
            quantized_int8 = isTRUE(quantize_assign)
          ),
          refinement = list(
            exact_passes = as.integer(refine_exact_iters),
            boundary_only = isTRUE(refine_boundary_only),
            temporal_stride = as.integer(refine_stride),
            spatial_weight = if (is.null(refine_alpha)) as.numeric(alpha) else as.numeric(refine_alpha)
          ),
          native_final = list(
            original_centers = centers,
            coord_centers_voxel = centers_xyz_voxel,
            label_ids = as.integer(core$label_ids),
            original_label_ids = as.integer(core$original_label_ids),
            counts = as.integer(core$counts)
          )
        )
      )
    ),
    class = c("cluster4d_result", "cluster_result", "list")
  )

  finalize_cluster4d_result(
    result, vec, mask, "corr_slic", result$parameters,
    input_contract = input_contract,
    data = data_prep
  )
}

#' Cluster4d using boundary-refined sketch SLIC
#'
#' Hybrid method: fast sketch-SLIC coarse assignment, then exact-correlation
#' refinement only on boundary voxels.
#'
#' @inheritParams cluster4d
#' @param embedding_dim Embedding dimension for coarse sketch stage.
#' @param sketch_repeats Number of independent sketches to average.
#' @param coarse_alpha Optional backward-compatible alias for the coarse spatial
#'   blend weight in `[0, 1]`.
#' @param boundary_passes Number of boundary-only refinement passes.
#' @param global_passes Number of local-global refinement passes after boundary passes.
#' @param refine_spatial_weight Optional spatial blend weight during exact
#'   refinement. `NULL` inherits `coarse_alpha`; both endpoints are supported.
#' @param refine_l2_weight Optional L2 prototype penalty during boundary refinement.
#' @param refine_stride Optional temporal stride for refinement (NULL = auto).
#' @param seed Random seed.
#' @param min_size Minimum component size for connectivity enforcement.
#' @param ... Reserved for future extensions.
#'
#' @details Both the coarse and exact-refinement spatial weights use the convex
#' blend `(1 - w) * correlation_distance + w * scaled_spatial_distance`.
#' `refine_l2_weight` is a separate optional feature penalty, so a strictly
#' spatial-only refinement uses `refine_spatial_weight = 1` and
#' `refine_l2_weight = 0`. The hash-sketch artifact is explicitly marked as a
#' pre-refinement embedding under `metadata$approximation$coarse_embedding`.
#' Public centers and coordinates are final-label means in original feature and
#' physical coordinate spaces, respectively.
#'
#' @return A cluster4d_result object
cluster4d_brsslic <- function(vec, mask, n_clusters = 100,
                              spatial_weight = 0.05,
                              max_iterations = 2,
                              connectivity = 6,
                              parallel = FALSE,
                              verbose = FALSE,
                              embedding_dim = 24L,
                              sketch_repeats = 1L,
                              coarse_alpha = NULL,
                              boundary_passes = 1L,
                              global_passes = 0L,
                              refine_spatial_weight = NULL,
                              refine_l2_weight = 0,
                              refine_stride = NULL,
                              seed = 1L,
                              min_size = NULL,
                              ...) {

  spatial_weight_missing <- missing(spatial_weight)
  input_contract <- validate_cluster4d_inputs(
    vec, mask, n_clusters, "cluster4d_brsslic"
  )
  n_clusters <- input_contract$n_clusters
  spatial_weight <- .cluster4d_scalar_number(
    spatial_weight, "spatial_weight", "cluster4d_brsslic", lower = 0, upper = 1
  )
  max_iterations <- .cluster4d_scalar_number(
    max_iterations, "max_iterations", "cluster4d_brsslic",
    lower = 1, upper = 1000, integer = TRUE
  )
  connectivity <- .cluster4d_scalar_number(
    connectivity, "connectivity", "cluster4d_brsslic", integer = TRUE
  )
  if (!connectivity %in% c(6L, 26L)) {
    stop("cluster4d_brsslic: connectivity must be 6 or 26", call. = FALSE)
  }
  parallel <- .cluster4d_scalar_logical(
    parallel, "parallel", "cluster4d_brsslic"
  )
  if (parallel) {
    stop(
      "cluster4d_brsslic: parallel is not supported by the package build",
      call. = FALSE
    )
  }
  verbose <- .cluster4d_scalar_logical(
    verbose, "verbose", "cluster4d_brsslic"
  )

  if (is.null(coarse_alpha)) {
    coarse_alpha <- spatial_weight
  }

  coarse_alpha <- .cluster4d_scalar_number(
    coarse_alpha, "coarse_alpha", "cluster4d_brsslic", lower = 0, upper = 1
  )
  if (!spatial_weight_missing && !isTRUE(all.equal(coarse_alpha, spatial_weight))) {
    stop(
      "cluster4d_brsslic: coarse_alpha and spatial_weight must agree when both are supplied",
      call. = FALSE
    )
  }
  spatial_weight <- coarse_alpha
  if (is.null(refine_spatial_weight)) refine_spatial_weight <- coarse_alpha
  refine_spatial_weight <- .cluster4d_scalar_number(
    refine_spatial_weight, "refine_spatial_weight", "cluster4d_brsslic",
    lower = 0, upper = 1
  )
  refine_l2_weight <- .cluster4d_scalar_number(
    refine_l2_weight, "refine_l2_weight", "cluster4d_brsslic", lower = 0
  )
  if (!is.null(refine_stride)) {
    refine_stride <- .cluster4d_scalar_number(
      refine_stride, "refine_stride", "cluster4d_brsslic",
      lower = 1, upper = 64, integer = TRUE
    )
  }
  embedding_dim <- .cluster4d_scalar_number(
    embedding_dim, "embedding_dim", "cluster4d_brsslic", lower = 8, integer = TRUE
  )
  sketch_repeats <- .cluster4d_scalar_number(
    sketch_repeats, "sketch_repeats", "cluster4d_brsslic",
    lower = 1, upper = 8, integer = TRUE
  )
  boundary_passes <- .cluster4d_scalar_number(
    boundary_passes, "boundary_passes", "cluster4d_brsslic",
    lower = 0, upper = 16, integer = TRUE
  )
  global_passes <- .cluster4d_scalar_number(
    global_passes, "global_passes", "cluster4d_brsslic",
    lower = 0, upper = 16, integer = TRUE
  )
  seed <- .cluster4d_scalar_number(
    seed, "seed", "cluster4d_brsslic",
    lower = -.Machine$integer.max, upper = .Machine$integer.max, integer = TRUE
  )
  if (!is.null(min_size)) {
    min_size <- .cluster4d_scalar_number(
      min_size, "min_size", "cluster4d_brsslic", lower = 0, integer = TRUE
    )
  }

  n_threads <- 1L
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
  label_values <- as.integer(core$label_ids)
  if (!identical(label_values, seq_len(as.integer(core$actual_k))) ||
      !identical(sort(unique(labels)), label_values)) {
    stop("cluster4d_brsslic: C++ core returned a non-contiguous final label map")
  }

  grp <- factor(labels, levels = label_values)
  counts <- as.numeric(tabulate(as.integer(grp), nbins = length(label_values)))
  counts[counts == 0] <- NA_real_
  coord_sum <- rowsum(data_prep$coords, grp, reorder = TRUE)
  coord_centers <- coord_sum / counts
  dimnames(coord_centers) <- NULL

  centers <- unname(as.matrix(core$original_centers))
  centers_xyz_voxel <- unname(as.matrix(core$centers_xyz))

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
        seed = as.integer(seed),
        min_size = if (is.null(min_size)) NULL else as.integer(min_size)
      ),
      metadata = list(
        cpp_params = core$params,
        label_values = label_values,
        centers_xyz_voxel = centers_xyz_voxel,
        approximation = list(
          schema_version = "1.0.0",
          family = "boundary_refined_correlation_sketch_slic",
          coarse_embedding = list(
            basis = "hash",
            dimension = as.integer(embedding_dim),
            repeats = as.integer(sketch_repeats),
            centers = unname(as.matrix(core$coarse_embedding$centers)),
            coord_centers_voxel = unname(as.matrix(core$coarse_embedding$centers_xyz)),
            row_to_label = as.integer(core$coarse_embedding$label_ids),
            original_label_ids = as.integer(core$coarse_embedding$original_label_ids),
            counts = as.integer(core$coarse_embedding$counts),
            stage = "pre_refinement"
          ),
          refinement = list(
            correlation = "exact_pearson",
            boundary_passes = as.integer(boundary_passes),
            global_passes = as.integer(global_passes),
            temporal_stride = as.integer(core$params$refine_stride),
            spatial_weight = as.numeric(refine_spatial_weight),
            l2_weight = as.numeric(refine_l2_weight)
          ),
          native_final = list(
            original_centers = centers,
            coord_centers_voxel = centers_xyz_voxel,
            label_ids = as.integer(core$label_ids),
            original_label_ids = as.integer(core$original_label_ids),
            counts = as.integer(core$counts)
          )
        )
      )
    ),
    class = c("cluster4d_result", "cluster_result", "list")
  )

  finalize_cluster4d_result(result, vec, mask, "brs_slic", result$parameters)
}

#' Cluster4d using slice_msf method
#'
#' Wrapper for slice_msf algorithm with standardized interface.
#'
#' @inheritParams cluster4d
#' @param num_runs Number of independent runs
#' @param consensus Use consensus fusion
#' @param stitch_z Whether graph edges may cross axial slices.
#' @param r Number of non-DC DCT modes.
#' @param gamma Reserved compatibility parameter; must be zero.
#' @param z_mult Optional axial sketch-smoothing fraction in `[0, 1]`.
#' @param min_size Minimum component size after connectivity enforcement. Default NULL (auto).
#' @param seed Integer seed for multi-run DCT subspaces.
#' @param ensemble_fraction Candidate-pool fraction for multi-run DCT subspaces.
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
                               r = 12,
                               gamma = 0,
                               z_mult = 0,
                               min_size = NULL,
                               seed = 1L,
                               ensemble_fraction = 0.8,
                               ...) {

  if (!missing(parallel)) {
    stop("cluster4d_slice_msf: parallel is not supported", call. = FALSE)
  }
  input <- validate_cluster4d_inputs(
    vec, mask, n_clusters, "cluster4d_slice_msf"
  )

  # Convert spatial_weight to compactness
  # slice_msf uses medium compactness values (1-10 typical)
  compactness <- spatial_weight * 10

  # Map connectivity (slice_msf uses 2D connectivity)
  nbhd <- if (connectivity == 26) 8 else 4

  # Compute a conservative default min_size from expected cluster size.
  # Keep this low so exact-K merging can preserve boundaries on heterogeneous data.
  if (is.null(min_size)) {
    n_voxels <- length(input$mask_idx)
    expected_size <- n_voxels / n_clusters
    min_size <- max(2L, min(20L, floor(expected_size * 0.1)))
  }

  result <- slice_msf(
    vec = vec,
    mask = mask,
    target_k_global = n_clusters,
    compactness = compactness,
    min_size = min_size,
    nbhd = nbhd,
    num_runs = num_runs,
    consensus = consensus,
    use_features = num_runs > 1L,
    stitch_z = stitch_z,
    r = r,
    gamma = gamma,
    z_mult = z_mult,
    seed = seed,
    ensemble_fraction = ensemble_fraction,
    ...
  )
  
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
      r = r,
      gamma = gamma,
      z_mult = z_mult,
      seed = seed,
      ensemble_fraction = ensemble_fraction
    ),
    dots  # Include any additional parameters passed through
  )
  
  finalize_cluster4d_result(result, vec, mask, "slice_msf", result$parameters)
}

#' Cluster4d using FLASH-3D method
#'
#' Wrapper for FLASH-3D algorithm with standardized interface.
#'
#' @inheritParams cluster4d
#' @param lambda_t Optional finite non-negative temporal weight for Hamming
#'   distance. When omitted, it is \code{1 - spatial_weight}, so the unified
#'   interface preserves both feature-only and spatial-only endpoints.
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
                             lambda_t = NULL,
                             bits = 64,
                             dctM = 16,
                             ...) {

  vec <- ensure_neurovec(vec)
  validate_cluster4d_inputs(vec, mask, n_clusters, "cluster4d_flash3d")
  spatial_weight <- .cluster4d_scalar_number(
    spatial_weight, "spatial_weight", "cluster4d_flash3d"
  )
  if (spatial_weight < 0 || spatial_weight > 1) {
    stop("cluster4d_flash3d: spatial_weight must be between 0 and 1", call. = FALSE)
  }
  normalized_weight_mixture <- is.null(lambda_t)
  if (normalized_weight_mixture) lambda_t <- 1 - spatial_weight
  
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
  result$parameters <- .cluster4d_merge_parameters(
    result$parameters,
    c(list(
      n_clusters_requested = n_clusters,
      spatial_weight = spatial_weight,
      lambda_s = spatial_weight,
      lambda_t = lambda_t,
      normalized_weight_mixture = normalized_weight_mixture,
      max_iterations = max_iterations,
      rounds = max_iterations,
      bits = bits,
      dctM = dctM
    ), dots)
  )

  finalize_cluster4d_result(result, vec, mask, "flash3d", result$parameters)
}

#' Cluster4d using ACSC method
#'
#' Wrapper for Adaptive Correlation Superclustering algorithm with standardized interface.
#' ACSC uses graph-based clustering with Louvain community detection and optional
#' boundary refinement.
#'
#' @inheritParams cluster4d
#' @param block_size Approximate side length of initial blocks. Default 2.
#' @param correlation_metric Correlation definition used consistently for graph
#'   construction and refinement: `"pearson"`, `"spearman"`, or `"robust"`.
#' @param spatial_weighting Spatial block-edge weights: `"gaussian"` or `"binary"`.
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
                           correlation_metric = c("pearson", "spearman", "robust"),
                           spatial_weighting = c("gaussian", "binary"),
                           refine = TRUE,
                           max_refine_iter = NULL,
                           ...) {

  validate_cluster4d_inputs(vec, mask, n_clusters, "cluster4d_acsc")
  correlation_metric <- match.arg(correlation_metric)
  spatial_weighting <- match.arg(spatial_weighting)

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
    correlation_metric = correlation_metric,
    spatial_weighting = spatial_weighting,
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
  result$parameters <- .cluster4d_merge_parameters(
    result$parameters,
    c(list(
      n_clusters_requested = n_clusters,
      spatial_weight = spatial_weight,
      alpha = alpha,
      block_size = block_size,
      correlation_metric = correlation_metric,
      spatial_weighting = spatial_weighting,
      refine = refine,
      max_refine_iter = max_refine_iter
    ), dots)
  )

  finalize_cluster4d_result(result, vec, mask, "acsc", result$parameters)
}
