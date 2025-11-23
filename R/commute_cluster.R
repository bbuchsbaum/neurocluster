#' Handle Zero-Variance Voxels
#'
#' Internal helper to replace zero-variance or NA voxels with small noise.
#' Uses a local seed to avoid affecting the global RNG state.
#'
#' @param X Numeric matrix (timepoints x voxels)
#' @param seed Optional seed for reproducibility. If NULL, uses non-deterministic noise.
#' @return Cleaned matrix
#' @keywords internal
#' @noRd
handle_bad_voxels <- function(X, seed = NULL) {
  csds <- matrixStats::colSds(X)
  bad <- which(is.na(csds) | csds < 1e-9)

  if (length(bad) > 0) {
    warning(sprintf(
      "commute_cluster: %d voxels have zero variance or NA. Injecting small noise to avoid singularity.",
      length(bad)
    ))

    # Use seed if provided for reproducibility
    if (!is.null(seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
        get(".Random.seed", envir = .GlobalEnv)
      } else {
        NULL
      }
      on.exit({
        if (!is.null(old_seed)) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(seed)
    }

    X[, bad] <- matrix(rnorm(length(bad) * nrow(X), sd = 0.01), nrow(X), length(bad))
  }

  X
}


#' Commute Time Clustering
#'
#' Performs spatially constrained clustering on a \code{NeuroVec} instance
#' using commute time distance (spectral embedding) and K-means clustering.
#'
#' @param bvec A \code{NeuroVec} instance supplying the data to cluster.
#' @param mask A \code{NeuroVol} mask defining the voxels to include in the clustering result.
#'   If the mask contains \code{numeric} data, nonzero values will define the included voxels.
#'   If the mask is a \code{\linkS4class{LogicalNeuroVol}}, then \code{TRUE} will define the set
#'   of included voxels.
#' @param K The number of clusters to find. Default is 100.
#' @param ncomp The number of components to use for the commute time embedding.
#'   Default is the ceiling of \code{sqrt(K*2)}.
#' @param alpha A numeric value controlling the balance between spatial and feature similarity.
#'   Default is 0.5 (balanced). Range: 0 (spatial only) to 1 (feature only).
#' @param sigma1 A numeric value controlling the spatial weighting function. Default is 0.73.
#' @param sigma2 A numeric value controlling the feature weighting function. Default is 5.
#' @param connectivity An integer representing the number of nearest neighbors to consider
#'   when constructing the similarity graph. Default is 27 (full 3D neighborhood).
#' @param weight_mode A character string indicating the type of weight function for the
#'   similarity graph. Options are "binary" and "heat". Default is "heat".
#' @param noise_seed Optional integer seed for reproducible noise injection when handling
#'   zero-variance voxels. If NULL (default), uses non-deterministic noise.
#'   **Note**: This only controls noise injection, not k-means initialization. For full
#'   reproducibility, wrap the entire call in `set.seed()`.
#' @param verbose Logical. If TRUE, print progress messages. Default is TRUE.
#'
#' @return A \code{list} of class \code{commute_time_cluster_result} (inheriting from
#'   \code{cluster_result}) with the following elements:
#' \describe{
#'   \item{clusvol}{An instance of type \linkS4class{ClusteredNeuroVol}.}
#'   \item{cluster}{A vector of cluster indices equal to the number of voxels in the mask.}
#'   \item{centers}{A matrix of cluster centers (K x T) where T is the number of timepoints.}
#'   \item{coord_centers}{A matrix of spatial coordinates (K x 3) with each row
#'     corresponding to a cluster centroid.}
#'   \item{embedding}{The spectral embedding coordinates (N x ncomp) from commute time distance.}
#'   \item{n_clusters}{The number of clusters (same as K).}
#'   \item{method}{Character string "commute_time".}
#' }
#'
#' @details
#' ## Algorithm Overview
#'
#' Commute time clustering uses spectral graph theory to embed voxels into a lower-dimensional
#' space where geodesic distances on the graph approximate commute times (expected random walk
#' return times). This embedding respects both spatial proximity and feature similarity.
#'
#' The algorithm has three main steps:
#'
#' 1. **Graph Construction**: Build a weighted adjacency matrix combining spatial and
#'    feature similarity using `neighborweights::weighted_spatial_adjacency()`.
#'
#' 2. **Spectral Embedding**: Compute commute time distances via eigendecomposition of
#'    the graph Laplacian using `neighborweights::commute_time_distance()`.
#'
#' 3. **Clustering**: Apply k-means to the embedded coordinates.
#'
#' ## Scalability Warning
#'
#' **This method is computationally expensive and NOT recommended for whole-brain clustering.**
#'
#' - **Complexity**: O(N³) for eigendecomposition where N = number of voxels
#' - **Memory**: O(N²) for adjacency matrix storage
#' - **Practical limit**: ~10,000 voxels (ROI-based analysis)
#' - **Whole-brain**: ~100,000+ voxels will likely crash or take hours
#'
#' For large-scale clustering, consider:
#' - `slice_msf()`: Slice-based minimum spanning forests
#' - `acsc()`: Adaptive correlation superclustering
#' - `supervoxels()` or `snic()`: Iterative spatial methods
#'
#' ## Parallelization Status
#'
#' **Currently NOT explicitly parallelized.** The algorithm runs sequentially,
#' but matrix operations may use multi-threaded BLAS/LAPACK libraries.
#'
#' ### Why Not Parallelized:
#'
#' - **External dependencies**: Uses `neighborweights` package functions
#' - **Eigendecomposition**: Difficult to parallelize efficiently in R
#' - **Already optimized**: BLAS/LAPACK typically use multiple threads automatically
#' - **Bottleneck**: Eigendecomposition dominates runtime regardless of parallelization
#'
#' ### Performance Tips:
#'
#' - **Use optimized BLAS**: OpenBLAS, Intel MKL, or Apple Accelerate
#' - **Reduce connectivity**: Smaller neighborhoods = sparser matrices (e.g., 6 or 18)
#' - **Increase alpha**: Higher values emphasize features over space, reducing graph density
#' - **Use fewer components**: Set `ncomp` lower (e.g., `ncomp = K/2`) for faster embedding
#' - **Pre-filter voxels**: Remove low-variance voxels before clustering
#' - **ROI analysis**: Apply to small regions of interest rather than whole brain
#'
#' ### Common Issues:
#'
#' **Eigenvalue Errors**: Often due to singular or near-singular weight matrices.
#'
#' Causes:
#' - Duplicate or perfectly correlated time series
#' - Disconnected graph components
#' - Insufficient connectivity parameter
#'
#' Solutions:
#' - Increase `connectivity` (e.g., 27 instead of 6)
#' - Adjust `alpha` to balance spatial/feature weights
#' - Use `noise_seed` for reproducible noise injection
#' - Check for and remove constant voxels beforehand
#'
#' **Memory Errors**: Adjacency matrix requires O(N²) memory.
#'
#' Solutions:
#' - Reduce number of voxels (subsample or use smaller ROI)
#' - Use alternative method for large N
#'
#' @examples
#' \dontrun{
#' # Small example with synthetic data
#' library(neuroim2)
#' mask <- NeuroVol(array(1, c(20, 20, 20)), NeuroSpace(c(20, 20, 20)))
#' vec <- replicate(10, NeuroVol(array(runif(20*20*20), c(20, 20, 20)),
#'   NeuroSpace(c(20, 20, 20))), simplify = FALSE)
#' vec <- do.call(concat, vec)
#'
#' # Run clustering (8000 voxels - feasible for this method)
#' commute_res <- commute_cluster(vec, mask, K = 50, verbose = TRUE)
#'
#' # Access results
#' print(commute_res$n_clusters)
#' plot(commute_res$clusvol)
#' }
#'
#' \dontrun{
#' # With reproducible noise injection (for zero-variance voxels)
#' commute_res <- commute_cluster(vec, mask, K = 50, noise_seed = 42)
#'
#' # For full reproducibility (including k-means), use set.seed() wrapper
#' set.seed(123)
#' commute_res <- commute_cluster(vec, mask, K = 50, noise_seed = 42)
#'
#' # ROI-based analysis (recommended workflow)
#' roi_mask <- mask  # In practice, use a smaller ROI
#' roi_mask[1:10, , ] <- 0  # Reduce voxels
#' commute_res <- commute_cluster(vec, roi_mask, K = 30)
#' }
#' @seealso
#' \code{\link{snic}} for a faster non-iterative method
#' \code{\link{slice_msf}} for scalable slice-based clustering
#' \code{\link{acsc}} for large-scale correlation-based clustering
#'
#' @importFrom neuroim2 NeuroVec NeuroVol series index_to_coord ClusteredNeuroVol
#' @importFrom matrixStats colSds
#' @importFrom neighborweights weighted_spatial_adjacency commute_time_distance
#' @importFrom stats kmeans
#' @importFrom Matrix t
#' @import assertthat
#'
#' @export
commute_cluster <- function(bvec,
                            mask,
                            K = 100,
                            ncomp = ceiling(sqrt(K * 2)),
                            alpha = 0.5,
                            sigma1 = 0.73,
                            sigma2 = 5,
                            connectivity = 27,
                            weight_mode = c("binary", "heat"),
                            noise_seed = NULL,
                            verbose = TRUE) {

  # Validate inputs
  weight_mode <- match.arg(weight_mode)
  assertthat::assert_that(K > 0, msg = "K must be positive")
  assertthat::assert_that(ncomp > 0, msg = "ncomp must be positive")
  assertthat::assert_that(alpha >= 0 && alpha <= 1, msg = "alpha must be in [0, 1]")
  assertthat::assert_that(connectivity > 0, msg = "connectivity must be positive")

  # Extract mask indices and coordinates
  mask.idx <- which(mask > 0)
  n_voxels <- length(mask.idx)

  if (n_voxels == 0) {
    stop("Mask is empty. No voxels to cluster.")
  }

  # Scalability warning
  if (n_voxels > 15000) {
    warning(sprintf(
      paste0(
        "commute_cluster: Processing %d voxels. This method has O(N^3) complexity and O(N^2) memory.\n",
        "  Expected runtime: Very slow (possibly hours)\n",
        "  Memory required: ~%.1f GB\n",
        "  Consider using slice_msf(), acsc(), or supervoxels() for large datasets."
      ),
      n_voxels,
      (n_voxels^2 * 8) / 1e9  # Rough memory estimate in GB
    ))
  }

  grid <- index_to_coord(mask, mask.idx)

  if (verbose) {
    message(sprintf("commute_cluster: Extracted %d voxels from mask", n_voxels))
  }

  # Extract feature matrix (timepoints x voxels)
  feature_mat <- neuroim2::series(bvec, mask.idx)

  # Handle zero-variance voxels
  feature_mat <- handle_bad_voxels(feature_mat, seed = noise_seed)

  # Scale features
  feature_mat <- base::scale(feature_mat)

  # Build adjacency graph
  if (verbose) {
    message("commute_cluster: Constructing weighted spatial adjacency graph...")
  }

  W <- neighborweights::weighted_spatial_adjacency(
    grid,
    t(feature_mat),  # neighborweights expects voxels x features
    dthresh = sigma2 * 6,
    nnk = connectivity,
    wsigma = sigma1,
    sigma = sigma2,
    alpha = alpha,
    weight_mode = weight_mode,
    include_diagonal = FALSE,
    stochastic = TRUE
  )

  # Ensure symmetry (only if needed to avoid unnecessary computation)
  if (!isSymmetric(W)) {
    if (verbose) {
      message("commute_cluster: Symmetrizing adjacency matrix...")
    }
    # Handle sparse matrices efficiently
    if (inherits(W, "Matrix")) {
      W <- (W + Matrix::t(W)) / 2
    } else {
      W <- (W + t(W)) / 2
    }
  }

  # Commute time embedding (spectral decomposition)
  if (verbose) {
    message(sprintf(
      "commute_cluster: Computing commute-time embedding (%d components)...",
      ncomp
    ))
  }

  ct <- tryCatch({
    neighborweights::commute_time_distance(W, ncomp = ncomp)
  }, error = function(e) {
    if (grepl("TridiagEigen|eigen", e$message, ignore.case = TRUE)) {
      stop(
        "Eigenvalue decomposition failed in commute clustering.\n",
        "  Common causes:\n",
        "  1. Disconnected graph (try increasing 'connectivity')\n",
        "  2. Singular weight matrix (try adjusting 'alpha')\n",
        "  3. Perfectly correlated time series (use 'noise_seed')\n",
        "  4. Too few neighbors relative to ncomp\n",
        "\n  Try:\n",
        "  - Increase connectivity (e.g., connectivity = 27)\n",
        "  - Reduce ncomp (e.g., ncomp = K/2)\n",
        "  - Adjust alpha (try 0.3 or 0.7)\n",
        "  - Set noise_seed for reproducible noise injection\n",
        "\n  Original error: ", e$message,
        call. = FALSE
      )
    } else {
      stop("Error in commute time distance calculation: ", e$message, call. = FALSE)
    }
  })

  # K-means clustering on embedded coordinates
  if (verbose) {
    message(sprintf("commute_cluster: Running k-means with K=%d clusters...", K))
  }

  # Use multiple starts for robustness
  # Note: If noise_seed is set, k-means will also be reproducible
  # But we don't re-set the seed here to avoid affecting global RNG state
  # User should wrap entire call in set.seed() if full reproducibility needed
  kres <- kmeans(ct$cds, centers = K, iter.max = 500, nstart = 10)

  # Create clustered volume
  logical_mask <- mask > 0
  kvol <- ClusteredNeuroVol(logical_mask, kres$cluster)

  # Compute centroids in feature and spatial space
  if (verbose) {
    message("commute_cluster: Computing final centroids...")
  }

  centroids <- compute_centroids(feature_mat, grid, kres$cluster, medoid = FALSE)

  # Construct result object
  ret <- structure(
    list(
      clusvol = kvol,
      cluster = kres$cluster,
      centers = centroids$center,
      coord_centers = centroids$centroid,
      embedding = ct$cds,  # Include for potential downstream analysis
      n_clusters = K,
      method = "commute_time"
    ),
    class = c("commute_time_cluster_result", "cluster_result", "list")
  )

  if (verbose) {
    message("commute_cluster: Done!")
  }

  ret
}
