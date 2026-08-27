#' Evaluate a randomized decomposition under a fixed seed.
#'
#' Randomized SVD backends draw their projections from the session RNG, which
#' makes an otherwise deterministic clustering pipeline irreproducible. Pin the
#' seed for the call and put the caller's stream back exactly as it was.
#'
#' @keywords internal
#' @noRd
.svd_with_fixed_seed <- function(expr, seed = 20090601L) {
  .cluster4d_with_fixed_seed(expr, seed)
}


#' Compress High-Dimensional Features via Randomized SVD
#'
#' Reduces the dimensionality of fMRI time series data using Singular Value
#' Decomposition (SVD) while preserving the majority of variance. This is a key
#' component of the G3S (Gradient-Guided Geodesic Supervoxels) algorithm,
#' enabling 10-20x speedup in similarity computations.
#'
#' @param feature_mat A numeric matrix with voxels in rows and timepoints in columns
#'   (N x T). The function will center and scale the data automatically.
#' @param n_components Target number of dimensions (default: 15). Higher values
#'   preserve more variance but reduce speed gains.
#' @param variance_threshold Minimum proportion of variance to preserve (0-1).
#'   If the requested n_components doesn't meet this threshold, more components
#'   will be retained. Default: 0.95 (95% variance).
#' @param use_irlba Logical; if TRUE and the irlba package is available, use
#'   fast randomized SVD for large datasets (>10,000 voxels). Default: TRUE.
#' @param use_rsvd Logical; if TRUE and the rsvd package is available, prefer the
#'   randomized SVD implementation from \pkg{rsvd} for large datasets
#'   (>10,000 voxels). This can outperform irlba on tall-and-skinny matrices.
#'   Default: TRUE.
#'
#' @section Reproducibility:
#' Both randomized backends draw random projections. They are run under a fixed
#' internal seed and the caller's RNG state is restored afterwards, so repeated
#' calls on the same data return the same decomposition and the function never
#' perturbs the surrounding stream. Without that, two identical
#' \code{\link{cluster4d_g3s}} runs could disagree on half their voxel labels.
#'
#' @return A list with components:
#'   \item{features}{Compressed feature matrix (N x M), starting from the
#'     requested \code{n_components} and increasing M when needed to satisfy
#'     \code{variance_threshold}. Each row is normalized to unit length for
#'     cosine similarity calculations.}
#'   \item{variance_explained}{Proportion of total variance explained by the
#'     retained components (0-1).}
#'   \item{n_components}{Actual number of components retained (may be > n_components
#'     input if needed to meet variance_threshold).}
#'   \item{rotation}{The right singular vectors (V matrix) for transforming new data.}
#'   \item{singular_values}{The singular values (d vector) for each component.}
#'   \item{center}{Column means used for centering (for transforming new data).}
#'   \item{scale}{Column standard deviations used for scaling (for transforming new data).}
#'
#' @details
#' ## Algorithm
#'
#' The compression follows these steps:
#'
#' 1. **Center and scale**: Each timepoint is centered (mean=0) and scaled (sd=1)
#'    to ensure equal contribution from all timepoints.
#'
#' 2. **SVD decomposition**: Computes \code{X = U * D * t(V)} where:
#'    - U (N x k): Left singular vectors (voxel loadings)
#'    - D (k x k): Diagonal matrix of singular values
#'    - V (T x k): Right singular vectors (time loadings)
#'
#' 3. **Compression**: Retains only the first M components where M is chosen to
#'    balance speed (fewer components) and accuracy (more variance explained).
#'
#' 4. **Normalization**: Each compressed feature vector is normalized to unit length,
#'    enabling fast cosine similarity via dot products: \code{cor(x, y) approx x * y}
#'
#' ## Performance Characteristics
#'
#' - **Memory**: Compressed data uses M/T of original size (e.g., 15/300 = 5%)
#' - **Speed**: Similarity calculations are M/T times faster (e.g., 20x for 15/300)
#' - **Accuracy**: Typically preserves 95%+ of signal with M=10-20 for fMRI data
#'
#' ## Choosing n_components
#'
#' - **Small (M=10)**: Maximum speed, good for noisy or heavily smoothed data
#' - **Medium (M=15)**: Balanced, recommended default for most fMRI data
#' - **Large (M=25)**: High accuracy, use for high SNR data or critical applications
#'
#' The function will automatically increase n_components if needed to meet
#' variance_threshold. If that requires the full available rank, the exact
#' full-rank decomposition is used.
#'
#' @examples
#' # Simulate fMRI time series data
#' n_voxels <- 1000
#' n_timepoints <- 300
#' latent <- matrix(rnorm(n_voxels * 8), n_voxels, 8)
#' loadings <- matrix(rnorm(n_timepoints * 8), n_timepoints, 8)
#' feature_mat <- latent %*% t(loadings) +
#'   matrix(rnorm(n_voxels * n_timepoints, sd = 0.05), n_voxels, n_timepoints)
#'
#' # Compress to 15 dimensions (typical for G3S)
#' compressed <- compress_features_svd(feature_mat, n_components = 15)
#' print(compressed$variance_explained)  # At least 0.95
#' print(ncol(compressed$features))      # 15 here; may increase when necessary
#'
#' # More aggressive compression
#' compressed_fast <- compress_features_svd(feature_mat, n_components = 10)
#'
#' # Conservative compression (preserve 98% variance)
#' compressed_accurate <- compress_features_svd(
#'   feature_mat,
#'   n_components = 15,
#'   variance_threshold = 0.98
#' )
#'
#' @seealso
#' \code{\link{cluster4d_g3s}} for the full G3S clustering algorithm.
#' \code{\link{transform_new_data_svd}} for applying the compression to new data.
#'
#' @export
#' @importFrom stats sd
compress_features_svd <- function(feature_mat,
                                   n_components = 15,
                                   variance_threshold = 0.95,
                                   use_irlba = TRUE,
                                   use_rsvd = TRUE) {

  # Input validation
  if (!is.matrix(feature_mat)) {
    feature_mat <- as.matrix(feature_mat)
  }

  n_voxels <- nrow(feature_mat)
  n_timepoints <- ncol(feature_mat)
  max_possible <- min(n_voxels, n_timepoints)

  if (max_possible < 1L) {
    stop("feature_mat must have at least one row and one column")
  }
  if (length(n_components) != 1L || is.na(n_components) ||
      !is.finite(n_components) || n_components < 1 ||
      n_components != floor(n_components)) {
    stop("n_components must be a positive finite integer")
  }
  n_components <- as.integer(n_components)
  requested_components <- n_components

  if (n_components > max_possible) {
    warning("n_components (", n_components, ") exceeds matrix dimensions. ",
            "Setting to ", max_possible)
    n_components <- max_possible
    requested_components <- n_components
  }

  if (length(variance_threshold) != 1L || is.na(variance_threshold) ||
      !is.finite(variance_threshold) || variance_threshold < 0 ||
      variance_threshold > 1) {
    stop("variance_threshold must be between 0 and 1")
  }

  # Step 1: Center and scale.
  # Column statistics come from colMeans plus one pass over the centred matrix.
  # apply(feature_mat, 2, sd) splits the matrix into a list of columns and runs
  # an R-level loop over timepoints, which is the slowest step here on real
  # 4D data.
  observed <- is.finite(feature_mat)
  all_observed <- all(observed)
  n_observed <- if (all_observed) {
    rep.int(n_voxels, n_timepoints)
  } else {
    .colSums(observed, n_voxels, n_timepoints)
  }
  col_means <- colMeans(feature_mat, na.rm = TRUE)
  col_means[!is.finite(col_means)] <- 0

  feature_mat_scaled <- base::scale(feature_mat, center = col_means,
                                    scale = FALSE)
  # Centering a finite matrix by finite means cannot introduce a non-finite
  # value, so re-scanning the whole N x T copy is only needed when the input
  # itself had gaps.
  if (!all_observed && any(!is.finite(feature_mat_scaled))) {
    warning("Non-finite values detected after scaling. Replacing with 0.")
    feature_mat_scaled[!is.finite(feature_mat_scaled)] <- 0
  }
  col_sds <- sqrt(
    .colSums(feature_mat_scaled * feature_mat_scaled, n_voxels, n_timepoints) /
      pmax(n_observed - 1, 1)
  )
  col_sds[!is.finite(col_sds) | col_sds == 0] <- 1
  feature_mat_scaled <- feature_mat_scaled /
    rep(col_sds, each = n_voxels)
  attr(feature_mat_scaled, "scaled:center") <- col_means
  attr(feature_mat_scaled, "scaled:scale") <- col_sds

  randomized_min_voxels <- 10000L
  rsvd_available <- isTRUE(use_rsvd) && requireNamespace("rsvd", quietly = TRUE)
  use_rsvd_backend <- rsvd_available && n_voxels > randomized_min_voxels
  irlba_available <- isTRUE(use_irlba) && requireNamespace("irlba", quietly = TRUE)
  use_irlba_backend <- isTRUE(use_irlba) &&
    n_voxels > randomized_min_voxels && irlba_available

  # Gram route for tall-and-skinny data (4D imaging is always N >> T): the
  # T x T cross-product is tiny, its eigendecomposition yields the *entire*
  # spectrum in one shot, and the scores are one GEMM away. That is several
  # times faster than a full SVD of the N x T matrix, and because the whole
  # spectrum is visible the rank needed to satisfy variance_threshold is read
  # off directly instead of being searched for by repeated decomposition.
  # Skipped when T is large enough that the T x T product dominates, where the
  # randomized backends win instead.
  gram_max_dim <- 512L
  use_gram <- n_voxels >= n_timepoints && n_timepoints <= gram_max_dim

  compute_svd <- function(k, announce = TRUE) {
    if (use_rsvd_backend && k < max_possible) {
      if (announce) {
        message("Using randomized SVD (rsvd) with k=", k)
      }
      return(.svd_with_fixed_seed(
        rsvd::rsvd(feature_mat_scaled, k = k, nu = k, nv = k)
      ))
    }

    if (use_irlba_backend && k < max_possible) {
      if (announce) {
        message("Using fast randomized SVD (irlba) for large dataset")
      }
      return(.svd_with_fixed_seed(
        irlba::irlba(feature_mat_scaled, nu = k, nv = k)
      ))
    }

    svd_full <- base::svd(feature_mat_scaled, nu = k, nv = k)
    list(
      u = svd_full$u,
      d = svd_full$d[1:k],
      v = svd_full$v
    )
  }

  # Step 2: SVD decomposition and bounded rank expansion. Partial randomized
  # decompositions do not reveal the omitted spectrum, so keep increasing the
  # requested rank until the observed retained energy meets the contract. The
  # full-rank endpoint is computed with base SVD to make the terminal check
  # exact rather than dependent on a randomized backend's truncation behavior.
  threshold_tolerance <- 100 * .Machine$double.eps
  initial_variance_ratio <- NA_real_
  announce <- TRUE

  if (use_gram) {
    gram <- crossprod(feature_mat_scaled)
    spectrum <- eigen(gram, symmetric = TRUE)
    eigenvalues <- pmax(spectrum$values, 0)
    total_variance <- sum(eigenvalues)
    prefix <- if (total_variance == 0) {
      rep(1, length(eigenvalues))
    } else {
      cumsum(eigenvalues) / total_variance
    }
    feasible <- which(prefix + threshold_tolerance >= variance_threshold)
    retained <- if (length(feasible)) feasible[[1L]] else length(eigenvalues)
    initial_variance_ratio <- min(1, prefix[[min(n_components, length(prefix))]])
    n_components <- as.integer(min(max_possible, max(n_components, retained)))
    variance_ratio <- min(1, prefix[[n_components]])
    keep <- seq_len(n_components)
    rotation <- spectrum$vectors[, keep, drop = FALSE]
    # Scores are U * D; forming them directly avoids ever materialising U.
    compressed <- feature_mat_scaled %*% rotation
    svd_result <- list(
      u = NULL,
      d = sqrt(eigenvalues[keep]),
      v = rotation
    )
  } else {

    total_variance <- sum(feature_mat_scaled^2)

    repeat {
      svd_result <- compute_svd(n_components, announce = announce)
      announce <- FALSE
      explained_variance <- sum(svd_result$d^2)
      variance_ratio <- if (total_variance == 0) {
        1
      } else {
        min(1, max(0, explained_variance / total_variance))
      }
      if (is.na(initial_variance_ratio)) {
        initial_variance_ratio <- variance_ratio
      }

      if (variance_ratio + threshold_tolerance >= variance_threshold ||
          n_components >= max_possible) {
        break
      }

      ratio_estimate <- if (variance_ratio > 0) {
        min(max_possible, ceiling(n_components * variance_threshold / variance_ratio))
      } else {
        max_possible
      }
      n_components <- as.integer(min(
        max_possible,
        max(n_components + 10L, ceiling(1.5 * n_components), ratio_estimate)
      ))
    }

    # A growth step may deliberately overshoot the first feasible rank. Retain
    # the smallest available prefix that satisfies the threshold (but never less
    # than the user's requested dimensionality) so the accuracy contract does
    # not silently defeat the purpose of compression by retaining extra noise.
    if (total_variance > 0) {
      prefix_ratio <- cumsum(svd_result$d^2) / total_variance
      feasible <- which(prefix_ratio + threshold_tolerance >= variance_threshold)
      if (length(feasible)) {
        retained_components <- max(requested_components, feasible[[1L]])
        if (retained_components < n_components) {
          keep <- seq_len(retained_components)
          svd_result$u <- svd_result$u[, keep, drop = FALSE]
          svd_result$d <- svd_result$d[keep]
          svd_result$v <- svd_result$v[, keep, drop = FALSE]
          n_components <- as.integer(retained_components)
          variance_ratio <- min(1, prefix_ratio[[retained_components]])
        }
      }
    }

    # Step 4: compressed features = U * diag(d)
    compressed <- sweep(svd_result$u, 2, svd_result$d, `*`)

  }

  if (variance_ratio + threshold_tolerance < variance_threshold) {
    stop(
      "variance_threshold could not be met at the full available rank: retained ",
      signif(variance_ratio, 6), " < requested ", variance_threshold,
      call. = FALSE
    )
  }

  if (n_components > requested_components) {
    message(
      "Requested n_components (", requested_components, ") explained ",
      round(initial_variance_ratio * 100, 1), "% variance; increased to ",
      n_components, " components to retain ",
      round(variance_ratio * 100, 1), "%."
    )
  }

  # Step 5: Normalize each row to unit length for cosine similarity
  # For normalized vectors: cor(x, y) approx dot(x, y) when x and y are unit length
  row_norms <- sqrt(rowSums(compressed^2))
  row_norms[row_norms == 0] <- 1  # Avoid division by zero

  compressed_normalized <- compressed / row_norms

  # Return results
  list(
    features = compressed_normalized,
    variance_explained = variance_ratio,
    n_components = n_components,
    rotation = svd_result$v,
    singular_values = svd_result$d,
    center = col_means,
    scale = col_sds
  )
}


#' Transform New Data Using Existing SVD Compression
#'
#' Applies a previously computed SVD compression to new data, ensuring consistency
#' across multiple runs or when adding new voxels.
#'
#' @param new_data Numeric matrix (N_new x T) with same number of columns as
#'   the original training data.
#' @param compression_result List returned by \code{\link{compress_features_svd}}
#'   containing the rotation matrix and scaling parameters.
#'
#' @return Compressed and normalized feature matrix (N_new x M) where M is the
#'   number of components in the original compression.
#'
#' @examples
#' # Train on subset of data
#' train_data <- matrix(rnorm(500 * 300), 500, 300)
#' compression <- compress_features_svd(train_data, n_components = 15)
#'
#' # Apply to new data
#' test_data <- matrix(rnorm(100 * 300), 100, 300)
#' compressed_test <- transform_new_data_svd(test_data, compression)
#'
#' @export
transform_new_data_svd <- function(new_data, compression_result) {
  if (!is.matrix(new_data)) {
    new_data <- as.matrix(new_data)
  }

  # Apply same centering and scaling as training data
  new_data_scaled <- base::scale(new_data,
                                 center = compression_result$center,
                                 scale = compression_result$scale)

  # Handle non-finite values using the same policy as training compression.
  if (any(!is.finite(new_data_scaled))) {
    new_data_scaled[!is.finite(new_data_scaled)] <- 0
  }

  # Project onto principal components: X_new * V. On the training data this is
  # exactly U * D, which is the score representation returned above. Applying
  # D again would produce U * D^2 and change cosine geometry component-wise.
  compressed <- new_data_scaled %*% compression_result$rotation

  # Normalize to unit length
  row_norms <- sqrt(rowSums(compressed^2))
  row_norms[row_norms == 0] <- 1

  compressed_normalized <- compressed / row_norms

  compressed_normalized
}


#' Estimate Optimal Number of SVD Components
#'
#' Uses the elbow method or cumulative variance to suggest an appropriate
#' number of components for SVD compression.
#'
#' @param feature_mat Numeric matrix (N x T)
#' @param max_components Maximum components to test (default: 50)
#' @param method Either "elbow" (find elbow in scree plot) or "variance"
#'   (return number needed for 95% variance). Default: "variance".
#' @param variance_target Target variance for "variance" method (default: 0.95)
#'
#' @return Integer; suggested number of components
#'
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(1000 * 300), 1000, 300)
#' n_opt <- estimate_optimal_components(data)
#' print(n_opt)
#' }
#'
#' @export
estimate_optimal_components <- function(feature_mat,
                                       max_components = 50,
                                       method = c("variance", "elbow"),
                                       variance_target = 0.95) {
  method <- match.arg(method)

  if (!is.matrix(feature_mat)) {
    feature_mat <- as.matrix(feature_mat)
  }

  # Limit max_components to matrix dimensions
  max_components <- min(max_components, nrow(feature_mat), ncol(feature_mat))

  # Center and scale
  feature_mat_scaled <- base::scale(feature_mat, center = TRUE, scale = TRUE)
  feature_mat_scaled[is.na(feature_mat_scaled)] <- 0

  # Compute SVD
  if (requireNamespace("irlba", quietly = TRUE) && nrow(feature_mat) > 1000) {
    svd_result <- irlba::irlba(feature_mat_scaled, nu = 0, nv = max_components)
  } else {
    svd_result <- base::svd(feature_mat_scaled, nu = 0, nv = max_components)
  }

  # Variance explained by each component
  var_explained <- svd_result$d^2 / sum(feature_mat_scaled^2)
  cum_var <- cumsum(var_explained)

  if (method == "variance") {
    # Find first component that exceeds target
    n_opt <- which(cum_var >= variance_target)[1]
    if (is.na(n_opt)) {
      warning("Target variance ", variance_target, " not achievable with ",
              max_components, " components. Using all.")
      n_opt <- max_components
    }
  } else {
    # Elbow method: find point of maximum curvature
    # Use second derivative of variance explained
    if (length(var_explained) < 3) {
      n_opt <- length(var_explained)
    } else {
      # Compute normalized second derivative
      x <- seq_along(var_explained[1:(length(var_explained)-2)])
      y <- var_explained[1:(length(var_explained)-2)]
      second_deriv <- diff(diff(y))

      # Find maximum (most negative) second derivative
      n_opt <- which.min(second_deriv) + 1

      # Ensure reasonable bounds
      n_opt <- max(5, min(n_opt, max_components))
    }
  }

  as.integer(n_opt)
}
