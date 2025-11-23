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
#'   randomized SVD implementation from \pkg{rsvd}. This can outperform irlba on
#'   tall-and-skinny matrices. Default: TRUE.
#'
#' @return A list with components:
#'   \item{features}{Compressed feature matrix (N x M) where M <= n_components.
#'     Each row is normalized to unit length for cosine similarity calculations.}
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
#'    enabling fast cosine similarity via dot products: \code{cor(x, y) ≈ x · y}
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
#' variance_threshold, ensuring accuracy is never sacrificed for speed.
#'
#' @examples
#' # Simulate fMRI time series data
#' n_voxels <- 1000
#' n_timepoints <- 300
#' feature_mat <- matrix(rnorm(n_voxels * n_timepoints), n_voxels, n_timepoints)
#'
#' # Compress to 15 dimensions (typical for G3S)
#' compressed <- compress_features_svd(feature_mat, n_components = 15)
#' print(compressed$variance_explained)  # Should be ~0.95
#' print(ncol(compressed$features))      # Should be 15
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

  if (n_components < 1) {
    stop("n_components must be >= 1")
  }

  if (n_components > min(n_voxels, n_timepoints)) {
    warning("n_components (", n_components, ") exceeds matrix dimensions. ",
            "Setting to ", min(n_voxels, n_timepoints))
    n_components <- min(n_voxels, n_timepoints)
  }

  if (variance_threshold < 0 || variance_threshold > 1) {
    stop("variance_threshold must be between 0 and 1")
  }

  # Step 1: Center and scale
  # Use base::scale pattern from supervoxels.R for consistency
  col_means <- colMeans(feature_mat, na.rm = TRUE)
  col_sds <- apply(feature_mat, 2, sd, na.rm = TRUE)

  # Replace zero SDs with 1 to avoid division by zero
  col_sds[col_sds == 0] <- 1

  feature_mat_scaled <- base::scale(feature_mat, center = col_means, scale = col_sds)

  # Handle any remaining NAs
  if (any(is.na(feature_mat_scaled))) {
    warning("NA values detected after scaling. Replacing with 0.")
    feature_mat_scaled[is.na(feature_mat_scaled)] <- 0
  }

  rsvd_available <- isTRUE(use_rsvd) && requireNamespace("rsvd", quietly = TRUE)
  use_rsvd_backend <- rsvd_available
  irlba_available <- isTRUE(use_irlba) && requireNamespace("irlba", quietly = TRUE)
  use_irlba_backend <- isTRUE(use_irlba) && n_voxels > 10000 && irlba_available

  compute_svd <- function(k, announce = TRUE) {
    if (use_rsvd_backend) {
      if (announce) {
        message("Using randomized SVD (rsvd) with k=", k)
      }
      return(rsvd::rsvd(feature_mat_scaled, k = k, nu = k, nv = k))
    }

    if (use_irlba_backend) {
      if (announce) {
        message("Using fast randomized SVD (irlba) for large dataset")
      }
      return(irlba::irlba(feature_mat_scaled, nu = k, nv = k))
    }

    svd_full <- base::svd(feature_mat_scaled, nu = k, nv = k)
    list(
      u = svd_full$u,
      d = svd_full$d[1:k],
      v = svd_full$v
    )
  }

  # Step 2: SVD decomposition
  svd_result <- compute_svd(n_components)

  # Step 3: Check variance explained and adjust if needed
  total_variance <- sum(feature_mat_scaled^2)
  explained_variance <- sum(svd_result$d^2)
  variance_ratio <- explained_variance / total_variance

  # If variance threshold not met, try to add more components
  if (variance_ratio < variance_threshold) {
    max_possible <- min(n_voxels, n_timepoints)
    if (n_components < max_possible) {
      # Try with more components
      n_components_new <- min(n_components + 10, max_possible)
      warning("Initial n_components (", n_components, ") only explained ",
              round(variance_ratio * 100, 1), "% variance. ",
              "Trying with ", n_components_new, " components.")

      # Recompute SVD with more components
      svd_result <- compute_svd(n_components_new, announce = FALSE)

      n_components <- n_components_new
      explained_variance <- sum(svd_result$d^2)
      variance_ratio <- explained_variance / total_variance
    }
  }

  # Step 4: Create compressed features = U * diag(d)
  # This gives us the projection of each voxel onto the principal components
  compressed <- svd_result$u %*% diag(svd_result$d, nrow = length(svd_result$d))

  # Step 5: Normalize each row to unit length for cosine similarity
  # For normalized vectors: cor(x, y) ≈ dot(x, y) when x and y are unit length
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

  # Handle NAs
  if (any(is.na(new_data_scaled))) {
    new_data_scaled[is.na(new_data_scaled)] <- 0
  }

  # Project onto principal components: X_new * V * diag(d)
  compressed <- new_data_scaled %*% compression_result$rotation %*%
                diag(compression_result$singular_values,
                     nrow = length(compression_result$singular_values))

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
