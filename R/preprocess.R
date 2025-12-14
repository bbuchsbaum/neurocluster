#' Normalize Volumes by Removing Mean Offset
#'
#' Fast removal of volume-wise mean offset from 4D neuroimaging data.
#' Each volume/timepoint is centered to have zero mean across voxels.
#'
#' @param vec A \code{NeuroVec} object (4D neuroimaging data)
#' @param mask Optional \code{NeuroVol} mask. If NULL, uses all non-zero voxels.
#'
#' @return A \code{NeuroVec} with each volume centered to zero mean
#'
#' @details
#' This function removes global signal offset that can vary across timepoints
#' due to scanner drift, physiological noise, or other factors. For each

#' timepoint t:
#'
#' \deqn{output[,,,t] = input[,,,t] - mean(input[,,,t])}
#'
#' Uses fast C++ parallel implementation via RcppParallel.
#'
#' @examples
#' \dontrun{
#' vec <- neuroim2::read_vec("functional.nii.gz")
#' vec_normalized <- normalize_volumes(vec)
#' }
#'
#' @export
normalize_volumes <- function(vec, mask = NULL) {
  stopifnot(inherits(vec, "NeuroVec"))

  # Extract data matrix (voxels x time)
  if (is.null(mask)) {
    # Use all voxels
    mat <- neuroim2::as.matrix(vec)
  } else {
    # Use masked voxels
    idx <- which(mask > 0)
    mat <- neuroim2::series(vec, idx)
    mat <- t(mat)  # transpose to voxels x time
  }

  # Call C++ implementation
  mat_norm <- normalize_volumes_cpp(mat)

  # Reconstruct NeuroVec
  if (is.null(mask)) {
    # Full reconstruction
    dims <- dim(vec)
    arr <- array(0, dim = dims)
    arr[] <- as.vector(mat_norm)
    neuroim2::NeuroVec(arr, neuroim2::space(vec))
  } else {
    # Masked reconstruction
    idx <- which(mask > 0)
    dims <- dim(vec)
    arr <- array(0, dim = dims)
    for (t in seq_len(ncol(mat_norm))) {
      vol <- array(0, dim = dims[1:3])
      vol[idx] <- mat_norm[, t]
      arr[,,,t] <- vol
    }
    neuroim2::NeuroVec(arr, neuroim2::space(vec))
  }
}


#' Detrend Voxel Timeseries
#'
#' Fast removal of trends from each voxel's timeseries in 4D neuroimaging data.
#' Supports linear, polynomial, and DCT (discrete cosine transform) detrending.
#'
#' @param vec A \code{NeuroVec} object (4D neuroimaging data)
#' @param mask Optional \code{NeuroVol} mask. If NULL, uses all non-zero voxels.
#' @param method Detrending method: "linear" (default), "poly", or "dct"
#' @param degree For method="poly", the polynomial degree (1=linear, 2=quadratic, etc.)
#' @param n_basis For method="dct", number of DCT basis functions to remove.
#'   Higher values = more aggressive high-pass filtering.
#'
#' @return A \code{NeuroVec} with trend removed from each voxel
#'
#' @details
#' Three detrending methods are available:
#'
#' \strong{Linear (method="linear"):}
#' Fastest option. Removes linear drift by fitting y = a + b*t.
#'
#' \strong{Polynomial (method="poly"):}
#' Removes polynomial trend up to specified degree. degree=1 is linear,
#' degree=2 adds quadratic, etc. Uses orthonormalized basis for stability.
#'
#' \strong{DCT (method="dct"):}
#' Removes low-frequency components using discrete cosine transform basis.
#' This is a standard high-pass filter for fMRI data. The n_basis parameter
#' controls the cutoff:
#' \itemize{
#'   \item n_basis=1: mean removal only
#'   \item n_basis=2: approximately linear detrend
#'   \item n_basis=4-6: typical for fMRI (removes ~0.01 Hz drift for TR=2s)
#' }
#'
#' All methods use fast C++ parallel implementation via RcppParallel.
#'
#' @examples
#' \dontrun{
#' vec <- neuroim2::read_vec("functional.nii.gz")
#'
#' # Linear detrend (fastest)
#' vec_lin <- detrend_time(vec)
#'
#' # Quadratic detrend
#' vec_quad <- detrend_time(vec, method = "poly", degree = 2)
#'
#' # DCT high-pass filter (common for fMRI)
#' vec_dct <- detrend_time(vec, method = "dct", n_basis = 5)
#' }
#'
#' @export
detrend_time <- function(vec, mask = NULL, method = c("linear", "poly", "dct"),
                         degree = 2, n_basis = 4) {
  method <- match.arg(method)
  stopifnot(inherits(vec, "NeuroVec"))

  # Extract data matrix (voxels x time)
  if (is.null(mask)) {
    mat <- neuroim2::as.matrix(vec)
  } else {
    idx <- which(mask > 0)
    mat <- neuroim2::series(vec, idx)
    mat <- t(mat)
  }

  # Call appropriate C++ implementation

  mat_detrend <- switch(method,
    linear = detrend_time_cpp(mat),
    poly = detrend_poly_cpp(mat, degree),
    dct = detrend_dct_cpp(mat, n_basis)
  )

  # Reconstruct NeuroVec
  if (is.null(mask)) {
    dims <- dim(vec)
    arr <- array(0, dim = dims)
    arr[] <- as.vector(mat_detrend)
    neuroim2::NeuroVec(arr, neuroim2::space(vec))
  } else {
    idx <- which(mask > 0)
    dims <- dim(vec)
    arr <- array(0, dim = dims)
    for (t in seq_len(ncol(mat_detrend))) {
      vol <- array(0, dim = dims[1:3])
      vol[idx] <- mat_detrend[, t]
      arr[,,,t] <- vol
    }
    neuroim2::NeuroVec(arr, neuroim2::space(vec))
  }
}


#' Normalize and Detrend Neuroimaging Data
#'
#' Combined fast preprocessing: removes volume mean offsets and linear
#' trends in a single efficient operation.
#'
#' @param vec A \code{NeuroVec} object (4D neuroimaging data)
#' @param mask Optional \code{NeuroVol} mask. If NULL, uses all non-zero voxels.
#'
#' @return A \code{NeuroVec} with volume offsets and linear trends removed
#'
#' @details
#' Combines \code{\link{normalize_volumes}} and \code{\link{detrend_time}}
#' in a single pass for efficiency. The operations are applied in order:
#' 1. Remove volume-wise mean offset
#' 2. Remove linear trend from each voxel's timeseries
#'
#' Uses fast C++ parallel implementation via RcppParallel.
#'
#' @examples
#' \dontrun{
#' vec <- neuroim2::read_vec("functional.nii.gz")
#' vec_clean <- normalize_detrend(vec)
#' }
#'
#' @export
normalize_detrend <- function(vec, mask = NULL) {
  stopifnot(inherits(vec, "NeuroVec"))

  # Extract data matrix (voxels x time)
  if (is.null(mask)) {
    mat <- neuroim2::as.matrix(vec)
  } else {
    idx <- which(mask > 0)
    mat <- neuroim2::series(vec, idx)
    mat <- t(mat)
  }

  # Call combined C++ implementation
  mat_clean <- normalize_detrend_cpp(mat)

  # Reconstruct NeuroVec
  if (is.null(mask)) {
    dims <- dim(vec)
    arr <- array(0, dim = dims)
    arr[] <- as.vector(mat_clean)
    neuroim2::NeuroVec(arr, neuroim2::space(vec))
  } else {
    idx <- which(mask > 0)
    dims <- dim(vec)
    arr <- array(0, dim = dims)
    for (t in seq_len(ncol(mat_clean))) {
      vol <- array(0, dim = dims[1:3])
      vol[idx] <- mat_clean[, t]
      arr[,,,t] <- vol
    }
    neuroim2::NeuroVec(arr, neuroim2::space(vec))
  }
}
