#' Generate a synthetic 4D neuroimaging volume with labeled clusters
#'
#' Convenient synthetic data for examples, vignettes, and tests. Each scenario
#' produces a `NeuroVec` (time series volume), a `NeuroVol` mask, and the
#' ground-truth cluster labels so you can benchmark algorithms or build
#' illustrations quickly.
#'
#' @param scenario Pattern to embed. One of `"gaussian_blobs"`, `"z_layers"`,
#'   or `"checkerboard"`.
#' @param dims Integer vector of length 3 giving the spatial grid dimensions.
#' @param n_clusters Number of latent clusters to embed.
#' @param n_time Number of time points.
#' @param spread Characteristic spread of each Gaussian blob (scalar or length 3).
#' @param noise_sd Standard deviation of additive Gaussian noise.
#' @param amplitude_range Range used when scaling each voxel's amplitude.
#' @param spacing_mm Physical voxel size supplied to the generated `NeuroSpace`.
#' @param pattern_type Type of temporal patterns to generate. One of
#'   `"orthogonal"` (QR-decomposition for maximally distinct signals, default),
#'   `"harmonic"` (phase/frequency-shifted sinusoids), or `"random"` (random
#'   linear combinations of harmonics).
#' @param seed Optional seed for reproducibility.
#'
#' @return A list with elements `vec`, `mask`, `truth`, `coords`,
#'   `patterns`, `weights`, `dims`, `n_clusters`, and `scenario`.
#' @examples
#' syn <- generate_synthetic_volume(
#'   scenario = "gaussian_blobs",
#'   dims = c(12, 12, 6),
#'   n_clusters = 4,
#'   seed = 1
#' )
#' str(syn)
#' truth_vol <- array(syn$truth, syn$dims)
#' clusvol <- neuroim2::ClusteredNeuroVol(syn$mask > 0, clusters = syn$truth)
#' plot(clusvol, slice = c(6, 6, 3), view = "axial")
#' @export
#' @importFrom neuroim2 NeuroVol NeuroVec NeuroSpace concat
generate_synthetic_volume <- function(
    scenario = c("gaussian_blobs", "z_layers", "checkerboard"),
    dims = c(16, 16, 8),
    n_clusters = 4,
    n_time = 30,
    spread = c(3.5, 3.5, 2.0),
    noise_sd = 0.05,
    amplitude_range = c(0.8, 1.2),
    spacing_mm = c(3, 3, 3),
    pattern_type = c("orthogonal", "harmonic", "random"),
    seed = NULL) {
  scenario <- match.arg(scenario)
  pattern_type <- match.arg(pattern_type)
  if (!is.null(seed)) set.seed(seed)
  coords <- as.matrix(expand.grid(x = seq_len(dims[1]),
                                  y = seq_len(dims[2]),
                                  z = seq_len(dims[3])))
  assignments <- switch(
    scenario,
    gaussian_blobs = assign_gaussian_clusters(coords, dims, n_clusters, spread),
    z_layers = assign_z_layers(coords, dims, n_clusters, jitter_sd = 0.3),
    checkerboard = assign_checkerboard(coords, n_clusters = n_clusters,
                                       block = max(2, floor(min(dims) / (2 * max(1, n_clusters)))))
  )
  truth <- assignments$labels
  weights <- assignments$weights
  nvox <- length(truth)
  patterns <- synthetic_time_patterns(n_clusters, n_time, pattern_types = pattern_type)
  signals <- matrix(0, nrow = nvox, ncol = n_time)
  amp <- runif(nvox, amplitude_range[1], amplitude_range[2]) * weights
  for (k in seq_len(n_clusters)) {
    idx <- which(truth == k)
    if (length(idx) == 0) next
    base_matrix <- matrix(rep(patterns[k, ], each = length(idx)), nrow = length(idx))
    signals[idx, ] <- amp[idx] * base_matrix
  }
  signals <- signals + matrix(rnorm(nvox * n_time, sd = noise_sd), nrow = nvox)
  space <- NeuroSpace(dims, spacing_mm)
  mask <- NeuroVol(array(TRUE, dims), space)
  vec_list <- lapply(seq_len(n_time), function(t) {
    NeuroVol(array(signals[, t], dim = dims), space)
  })
  vec <- do.call(concat, vec_list)
  list(
    vec = vec,
    mask = mask,
    truth = truth,
    coords = coords,
    patterns = patterns,
    weights = weights,
    dims = dims,
    n_clusters = n_clusters,
    scenario = scenario,
    pattern_type = pattern_type
  )
}

#' Generate synthetic temporal patterns for cluster signals
#'
#' Creates distinct temporal patterns for each cluster. Used internally by
#' \code{\link{generate_synthetic_volume}} but exported for custom synthetic
#' data generation.
#'
#' @param n_clusters Number of clusters/patterns to generate.
#' @param n_time Number of time points.
#' @param pattern_types Type of patterns: \code{"orthogonal"} (QR-decomposition
#'   for maximally distinct signals), \code{"harmonic"} (phase/frequency-shifted
#'   sinusoids), or \code{"random"} (random linear combinations of harmonics).
#'
#' @return A matrix of dimension \code{n_clusters x n_time} with one temporal
#'   pattern per row.
#' @export
synthetic_time_patterns <- function(n_clusters, n_time,
                                    pattern_types = c("orthogonal", "harmonic", "random")) {

  pattern_types <- match.arg(pattern_types)
  t_seq <- seq(0, 2 * pi, length.out = n_time)
  patterns <- matrix(0, nrow = n_clusters, ncol = n_time)

  if (pattern_types == "orthogonal") {
    # Use QR decomposition to create truly orthogonal signals
    # This ensures each cluster has a maximally distinct temporal pattern
    random_mat <- matrix(rnorm(n_time * n_clusters), nrow = n_time, ncol = n_clusters)
    qr_decomp <- qr(random_mat)
    Q <- qr.Q(qr_decomp)  # n_time x n_clusters orthonormal matrix
    # Scale to reasonable amplitude and transpose to n_clusters x n_time
    patterns <- t(Q) * sqrt(n_time)
  } else {
    # Legacy harmonic/random patterns (kept for backward compatibility)
    for (k in seq_len(n_clusters)) {
      type <- if (pattern_types == "harmonic") "harmonic" else "random"
      if (type == "random") {
        coeffs <- rnorm(3)
        patterns[k, ] <- coeffs[1] * sin(t_seq) + coeffs[2] * cos(2 * t_seq) + coeffs[3]
      } else {
        freq <- runif(1, 0.6, 1.6)
        phase <- runif(1, 0, 2 * pi)
        patterns[k, ] <- sin(freq * t_seq + phase)
      }
      patterns[k, ] <- patterns[k, ] / max(abs(patterns[k, ]))
    }
  }
  patterns
}

assign_gaussian_clusters <- function(coords, dims, n_clusters, spread) {
  nvox <- nrow(coords)
  if (length(spread) == 1) {
    spread <- rep(spread, 3)
  }
  spread_mat <- matrix(rep(spread, length.out = n_clusters * 3), nrow = n_clusters, ncol = 3, byrow = TRUE)
  pad <- pmax(matrix(2, n_clusters, 3), ceiling(spread_mat))
  centers <- matrix(0, nrow = n_clusters, ncol = 3)
  for (k in seq_len(n_clusters)) {
    centers[k, ] <- runif(3, min = pad[k, ], max = dims - pad[k, ] + 1)
  }
  d2 <- matrix(0, nrow = nvox, ncol = n_clusters)
  for (k in seq_len(n_clusters)) {
    diff <- sweep(coords, 2, centers[k, ], FUN = "-")
    scaled <- sweep(diff, 2, spread_mat[k, ], FUN = "/")
    d2[, k] <- rowSums(scaled^2)
  }
  weights <- exp(-0.5 * d2)
  labels <- max.col(weights, ties.method = "first")
  peak <- weights[cbind(seq_len(nvox), labels)]
  peak <- peak / max(peak)
  peak <- pmax(peak, 0.2)
  list(labels = labels, weights = peak, centers = centers)
}

assign_z_layers <- function(coords, dims, n_clusters, jitter_sd = 0.1) {
  z_vals <- coords[, 3] + rnorm(nrow(coords), sd = jitter_sd)
  breaks <- seq(0.5, dims[3] + 0.5, length.out = n_clusters + 1)
  labels <- cut(z_vals, breaks = breaks, labels = FALSE, include.lowest = TRUE)
  weights <- rep(1, length(labels))
  list(labels = labels, weights = weights)
}

assign_checkerboard <- function(coords, n_clusters, block = 2) {
  block_ids <- ((coords[, 1] - 1) %/% block + (coords[, 2] - 1) %/% block + (coords[, 3] - 1) %/% block) %% n_clusters
  labels <- block_ids + 1
  weights <- rep(1, length(labels))
  list(labels = labels, weights = weights)
}
