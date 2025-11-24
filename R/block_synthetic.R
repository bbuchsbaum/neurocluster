#' Make a tiny block-structured synthetic volume
#'
#' A scikit-learn style "make_blobs" for neuroimaging demos. Builds a small
#' 2D/3D mask with three spatial bands (left/mid/right) and clearly separated
#' time-series patterns plus mild Gaussian noise. Designed to be fast,
#' deterministic, and visually intuitive in vignettes and tests.
#'
#' @param dims Integer vector length 3, spatial dimensions. Default `c(12,12,1)`.
#' @param ntime Number of time points. Default 80.
#' @param noise Standard deviation of added Gaussian noise. Default 0.15.
#' @param seed Optional integer seed for reproducibility. Default `NULL` (no change).
#'
#' @return List with components:
#' \describe{
#'   \item{vec}{\code{NeuroVec} 4D volume}
#'   \item{mask}{\code{NeuroVol} mask}
#'   \item{truth}{Integer ground-truth labels (length = prod(dims))}
#' }
#' @export
make_block_synthetic <- function(dims = c(12, 12, 1), ntime = 80, noise = 0.15, seed = NULL) {
  stopifnot(length(dims) == 3)
  if (!is.null(seed)) {
    old_seed <- .Random.seed
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }

  mask <- NeuroVol(array(1L, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  coords <- arrayInd(seq_len(nvox), dims)
  ts_mat <- matrix(0, nrow = nvox, ncol = ntime)

  left  <- coords[, 1] <= dims[1] / 3
  mid   <- coords[, 1] > dims[1] / 3 & coords[, 1] <= 2 * dims[1] / 3
  right <- coords[, 1] > 2 * dims[1] / 3

  t <- seq(0, 2 * pi, length.out = ntime)
  sig_left  <- sin(t)
  sig_mid   <- sign(sin(3 * t)) * 0.8
  sig_right <- sin(2 * t + pi / 4) * 0.6

  add_region <- function(sel, sig) {
    idx <- which(sel)
    ts_mat[idx, ] <<- sig + matrix(rnorm(length(idx) * ntime, sd = noise),
                                   nrow = length(idx))
  }
  add_region(left,  sig_left)
  add_region(mid,   sig_mid)
  add_region(right, sig_right)

  vols <- lapply(seq_len(ntime), function(tt) {
    v <- array(0, dims); v[mask > 0] <- ts_mat[, tt]; NeuroVol(v, NeuroSpace(dims))
  })

  truth <- integer(nvox)
  truth[left]  <- 1L
  truth[mid]   <- 2L
  truth[right] <- 3L

  list(
    vec = do.call(concat, vols),
    mask = mask,
    truth = truth,
    dims = dims
  )
}
