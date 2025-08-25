# FLASH-3D: Fast 3D superclustering for fMRI
# Drop-in wrapper that accepts NeuroVec + mask (NeuroVol) like your previous method.
# It extracts masked time series, calls the C++ core, and returns a labeled volume.
#
# Dependencies: Rcpp, RcppParallel
# Compile once per session:
#   Rcpp::sourceCpp("flash3d.cpp")

#' FLASH-3D superclusters (native 3D, fMRI-aware, ultra-fast)
#'
#' @param vec A NeuroVec (4D) or a numeric 4D array [x,y,z,t].
#' @param mask A NeuroVol-like object or logical/int 3D array [x,y,z] where TRUE/nonzero marks voxels to cluster.
#' @param K Number of superclusters.
#' @param lambda_s Spatial weight (start). Will be annealed slightly upward over rounds.
#' @param lambda_t Temporal weight for Hamming distance of time-series hashes.
#' @param lambda_g Optional barrier weight (set 0 to disable). If >0, supply `barrier` volume.
#' @param rounds Number of outer rounds (seed→flood→recenter). 2–3 is plenty.
#' @param bits Length of temporal hash (64 or 128).
#' @param dctM Number of low DCT coefficients to rank-hash (default 12; 8–16 typical).
#' @param vox_scale Voxel size scaling for spatial distance (e.g., c(dx,dy,dz)); defaults to c(1,1,1).
#' @param barrier Optional 3D numeric array same dims as mask volume; higher values resist region growth.
#' @param verbose Print progress.
#' @return A list with components:
#'   \item{labels_vol}{integer 3D array of cluster ids (1..K) within mask, 0 outside).}
#'   \item{labels_mask}{integer vector of length sum(mask), the labels for masked voxels in mask order.}
#'   \item{K}{number of clusters}
#'   \item{method}{string, "FLASH-3D"}
#' @export
supervoxels_flash3d <- function(vec, mask, K,
                                lambda_s = 0.6,
                                lambda_t = 1.0,
                                lambda_g = 0.0,
                                rounds = 2L,
                                bits = 64L,
                                dctM = 12L,
                                vox_scale = c(1,1,1),
                                barrier = NULL,
                                verbose = TRUE) {
  dims <- infer_dims(vec, mask)
  nx <- dims[1]; ny <- dims[2]; nz <- dims[3]; T <- dims[4]
  # logical/int mask array [x,y,z]
  mask_arr <- as_mask_array(mask, c(nx,ny,nz))
  nmask <- sum(mask_arr != 0L)

  # Extract masked time series: T x Nmask
  ts <- extract_mask_timeseries(vec, mask_arr)
  # Linear indices of masked voxels, 1-based
  lin <- which(as.vector(mask_arr != 0L))

  # Optional barrier: full-grid numeric vector (length nx*ny*nz)
  barr_vec <- NULL
  if (!is.null(barrier)) {
    stopifnot(all(dim(barrier)[1:3] == c(nx,ny,nz)))
    barr_vec <- as.numeric(barrier)
  }

  labs_mask <- flash3d_supervoxels_cpp(ts, lin, c(nx,ny,nz), as.integer(K),
                                       c(lambda_s, lambda_t, lambda_g),
                                       rounds = as.integer(rounds),
                                       bits = as.integer(bits),
                                       dctM = as.integer(dctM),
                                       vox_scale = as.numeric(vox_scale),
                                       barrier_opt = barr_vec,
                                       verbose = isTRUE(verbose))

  # Build label volume
  labels_vol <- integer(nx*ny*nz)
  labels_vol[lin] <- labs_mask
  labels_vol <- array(labels_vol, dim = c(nx,ny,nz))

  structure(list(labels_vol = labels_vol,
                 labels_mask = labs_mask,
                 K = as.integer(K),
                 method = "FLASH-3D",
                 dims = c(nx,ny,nz,T)),
            class = c("ClusteredNeuroVol","FLASH3D","list"))
}

# --- Helpers to adapt to NeuroVec/NeuroVol or plain arrays --------------------

infer_dims <- function(vec, mask) {
  # vec: NeuroVec-like or 4D array [x,y,z,t]
  # mask: NeuroVol-like or 3D array [x,y,z]
  if (is.array(vec) && length(dim(vec)) == 4) {
    nx <- dim(vec)[1]; ny <- dim(vec)[2]; nz <- dim(vec)[3]; T <- dim(vec)[4]
  } else if (!is.null(attr(vec, "dim"))) {
    d <- attr(vec, "dim")
    if (length(d) == 4) { nx <- d[1]; ny <- d[2]; nz <- d[3]; T <- d[4] } else {
      stop("Unsupported vec type: need 4D array-like with dims [x,y,z,t].")
    }
  } else if (inherits(vec, "NeuroVec")) {
    # Best effort: try common accessors
    nx <- vec@dim[1]; ny <- vec@dim[2]; nz <- vec@dim[3]; T <- vec@dim[4]
  } else {
    stop("Unsupported vec type; pass a 4D array [x,y,z,t] or an object with @dim slot.")
  }
  c(nx,ny,nz,T)
}

as_mask_array <- function(mask, dims3) {
  if (is.logical(mask) || is.integer(mask) || is.numeric(mask)) {
    arr <- array(as.integer(mask != 0), dim = dims3)
  } else if (is.array(mask) && length(dim(mask))==3) {
    arr <- array(as.integer(mask != 0), dim = dim(mask))
  } else if (inherits(mask, "NeuroVol")) {
    # best effort; try to coerce via as.array
    arr <- as.integer(as.array(mask) != 0)
    arr <- array(arr, dim = dims3)
  } else {
    stop("Unsupported mask type; pass a 3D logical/int array or NeuroVol.")
  }
  arr
}

extract_mask_timeseries <- function(vec, mask_arr) {
  # Returns T x Nmask matrix
  if (is.array(vec) && length(dim(vec))==4) {
    nx <- dim(vec)[1]; ny <- dim(vec)[2]; nz <- dim(vec)[3]; T <- dim(vec)[4]
    inds <- which(as.vector(mask_arr != 0L))
    # gather
    # Convert 1D linear indices to 3D subscripts vectorized
    get_ts <- function(idx) {
      z <- (idx-1) %/% (nx*ny)
      rem <- (idx-1) %% (nx*ny)
      y <- rem %/% nx
      x <- rem %% nx
      vec[x+1, y+1, z+1, ]
    }
    # Build matrix column by column
    N <- length(inds)
    out <- matrix(0.0, nrow = T, ncol = N)
    for (j in seq_len(N)) out[,j] <- get_ts(inds[j])
    out
  } else if (inherits(vec, "NeuroVec")) {
    # Best effort: assume there is a values() accessor that returns a T x Nmask matrix
    if (!is.null(getGeneric("values"))) {
      values(vec, mask_arr != 0L) # may already be T x Nmask
    } else {
      stop("Don't know how to extract time series from NeuroVec; pass a 4D array instead.")
    }
  } else {
    stop("Unsupported vec type; pass a 4D array [x,y,z,t].")
  }
}
