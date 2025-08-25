
#' Fast 4D SLIC supervoxels (mask-aware, gradient relocation, preserve-K)
#'
#' @param bvec NeuroVec (X,Y,Z,T)
#' @param mask NeuroVol or logical 3D array
#' @param K target supervoxels
#' @param compactness SLIC 'm' tradeoff (higher = more compact)
#' @param max_iter outer SLIC iters
#' @param n_threads threads (0=auto)
#' @param step_mm seed spacing in mm (NULL=auto from bbox/K)
#' @param n_components JL random projection of time series (0 = no projection)
#' @param feature_norm 'zscale','l2','none'
#' @param seed_method 'mask_poisson','mask_grid','grid','farthest'
#' @param seed_relocate 'correlation','intensity','spatial','none'
#' @param seed_relocate_radius search radius (voxels) for relocation
#' @param connectivity 6 or 26
#' @param strict_connectivity enforce one connected component per label
#' @param enforce_connectivity alias for strict_connectivity (back-compat)
#' @param preserve_k ensure exactly K non-empty labels (post-connectivity split+refine)
#' @param topup_iters small number of refine iterations after splitting
#' @param min_size minimum component size for merges (currently not used in strict pass)
#' @param verbose logical
#' @export
slic4d_supervoxels <- function(bvec, mask,
                               K,
                               compactness = 10,
                               max_iter = 10,
                               n_threads = 0,
                               step_mm = NULL,
                               n_components = 0,
                               feature_norm = c("zscale", "l2", "none"),
                               seed_method = c("mask_poisson","mask_grid","grid","farthest"),
                               seed_relocate = c("correlation","intensity","spatial","none"),
                               seed_relocate_radius = 1L,
                               connectivity = c(26L,6L),
                               strict_connectivity = TRUE,
                               enforce_connectivity = TRUE,
                               preserve_k = FALSE,
                               topup_iters = 2L,
                               min_size = 0L,
                               verbose = FALSE) {
  feature_norm <- match.arg(feature_norm)
  seed_method  <- match.arg(seed_method)
  seed_relocate <- match.arg(seed_relocate)
  connectivity <- as.integer(match.arg(as.character(connectivity), c("26","6")))

  stopifnot(inherits(bvec, "NeuroVec"))
  if (inherits(mask, "NeuroVol")) {
    mask_arr <- as.logical(mask)
    sp <- neuroim2::space(mask)
  } else {
    mask_arr <- as.logical(mask)
    sp <- neuroim2::space(neuroim2::vol(mask_arr))
  }

  mask_idx <- which(mask_arr)
  if (length(mask_idx) == 0) stop("Mask is empty")

  dims <- dim(mask_arr)
  voxmm <- neuroim2::voxdim(sp)

  # coords in mm for masked voxels
  coords <- neuroim2::index_to_coord(sp, mask_idx)
  coords <- as.matrix(coords)

  # features (time series per voxel)
  feat <- neuroim2::series(bvec, mask_idx)
  feat <- as.matrix(feat)

  # normalization
  if (feature_norm == "zscale") {
    mu <- colMeans(feat)
    sdv <- apply(feat, 2, sd); sdv[sdv==0] <- 1
    feat <- sweep(sweep(feat, 2, mu, "-"), 2, sdv, "/")
  } else if (feature_norm == "l2") {
    nrm <- sqrt(rowSums(feat*feat)); nrm[nrm==0] <- 1
    feat <- feat / nrm
  }

  # optional random projection
  if (n_components > 0 && n_components < ncol(feat)) {
    set.seed(123)
    R <- matrix(rnorm(ncol(feat) * n_components) / sqrt(n_components), nrow = ncol(feat), ncol = n_components)
    feat <- feat %*% R
  }

  # default step in mm
  if (is.null(step_mm)) {
    mins <- apply(coords, 2, min); maxs <- apply(coords, 2, max)
    vol  <- prod(pmax(1e-6, maxs - mins))
    step_mm <- (vol / K)^(1/3)
  }

  # gradient volume for relocation
  grad_masked <- numeric(0)
  if (seed_relocate != "none") {
    if (seed_relocate == "correlation") {
      img4d <- as.array(bvec)
      grad3d <- correlation_gradient_cpp(img4d, as.numeric(mask_arr))
      dim(grad3d) <- dims
      grad_masked <- grad3d[mask_idx]
    } else if (seed_relocate == "intensity") {
      img4d <- as.array(bvec)
      mean3d <- apply(img4d, c(1,2,3), mean)
      grad3d <- .grad3d_fdiff(mean3d)
      grad_masked <- grad3d[mask_idx]
    } else if (seed_relocate == "spatial") {
      if (!requireNamespace("neighborweights", quietly = TRUE)) {
        warning("neighborweights not installed; falling back to intensity gradient")
        img4d <- as.array(bvec); mean3d <- apply(img4d, c(1,2,3), mean)
        grad3d <- .grad3d_fdiff(mean3d)
      } else {
        grad3d <- spatial_gradient(neuroim2::vol(apply(as.array(bvec), c(1,2,3), mean), space=sp),
                                   neuroim2::vol(mask_arr, space=sp))
        if (inherits(grad3d, "NeuroVol")) grad3d <- as.array(grad3d)
      }
      grad_masked <- grad3d[mask_idx]
    }
  }

  core <- slic4d_core(feat, coords,
                      mask_lin_idx = as.integer(mask_idx) - 1L,
                      dims = as.integer(dims),
                      voxmm = as.numeric(voxmm),
                      K = as.integer(K),
                      compactness = as.numeric(compactness),
                      max_iter = as.integer(max_iter),
                      step_mm = as.numeric(step_mm),
                      n_threads = as.integer(n_threads),
                      seed_method = seed_method,
                      enforce_connectivity = isTRUE(enforce_connectivity) || isTRUE(strict_connectivity),
                      min_size = as.integer(min_size),
                      connectivity = as.integer(connectivity),
                      strict_connectivity = isTRUE(strict_connectivity),
                      preserve_k = isTRUE(preserve_k),
                      topup_iters = as.integer(topup_iters),
                      grad_masked = as.numeric(grad_masked),
                      seed_relocate_radius = as.integer(seed_relocate_radius),
                      verbose = verbose)

  kvol <- ClusteredNeuroVol(mask_arr, clusters = core$labels)

  structure(list(
    clusvol       = kvol,
    cluster       = core$labels,
    centers       = core$center_feats,
    coord_centers = core$center_coords
  ), class = c("cluster_result", "list"))
}

.grad3d_fdiff <- function(arr) {
  dx <- arr * 0; dy <- dx; dz <- dx
  dx[-1,,] <- arr[-1,,] - arr[-nrow(arr),,]
  dy[,-1,] <- arr[,-1,] - arr[,-ncol(arr),]
  dz[,,-1] <- arr[,,-1] - arr[,,-dim(arr)[3]]
  sqrt(dx*dx + dy*dy + dz*dz)
}

#' Visualize gradient used for relocation
#' @export
slic4d_grad_summary <- function(bvec, mask, method=c("correlation","intensity","spatial")) {
  method <- match.arg(method)
  dims <- if (inherits(mask,"NeuroVol")) neuroim2::dim(mask) else dim(mask)
  sp <- if (inherits(mask,"NeuroVol")) neuroim2::space(mask) else neuroim2::space(neuroim2::vol(mask))
  if (method=="correlation") {
    img4d <- as.array(bvec)
    grad3d <- correlation_gradient_cpp(img4d, as.numeric(as.logical(mask)))
    dim(grad3d) <- dims
    return(grad3d)
  } else if (method=="intensity") {
    img4d <- as.array(bvec)
    mean3d <- apply(img4d, c(1,2,3), mean)
    return(.grad3d_fdiff(mean3d))
  } else {
    if (!requireNamespace("neighborweights", quietly = TRUE))
      stop("neighborweights not installed for 'spatial' method.")
    spatial_gradient(neuroim2::vol(apply(as.array(bvec), c(1,2,3), mean), space=sp),
                     neuroim2::vol(as.logical(mask), space=sp))
  }
}
