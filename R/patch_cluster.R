

patch_cluster <- function(bvec, mask, K=500, radius=8, knn=5, filter=list(lp=-1, hp=-1)) {

  mask.idx <- which(mask > 0)
  grid <- indexToCoord(mask, mask.idx)
  vgrid <- indexToGrid(mask, mask.idx)

  valmat <- series(bvec, mask.idx)

  if (any(filter) > 0) {
    message("turbo_cluster: filtering time series")
    valmat <- filter_mat(valmat, filter$lp, filter$hp)
  }

  slight <- neuroim::Searchlight(mask, radius=radius,eager=TRUE)


}
