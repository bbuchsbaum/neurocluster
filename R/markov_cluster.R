

#' @importFrom neighborweights
markov_cluster <- function(mask, bvec, inflation=2, expansion=2, wsigma=.7, dsigma=2, ...) {
  mask.idx <- which(mask > 0)
  grid <- indexToCoord(mask, mask.idx)
  vgrid <- indexToGrid(mask, mask.idx)

  fmat <- series(bvec, mask.idx)

  A <- weighted_spatial_adjacency(grid, t(fmat), wsigma=wsigma, sigma=dsigma, weight_mode="heat",...)
  D <- colSums(A)
  An <- A %*% Diagonal(x=1/D)

  iter = 1
  while (iter < 10) {
    An <- An %*% An
    An <- An * An

    D <- colSums(An)
    conv = sum(D) - nrow(A)
    print(conv)
    An <- An %*% Diagonal(x=1/D)
    iter <- iter+1


  }
}
