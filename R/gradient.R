# ThreeByThreeOffset <- rbind(c(1,0,0),
#                             c(-1,0,0),
#                             c(0,1,0),
#                             c(0,-1,0),
#                             c(0,0,1),
#                             c(0,0,-1))
#
#
#
# TwoByTwoOffset <- rbind(c(1,0,0),
#                         c(-1,0,0),
#                         c(0,1,0),
#                         c(0,-1,0))
#






#' Spatial Gradient Calculation
#'
#' The spatial_gradient function calculates the spatial gradient of a \code{NeuroVol} instance within the specified mask.
#'
#' @param vol A \code{NeuroVol} instance for which the spatial gradient should be calculated.
#' @param mask A \code{NeuroVol} mask defining the voxels to include in the spatial gradient calculation.
#' If the mask contains \code{numeric} data, nonzero values will define the included voxels.
#' If the mask is a \code{\linkS4class{LogicalNeuroVol}}, then \code{TRUE} will define the set
#' of included voxels.
#' @param sigma A numeric value controlling the spatial weighting function. Default is 0.5.
#'
#' @return A \code{NeuroVol} instance containing the spatial gradient values for the input \code{vol}.
#'
#' @examples
#' mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
#' input_vol <- NeuroVol(array(runif(202020), c(20,20,20)),
#' NeuroSpace(c(20,20,20)))
#'
#' gradient_vol <- spatial_gradient(input_vol, mask)
#'
#' @seealso \code{\link{spatial_laplacian}}, \code{\link{weighted_spatial_adjacency}}
#' @importFrom neighborweights spatial_adjacency spatial_laplacian
#' @importFrom neuroim2 NeuroVol
#' @import assertthat
#'
#' @export
spatial_gradient <- function(vol, mask, sigma=.5) {
  mask.idx <- which(mask>0)

  G <- neighborweights::spatial_adjacency(index_to_coord(mask, mask.idx),
                                          dthresh=9, nnk=9,
                                          weight_mode="heat",
                                          sigma=.8, stochastic=TRUE)
  v <- vol[mask.idx]
  vs <- G %*% vol[mask.idx]

  S = neighborweights::spatial_laplacian(index_to_coord(mask, mask.idx),
                                         weight_mode="heat",
                                         nnk=27,
                                         dthresh=6,
                                         sigma=sigma,
                                         normalize=FALSE,
                                         stochastic=FALSE)

  grad <- S %*% vs
  NeuroVol(as.vector(grad), space(mask), indices=mask.idx)


}



#' @keywords internal
#' @noRd
# Function to check if any neighboring 4 voxels have a lower gradient value and update seeds
perturb_seeds <- function(grad, seeds) {
  result <- logical(nrow(seeds))
  new_seeds <- seeds

  for (i in seq_along(result)) {
    seed_i <- seeds[i, 1]
    seed_j <- seeds[i, 2]
    seed_k <- seeds[i, 3]
    seed_gradient <- grad[seed_i, seed_j, seed_k]

    neighbors <- matrix(c(seed_i - 1, seed_j, seed_k,
                          seed_i + 1, seed_j, seed_k,
                          seed_i, seed_j - 1, seed_k,
                          seed_i, seed_j + 1, seed_k), ncol = 3, byrow = TRUE)

    for (n in 1:nrow(neighbors)) {
      ni <- neighbors[n, 1]
      nj <- neighbors[n, 2]
      nk <- neighbors[n, 3]

      if (ni >= 1 && nj >= 1 && nk >= 1 &&
          ni <= dim(grad)[1] && nj <= dim(grad)[2] && nk <= dim(grad)[3]) {

        if (grad[ni, nj, nk] < seed_gradient) {
          result[i] <- TRUE
          new_seeds[i, ] <- c(ni, nj, nk)
          break
        }
      }
    }
  }

  return(list(result = result, new_seeds = new_seeds))
}



