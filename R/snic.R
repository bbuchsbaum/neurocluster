#' Find Initial Cluster Centers for Supervoxel Algorithm
#'
#' This function finds the initial cluster centers for a supervoxel algorithm.
#' Supervoxels are used to partition 3D image data into volumetric regions,
#' grouping similar voxels together. The initial cluster centers are crucial for
#' the performance and quality of the final supervoxels.
#'
#' @param cds A matrix or data frame representing the spatial coordinates of the voxels.
#' @param grad A vector representing the gradient values of the voxels.
#' @param K The desired number of supervoxels (clusters) in the output (default: 100).
#' @return A list containing two elements:
#'         selected - a vector of the selected indices corresponding to the initial cluster centers,
#'         coords - a matrix or data frame with the spatial coordinates of the initial cluster centers.
find_initial_points <- function(cds, grad, K=100) {
  # Validate inputs
  if (K > nrow(cds)) {
    warning("K larger than available points, using all available points")
    return(list(selected = 1:nrow(cds), coords = cds))
  }
  
  nfound=0
  batchsize=min(10, K)  # Don't exceed K
  sample_size=1000
  cds_cand <- cds
  sel <- c()
  iter <- 1
  sample_size <- min(sample_size, nrow(cds_cand))
  while (nfound < K) {

    cand <- sort(sample(1:nrow(cds_cand), sample_size))
    gvals <- grad[cds[cand,]]
    gvals <- (gvals - min(gvals))/ diff(range(gvals))
    if (iter == 1) {
      d <- RANN::nn2(cds_cand[cand,])
      ord <- order(d$nn.dists[,2] * (1-gvals), decreasing=TRUE)
      selected <- cand[ord[1:batchsize]]
    } else {

      d <- RANN::nn2(rbind(cds_cand[cand,], cds[sel,]), cds_cand[cand,])
      ord <- order(d$nn.dists[,2] * (1-gvals), decreasing=TRUE)

      selected <- cand[ord[1:batchsize]]
    }

    sel <- c(sel, selected)
    nfound <- length(sel)
    
    # Prevent selection of more points than needed
    if (nfound >= K) {
      sel <- sel[1:K]
      break
    }
    
    iter <- iter+1
  }

  list(selected=sel[1:min(K, length(sel))], coords=cds[sel[1:min(K, length(sel))],])

}


#' SNIC: Simple Non-Iterative Clustering
#'
#' The SNIC function performs a spatially constrained clustering on a \code{NeuroVec} instance
#' using the Simple Non-Iterative Clustering (SNIC) algorithm.
#'
#' @param vec A \code{NeuroVec} instance supplying the data to cluster.
#' @param mask A \code{NeuroVol} mask defining the voxels to include in the clustering result.
#' If the mask contains \code{numeric} data, nonzero values will define the included voxels.
#' If the mask is a \code{\linkS4class{LogicalNeuroVol}}, then \code{TRUE} will define the set
#' of included voxels.
#' @param compactness A numeric value controlling the compactness of the clusters, with larger values resulting
#' in more compact clusters. Default is 5.
#' @param K The number of clusters to find. Default is 500.
#' @param max_iter Maximum number of iterations for the SNIC algorithm. Default is 100.
#'   Currently ignored as SNIC algorithm uses internal convergence criteria.
#'
#' @return A \code{list} of class \code{snic_cluster_result} with the following elements:
#' \describe{
#' \item{clusvol}{An instance of type \linkS4class{ClusteredNeuroVol}.}
#' \item{gradvol}{A \code{NeuroVol} instance representing the spatial gradient of the reference volume.}
#' \item{cluster}{A vector of cluster indices equal to the number of voxels in the mask.}
#' \item{centers}{A matrix of cluster centers with each column representing the feature vector for a cluster.}
#' \item{coord_centers}{A matrix of spatial coordinates with each row corresponding to a cluster.}
#' }
#' 
#' @details
#' ## Parallelization Status
#' 
#' **Currently NOT parallelized.** SNIC uses a sequential priority queue-based algorithm
#' that processes voxels in order of their distance from cluster centers.
#' 
#' ### Sequential Operations:
#' 
#' 1. **Initialization**: Gradient-based seed selection using `find_initial_points()`
#'    - Finds K seed points with high gradient and spatial separation
#'    - Sequential search through candidate voxels
#' 
#' 2. **Priority Queue Processing** (C++ implementation):
#'    - Maintains a global priority queue of voxels to be assigned
#'    - Each voxel assignment depends on previously processed neighbors
#'    - Voxels are processed in order of their combined distance metric
#' 
#' 3. **Distance Computation**: For each voxel, calculates:
#'    - Feature distance to nearest cluster center
#'    - Spatial distance weighted by compactness parameter
#'    - Combined into single priority score
#' 
#' ### Why Not Parallelized:
#' 
#' - **Sequential dependency**: Priority queue enforces strict processing order
#' - **Global state**: Each voxel assignment affects subsequent assignments
#' - **Algorithm design**: SNIC's key innovation is its non-iterative, sequential nature
#' - **Coherent clusters**: Sequential processing ensures connected components
#' 
#' ### Performance Characteristics:
#' 
#' - **Complexity**: O(N log N) where N = number of voxels
#' - **Memory**: O(N) for priority queue and assignments
#' - **Speed**: Generally faster than iterative methods (supervoxels)
#' - **Single pass**: Processes each voxel exactly once
#' 
#' ### Performance Tips:
#' 
#' - **Reduce K**: Fewer clusters means less competition for voxels
#' - **Adjust compactness**: Higher values create more local clusters, faster processing
#' - **Pre-smooth data**: Reduce noise to improve gradient-based initialization
#' - **Use smaller masks**: Process ROIs separately if possible
#' - **Alternative**: Consider `slice_msf()` or `acsc()` for parallel execution
#' 
#' ### Comparison with Other Methods:
#' 
#' - **Faster than**: `supervoxels()` due to non-iterative nature
#' - **Slower than**: `slice_msf()` with parallel slices, `acsc()` with future backend
#' - **More coherent than**: Methods without spatial priority (ensures connectivity)
#'
#' @examples
#' mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
#' vec <- replicate(10, NeuroVol(array(runif(202020), c(20,20,20)),
#' NeuroSpace(c(20,20,20))), simplify=FALSE)
#' vec <- do.call(concat, vec)
#'
#' snic_res <- snic(vec, mask, compactness=5, K=100)
#'
#' @seealso \code{\link{supervoxels}}
#' @importFrom neuroim2 NeuroVec NeuroVol
#' @importFrom RANN nn2
#' @import assertthat
#'
#' @export
snic <- function(vec, mask, compactness=5, K=500, max_iter=100) {
  # Validate inputs
  if (K <= 0) {
    stop("K must be positive (K > 0), got: ", K)
  }
  
  mask.idx <- which(mask>0)
  valid_coords <- index_to_grid(mask, mask.idx)
  norm_coords <- sweep(valid_coords, 2, spacing(mask), "/")

  mask.grid <- index_to_grid(mask, mask.idx)
  mask_lookup <- array(0, dim(mask))
  mask_lookup[mask.grid] <- 0:(length(mask.idx)-1)

  vecmat <- series(vec, mask.idx)
  vecmat <- base::scale(vecmat)
  sam <- sample(1:ncol(vecmat), .05*ncol(vecmat))
  sf <- mean(as.vector(dist(t(vecmat[,sam])))^2)
  #vecmat <- vecmat/sf

  #ncd <- quantile(as.vector(dist(norm_coords[sam,])^2),.01)
  #norm_coords <- norm_coords/ncd


  refvol <- vols(vec)[[1]]
  grad <- spatial_gradient(refvol, mask)


  init <- find_initial_points(valid_coords, grad, K)
  centroid_idx = init$selected
  centroids <- norm_coords[init$selected,]


  s <- sqrt(nrow(valid_coords)/K)
  L <- array(0, dim(mask))

  ret = snic_main(L, mask@.Data,
                  centroids,
                  centroid_idx-1,
                  valid_coords-1,
                  norm_coords,
                  vecmat,
                  K, s,
                  compactness=compactness,
                  mask_lookup)


  # Create ClusteredNeuroVol with consistent logical mask (only positive values are TRUE)
  logical_mask <- mask > 0
  kvol <- ClusteredNeuroVol(logical_mask, clusters=ret[mask.idx])


  structure(
    list(clusvol=kvol,
         gradvol=grad,
         cluster=ret[mask.idx],
         centers=NULL,  # SNIC doesn't compute explicit centers
         coord_centers=NULL),  # SNIC doesn't compute explicit centers
    class=c("snic_cluster_result", "cluster_result", "list"))

}



