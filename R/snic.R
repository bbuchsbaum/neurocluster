# Load required libraries
library(collections)
library(Rcpp)

# cppFunction('
#   List update_centroid_online_rcpp(NumericVector centroid_x, NumericVector centroid_c, int centroid_n, NumericVector x_i, NumericVector c_i) {
#     // Calculate the new spatial position
#     NumericVector new_x = (centroid_x * centroid_n + x_i) / (centroid_n + 1);
#
#     // Calculate the new color value
#     NumericVector new_c = (centroid_c * centroid_n + c_i) / (centroid_n + 1);
#
#     // Update the number of pixels in the superpixel
#     int new_n = centroid_n + 1;
#
#     // Return the updated centroid as a list
#     return List::create(Named("x") = new_x, Named("c") = new_c, Named("n") = new_n);
#   }
# ')
#
#
#
# cppFunction('
#   List get_26_connected_neighbors_rcpp(int i, int j, int k, IntegerVector array_dims) {
#     List neighbors;
#
#     int max_i = array_dims[0];
#     int max_j = array_dims[1];
#     int max_k = array_dims[2];
#
#     for (int i_shift = -1; i_shift <= 1; ++i_shift) {
#       for (int j_shift = -1; j_shift <= 1; ++j_shift) {
#         for (int k_shift = -1; k_shift <= 1; ++k_shift) {
#           // Skip the center point (0, 0, 0) shift
#           if (i_shift == 0 && j_shift == 0 && k_shift == 0) {
#             continue;
#           }
#
#           int i_neighbor = i + i_shift;
#           int j_neighbor = j + j_shift;
#           int k_neighbor = k + k_shift;
#
#           // Check if the neighbor is within the array bounds
#           if (i_neighbor >= 1 && i_neighbor <= max_i &&
#               j_neighbor >= 1 && j_neighbor <= max_j &&
#               k_neighbor >= 1 && k_neighbor <= max_k) {
#             IntegerVector neighbor_coords = IntegerVector::create(i_neighbor, j_neighbor, k_neighbor);
#             neighbors.push_back(neighbor_coords);
#           }
#         }
#       }
#     }
#
#     return neighbors;
#   }
# ')
#
#
# cppFunction('
#   double compute_distance_rcpp(NumericVector x_i, NumericVector x_k,
#                                NumericVector c_i, NumericVector c_k,
#                                double s = 5, double compactness = 5) {
#     double dc = sum(pow(c_i - c_k, 2)); // squared Euclidean color distance
#     double ds = sum(pow(x_i - x_k, 2)); // squared Euclidean spatial distance
#
#     // Incorporate compactness parameter into the distance calculation
#     //double distance = sqrt(pow(dc, 2) + pow(ds / s * compactness, 2));
#     double distance = sqrt( (ds/s) + (dc/compactness));
#     return distance;
#   }
# ')




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
  nfound=0
  batchsize=10
  sample_size=1000
  cds_cand <- cds
  sel <- c()
  iter <- 1
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
    iter <- iter+1
  }

  list(selected=sel, coords=cds[sel[1:K],])

}

#
# get_26_connected_neighbors <- function(i, j, k, array_3d) {
#   neighbors <- list()
#   max_i <- dim(array_3d)[1]
#   max_j <- dim(array_3d)[2]
#   max_k <- dim(array_3d)[3]
#
#   for (i_shift in -1:1) {
#     for (j_shift in -1:1) {
#       for (k_shift in -1:1) {
#         # Skip the center point (0, 0, 0) shift
#         if (i_shift == 0 && j_shift == 0 && k_shift == 0) {
#           next
#         }
#
#         i_neighbor <- i + i_shift
#         j_neighbor <- j + j_shift
#         k_neighbor <- k + k_shift
#
#         # Check if the neighbor is within the array bounds
#         if (i_neighbor >= 1 && i_neighbor <= max_i &&
#             j_neighbor >= 1 && j_neighbor <= max_j &&
#             k_neighbor >= 1 && k_neighbor <= max_k) {
#           neighbors <- append(neighbors, list(c(i_neighbor, j_neighbor, k_neighbor)))
#         }
#       }
#     }
#   }
#
#
#
#   neighbors[sidx,,drop=FALSE]
# }



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
#' @import assertthat
#'
#' @export
snic <- function(vec, mask, compactness=5, K=500) {
  mask.idx <- which(mask>0)
  valid_coords <- index_to_grid(mask, mask.idx)
  norm_coords <- sweep(valid_coords, 2, spacing(mask), "/")

  mask.grid <- index_to_grid(mask, mask.idx)
  mask_lookup <- array(0, dim(mask))
  mask_lookup[mask.grid] <- 0:(length(mask.idx)-1)

  vecmat <- series(vec, mask.idx)
  vecmat <- scale(vecmat)
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


  kvol <- ClusteredNeuroVol(as.logical(mask), clusters=ret[mask.idx])

  browser()

  structure(
    list(clusvol=kvol,
         gradvol=grad,
         cluster=ret[mask.idx],
         centers=ret$center,
         coord_centers=ret$coord_centers),
    class=c("snic_cluster_result", "cluster_result", "list"))

}


#
# # Define the SNIC function
# snic <- function(vec, mask, K=300, compactness=20) {
#   # Preprocess the image, compute the normalizing factors, and obtain the initial centroids
#   # ...
#
#   pushcount = 1
#
#   mask.idx <- which(mask>0)
#   valid_coords <- index_to_grid(mask, mask.idx)
#   init <- find_initial_points(valid_coords, K)
#
#   #centroids <- initialize_centroids(valid_coords,K )
#   #pres <- kmeans(scale(valid_coords), K, iter.max =500)
#   #out <- RANN::nn2(scale(valid_coords), scale(pres$centers), 2)
#   #centroid_idx <- out$nn.idx[,1]
#   centroids <- valid_coords[init$selected,]
#
#   mask.grid <- index_to_grid(mask, mask.idx)
#   mask_lookup <- array(0, dim(mask))
#   mask_lookup[mask.grid] <- 1:length(mask.idx)
#
#   vecmat <- series(vec, mask.idx)
#   vecmat <- apply(vecmat, 2, function(x) x/norm(x, "2"))
#
#   #centroids <- pres$medoids
#   # Initialize the label map and priority queue
#
#   pushcount=0
#   L <- array(0, dim(mask))
#
#   Q <- priority_queue()
#   C <- vector("list", K)
#
#   s <- sqrt(nrow(valid_coords)/K)
#   # Populate the priority queue with initial elements
#   for (k in 1:K) {
#     x_k <- centroids[k,,drop=FALSE]
#     c_k <- vecmat[,centroid_idx[k]]
#     #e <- list(x = x_k, c = c_k, label = k, distance = 0)
#     e <- list(x = x_k, i_k = centroid_idx[k], label = k, distance = 0)
#
#     C[[k]] <- list(x = x_k, c = c_k, n = 1)
#     Q$push(e, priority = .Machine$integer.max-1)
#     pushcount = pushcount+1
#   }
#
#   nassigned <- 0
#   nvoxels <- length(mask.idx)
#   # Main loop
#   while (Q$size() > 0 && nassigned < nvoxels) {
#     #print(paste("Q:", Q$size()))
#     print(paste("L assigned", nassigned))
#     # Pop the element with the smallest distance
#     e_i <- Q$pop()
#     x_i <- as.vector(e_i$x)
#     #c_i <- as.vector(e_i$c)
#     #c_i <- scale(series(vec, e_i$x))
#     c_i <- vecmat[,e_i$i_k]
#     k_i <- e_i$label
#
#     # Check if the pixel has not been labeled yet
#     if (L[rbind(x_i)] == 0) {
#       L[rbind(x_i)] <- k_i
#       nassigned <- nassigned + 1
#       # Update the centroid with the current pixel
#       C[[k_i]] <- update_centroid_online(C[[k_i]], x_i, c_i)
#
#       cen_k <- C[[k_i]]$x
#       # Iterate through the 4-connected or 8-connected neighbors of the current pixel
#       neighbors <- get_26_connected_neighbors_rcpp(x_i[1], x_i[2], x_i[3], dim(mask))
#       neighbors <- do.call(rbind, neighbors)
#
#       #neighbor_idx <- grid_to_index(mask, neighbors)
#
#       neighbors <- neighbors[mask[neighbors]==TRUE & L[neighbors] == 0,,drop=FALSE]
#       print(paste("neighbors:", nrow(neighbors)))
#
#       if (nrow(neighbors) == 0) {
#         next
#       }
#
#       neighb_idx <- mask_lookup[neighbors]
#       compactness=10
#       # tmp <- lapply(1:nrow(neighbors), function(i) {
#       #   x_j <- neighbors[i,,drop=FALSE]
#       #   c_j <- vecmat[,neighb_idx[i]]
#       #   ds <- sum((x_i - x_j)^2)
#       #   dc <- 1- cor(C[[k_i]]$c, c_j)
#       #   d_j_ki <- round(sqrt((ds/s + dc/compactness)) * 1000000)
#       #   e_j <- list(x = x_j, i_k = neighb_idx[i], label = k_i, distance = d_j_ki)
#       #   list(x_j=x_j, c_j=c_j,ds=ds, dc=dc,d_j_ki, e_j=e_j, priority = round(.Machine$integer.max- d_j_ki))
#       # })
#
#       for (i in 1:nrow(neighbors)) {
#         x_j <- neighbors[i,,drop=FALSE]
#         # Check if the neighbor pixel has not been labeled yet
#         if (L[x_j] == 0) {
#           c_j <- vecmat[,neighb_idx[i]]
#           #c_j <- series(vec, x_j)
#           #c_j <- c_j / norm(c_j, "2")
#
#           # Compute the distance between the neighbor pixel and the current centroid
#           #ds <- sum((cen_k - x_j)^2)
#           #dc <- 1- cor(C[[k_i]]$c, c_j)
#
#           #print(paste("ds = ", ds))
#           #print(paste("dc = ", dc))
#
#           #d_j_ki <- round(sqrt((ds/s + dc/compactness)) * 1000000)
#           #d_j_ki <- ds
#
#
#           d_j_ki <- compute_distance_rcpp(cen_k/norm(x_i, "2"),
#                                           x_j/norm(x_k, "2"),
#                                           C[[k_i]]$c,
#                                           c_j, s=s, compactness) * 10000
#           print(paste("djki = ", d_j_ki))
#
#           # Create a new element for the neighbor pixel
#           #e_j <- list(x = x_j, c = c_j, label = k_i, distance = d_j_ki)
#           #e <- list(x = x_k, i_k = centroid_idx[k], label = k, distance = 0)
#           e_j <- list(x = x_j, i_k = neighb_idx[i], label = k_i, distance = round(d_j_ki))
#           #print(paste("pushing", x_j, "at distance ", d_j_ki))
#           # Push the new element into the priority queue
#           pushcount = pushcount+1
#           if (pushcount %% 1000 == 0) {
#             image(L[,,15], col=rainbow(256))
#             pushcount=0
#           }
#           #print(paste("pushcount", pushcount))
#           Q$push(e_j, priority = round(.Machine$integer.max- d_j_ki))
#         }
#       }
#     }
#   }
#
#   return(L)
# }
#
# # Run the SNIC algorithm
# #label_map <- snic(input_image, K = 200, compactness = 20)
#
# # Define the update_centroid_online function
# update_centroid_online <- function(centroid, x_i, c_i) {
#   # centroid: A list containing the current centroid information with elements x (spatial position), c (color), and n (number of pixels)
#   # x_i: Spatial position of the current pixel
#   # c_i: Color value of the current pixel
#
#   # Calculate the new spatial position
#   new_x <- (centroid$x * centroid$n + x_i) / (centroid$n + 1)
#
#   # Calculate the new color value
#   new_c <- (centroid$c * centroid$n + c_i) / (centroid$n + 1)
#   #new_c <- new_c / norm(new_c, "2")
#
#   # Update the number of pixels in the superpixel
#   new_n <- centroid$n + 1
#
#   # Return the updated centroid as a list
#   return(list(x = new_x, c = new_c/norm(new_c, "2"), n = new_n))
# }

