

computeCentroids <- function(valmat, grid, clusters, assignment, medoid=FALSE) {
  csplit <- split(1:length(assignment), assignment)
  if (!medoid) {

    lapply(csplit, function(id) {
      mat <- valmat[, id]
      coords <- grid[id,,drop=FALSE]
      list(center=rowMeans(mat), centroid=colMeans(coords))
    })
  } else {
    lapply(csplit, function(id) {

      mat <- valmat[, id]
      coords <- grid[idx,,drop=FALSE]
      coords_dist <- as.matrix(dist(coords))
      coords_medoid_ind <- which.min(rowSums(coords_dist))
      Dmat <- 1-cor(mat)
      mat_medoid_ind <- which.min(rowSums(Dmat))
      list(center=mat[,mat_medoid_ind], centroid=coords[coords_medoid_ind,])
    })

  }
}

## try multiple kmeans initializations, choose one with best intra-cluster correlation.

#' turbo_cluster
#'
#' @export
#' @import FNN
#' @import assertthat
#'
turbo_clusterR <- function(mask, bvec, K=500, lambda=.5, iterations=25, connectivity=27, shrink=0, use_medoid=FALSE) {
  assert_that(lambda >= 0 && lambda <= 1)
  assert_that(connectivity > 1 & connectivity <= 27)

  mask.idx <- which(mask > 0)
  grid <- indexToCoord(mask, mask.idx)
  vgrid <- indexToGrid(mask, mask.idx)

  if (shrink > 0) {
    message("running ", shrink, "shrinkage iterations with k = ", 5)
    bvec <- shrink_vector(mask, vgrid, bvec, k=5, iter=shrink)
  }

  kres <- kmeans(grid, K, iter.max=500)
  kvol <- BrainVolume(kres$cluster, space(mask), indices=mask.idx)

  clusid <- sort(unique(kres$cluster))
  neib <- get.knn(grid, k=connectivity)
  dthresh <- min(neib$nn.dist[,connectivity])

  iter <- 1
  switches <- 1
  iter.max <- iterations

  centroids <- computeCentroids(bvec, mask.idx, grid, sort(unique(kres$cluster)), kres$cluster, medoid=use_medoid)
  sp_centroids <- do.call(rbind, lapply(centroids, "[[", "centroid"))
  num_centroids <- do.call(rbind, lapply(centroids, "[[", "center"))

  denom <- max(get.knn(sp_centroids, k=1)$nn.dist[,1])
  valmat <- series(bvec, mask.idx)

  curclus <- kvol[mask.idx]


  ## port iteration to rcpp
  while (iter < iter.max && switches > 0) {
    message("turboclust, iteration: ", iter)

    newclus <- sapply(1:nrow(vgrid), function(i) {

      ind <- neib$nn.index[i,]
      D <- neib$nn.dist[i,]
      keep <- which(D < dthresh)
      ind <- ind[keep]

      oclus <- unique(kvol[mask.idx[ind]])
      diffclus <- which(oclus != curclus[i])

      if (length(diffclus) > 0) {

        candclus <- c(curclus[i], oclus[diffclus])

        cval <- sapply(candclus, function(cc) {
          (1 - cor(valmat[,i], centroids[[cc]]$center))/2
        })

        dvals <- sapply(candclus, function(cc) {
          sqrt(sum((grid[i,] - centroids[[cc]]$centroid)^2))
        })/denom

        cost <- lambda*cval + (1-lambda)*dvals
        #cost <- cvals
        newc <- candclus[which.min(cost)]
      } else {
        curclus[i]
      }
    })



    centroids <- computeCentroids(bvec, mask.idx, grid, sort(unique(newclus)), newclus, medoid=use_medoid)
    switches <- sum(newclus != curclus)

    curclus <- newclus
    message("tuboclust: nswitches ", switches)

    iter <- iter + 1
  }

  kvol <- BrainVolume(newclus, space(mask), indices=mask.idx)
  voxsplit <- split(data.frame(vgrid), newclus)

  list(clusvol=kvol,
       clusters=newclus,
       centers = do.call(rbind, lapply(centroids, "[[", "center")),
       coordCenters=do.call(rbind, lapply(centroids, "[[", "centroid")),
       voxelSets=voxsplit)
}

turbo_cluster <- function(mask, bvec, K=500, sigma1=1, sigma2=10, iterations=25, connectivity=27, shrink=0, use_medoid=FALSE) {
  assert_that(sigma1 >= 0 && sigma2 <= 1)
  assert_that(connectivity > 1 & connectivity <= 27)

  mask.idx <- which(mask > 0)
  grid <- indexToCoord(mask, mask.idx)
  vgrid <- indexToGrid(mask, mask.idx)

  if (shrink > 0) {
    message("running ", shrink, "shrinkage iterations with k = ", 5)
    bvec <- shrink_vector(mask, vgrid, bvec, k=5, iter=shrink)
  }

  kres <- kmeans(grid, K, iter.max=500)
  kvol <- BrainVolume(kres$cluster, space(mask), indices=mask.idx)

  clusid <- sort(unique(kres$cluster))
  neib <- FNN::get.knn(grid, k=connectivity)
  dthresh <- min(neib$nn.dist[,connectivity])

  iter <- 1
  switches <- 1
  iter.max <- iterations

  valmat <- series(bvec, mask.idx)
  centroids <- computeCentroids(valmat, grid, clusid, kres$cluster, medoid=use_medoid)
  sp_centroids <- do.call(rbind, lapply(centroids, "[[", "centroid"))
  num_centroids <- do.call(rbind, lapply(centroids, "[[", "center"))

  curclus <- as.integer(kvol[mask.idx])


  while (iter< iter.max && switches > 0) {
    candlist <- find_candidates(neib$nn.index-1, neib$nn.dist, curclus, dthresh)

    newclus <- best_candidate(candlist, curclus, t(grid), t(num_centroids), t(sp_centroids), valmat, sigma1, sigma2)

    switches <- attr(newclus, "nswitches")

    centroids <- computeCentroids(valmat, grid, sort(unique(newclus)), newclus, medoid=use_medoid)

    sp_centroids <- do.call(rbind, lapply(centroids, "[[", "centroid"))
    num_centroids <- do.call(rbind, lapply(centroids, "[[", "center"))


    message("switches: ", switches)
    curclus <- newclus
    iter <- iter + 1
  }

  kvol <- BrainVolume(newclus, space(mask), indices=mask.idx)
  voxsplit <- split(data.frame(vgrid), newclus)

  list(clusvol=kvol,
       clusters=newclus,
       centers = do.call(rbind, lapply(centroids, "[[", "center")),
       coordCenters=do.call(rbind, lapply(centroids, "[[", "centroid")),
       voxelSets=voxsplit)
}
