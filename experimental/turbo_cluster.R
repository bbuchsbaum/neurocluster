


#' tesselate a mask volume into K clusters with kmeans applied to spatial coordinates
#'
#' @param mask an instance of type \code{NeuroVol}
#' @param K the number of clusters
#' @export
tesselate <- function(mask, K=100) {
  mask.idx <- which(mask>0)
  grid <- index_to_coord(mask, mask.idx)

  gcen <- grid[as.integer(seq(1, nrow(grid), length.out=K)),]
  kres <- kmeans(grid, centers=gcen, iter.max=500)
  clusvol <- ClusteredNeuroVol(mask, kres$cluster)
}



# supervoxelR <- function(mask, bvec, K=500, lambda=.5, iterations=25, connectivity=27, use_medoid=FALSE) {
#   assert_that(lambda >= 0 && lambda <= 1)
#   assert_that(connectivity > 1 & connectivity <= 27)
#
#   mask.idx <- which(mask > 0)
#   grid <- indexToCoord(mask, mask.idx)
#   vgrid <- indexToGrid(mask, mask.idx)
#
#   #if (shrink > 0) {
#   #  message("running ", shrink, "shrinkage iterations with k = ", 5)
#   #  bvec <- shrink_vector(mask, vgrid, bvec, k=5, iter=shrink)
#   #}
#
#   kres <- kmeans(grid, K, iter.max=500)
#   kvol <- NeuroVol(kres$cluster, space(mask), indices=mask.idx)
#
#   clusid <- sort(unique(kres$cluster))
#   neib <- get.knn(grid, k=connectivity)
#
#   iter <- 1
#   switches <- 1
#   iter.max <- iterations
#
#   feature_mat <- series(bvec, mask.idx)
#
#   centroids <- compute_centroids(feature_mat, grid, kres$cluster, medoid=use_medoid)
#   sp_centroids <- do.call(rbind, lapply(centroids, "[[", "centroid"))
#   num_centroids <- do.call(rbind, lapply(centroids, "[[", "center"))
#
#   denom <- max(get.knn(sp_centroids, k=1)$nn.dist[,1])
#
#
#   curclus <- kvol[mask.idx]
#
#
#   ## port iteration to rcpp
#   while (iter < iter.max && switches > 0) {
#     message("turboclust, iteration: ", iter)
#
#     newclus <- sapply(1:nrow(vgrid), function(i) {
#
#       ind <- neib$nn.index[i,]
#       D <- neib$nn.dist[i,]
#       keep <- which(D < dthresh)
#       ind <- ind[keep]
#
#       oclus <- unique(kvol[mask.idx[ind]])
#       diffclus <- which(oclus != curclus[i])
#
#       if (length(diffclus) > 0) {
#
#         candclus <- c(curclus[i], oclus[diffclus])
#
#         cval <- sapply(candclus, function(cc) {
#           (1 - cor(feature_mat[,i], centroids[[cc]]$center))/2
#         })
#
#         dvals <- sapply(candclus, function(cc) {
#           sqrt(sum((grid[i,] - centroids[[cc]]$centroid)^2))
#         })/denom
#
#         cost <- lambda*cval + (1-lambda)*dvals
#         #cost <- cvals
#         newc <- candclus[which.min(cost)]
#       } else {
#         curclus[i]
#       }
#     })
#
#
#
#     centroids <- compute_centroids(bvec, mask.idx, grid, newclus, medoid=use_medoid)
#     switches <- sum(newclus != curclus)
#
#     curclus <- newclus
#     message("tuboclust: nswitches ", switches)
#
#     iter <- iter + 1
#   }
#
#   kvol <- NeuroVol(newclus, space(mask), indices=mask.idx)
#   voxsplit <- split(data.frame(vgrid), newclus)
#
#   list(clusvol=kvol,
#        clusters=newclus,
#        centers = do.call(rbind, lapply(centroids, "[[", "center")),
#        coordCenters=do.call(rbind, lapply(centroids, "[[", "centroid")),
#        voxelSets=voxsplit)
# }




#' @keywords internal
init_cluster <- function(bvec, mask, grid, K) {
  mask.idx <- which(mask>0)
  gcen <- grid[as.integer(seq(1, nrow(grid), length.out=K)),]
  kres <- kmeans(grid, centers=gcen, iter.max=500)
  clusvol <- NeuroVol(kres$cluster, space(mask), indices=mask.idx)
  kres$cluster
}


#' @inheritParams supervoxels
#' @param initclus
#' @importFrom assertthat assert_that
#' @import neuroim2
supervoxel_cluster_fit <- function(feature_mat, grid, K=min(500, nrow(grid)),sigma1=1, sigma2=5,
                              alpha=.5, iterations=25, connectivity=26, use_medoid=FALSE,
                              initclus, use_gradient=TRUE) {

  message("supervoxel_fit: sigma1 = ", sigma1, " sigma2 = ", sigma2)

  assert_that(connectivity > 1 & connectivity <= 27)
  assert_that(alpha >= 0 && alpha <= 1)

  feature_mat <- scale(feature_mat, center=TRUE, scale=TRUE)

  if (is.null(initclus)) {
    ## kmeans using coordinates only
    gcen <- grid[as.integer(seq(1, nrow(grid), length.out=K)),]
    kres <- kmeans(grid, centers=gcen, iter.max=500)

    clusid <- sort(unique(kres$cluster))
    curclus <- kres$cluster

    seeds <- FNN::get.knnx(grid, kres$centers)$nn.index[,1]
    sp_centroids <- grid[seeds,]
    num_centroids <- feature_mat[, seeds]

  } else {
    assert_that(length(initclus) == nrow(grid))
    clusid <- sort(unique(initclus))
    assert_that(length(clusid) == K)
    curclus <- initclus
    centroids <- compute_centroids(feature_mat, grid, curclus, medoid=use_medoid)
    sp_centroids <- do.call(rbind, centroids$centroid)
    num_centroids <- do.call(rbind, centroids$center)

    #browser()
  }

  ## find k neighbors within 'connectivity' radius
  neib <- FNN::get.knn(grid, k=connectivity)

  ## has to be changed for surface... pass in neighbor_fun?
  dthresh <- median(neib$nn.dist[,connectivity,drop=FALSE])
  message("dthresh: ", dthresh)

  iter <- 1
  switches <- 1
  iter.max <- iterations

  while (iter < iter.max && switches > 0) {

    candlist <- find_candidates(neib$nn.index-1, neib$nn.dist, curclus, dthresh)
    #browser()
    newclus <- best_candidate(candlist, curclus, t(grid),
                              t(num_centroids), t(sp_centroids),
                              feature_mat, sigma1, sigma2, alpha)

    switches <- attr(newclus, "nswitches")

    if (switches > 0) {
      centroids <- compute_centroids(feature_mat, grid, newclus, medoid=use_medoid)
      sp_centroids <- do.call(rbind, centroids$centroid)
      num_centroids <- do.call(rbind, centroids$center)
      curclus <- newclus
    }

    message("supervoxels_fit: iter ", iter, " -- num  switches =  ", switches)
    iter <- iter + 1
  }

  list(clusters=newclus,
       centers = do.call(rbind, lapply(centroids, "[[", "center")),
       coord_centers=do.call(rbind, lapply(centroids, "[[", "centroid")))
}

#' knn_shrink
#'
#' Replace each voxel by the mean of its k nearest neighbors in its local spatial neighborhood
#'
#' @param bvec bvec a \code{\linkS4class{NeuroVec}} instance
#' @param mask a mask defining the voxels to include in the clustering result.
#' @param k the number of nearest neighbors to average over.
#' @param connectivity the number of spatial neighbors to include in the search space around each voxel.
#' @export
#'
#' @examples
#' mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
#' bvec <- replicate(10, NeuroVol(array(runif(20*20*20), c(20,20,20)),
#' NeuroSpace(c(20,20,20))), simplify=FALSE)
#' bvec <- do.call(concat, bvec)
#'
#' sbvec <- knn_shrink(bvec, mask,k=3)
#'
knn_shrink <- function(bvec, mask, k=5, connectivity=27) {
  assert_that(inherits(bvec, "NeuroVec"))
  assert_that(inherits(mask, "NeuroVol"))
  assert_that(k >= 1)
  assert_that(connectivity >= k)

  mask <- as.logical(mask)
  mask.idx <- which(mask)
  grid <- index_to_coord(mask, mask.idx)

  feature_mat <- series(bvec, mask.idx)

  ## find k neighbors within 'connectivity' radius
  neib <- FNN::get.knn(grid, k=connectivity)

  sfeature_mat <- t(do.call(rbind, map(1:nrow(neib$nn.index), function(i) {
    rowMeans(feature_mat[, c(i, neib$nn.index[i,1:(k-1)])])
  })))

  SparseNeuroVec(sfeature_mat, space(bvec), mask=mask)
}

#' supervoxels
#'
#' Cluster a \code{NeuroVec} instance into a set of spatially constrained clusters.
#'
#'
#' @param bvec a \code{\linkS4class{NeuroVec}} instance supplying the data to cluster.
#' @param mask a \code{\linkS4class{NeuroVol}} mask defining the voxels to include in the clustering result.
#'    If the mask contains \code{numeric} data, nonzero values will define the included voxels. If the mask
#'    is a \code{\linkS4class{LogicalNeuroVol}} then \code{TRUE} will define the set of included voxels.
#' @param K the number of clusters to find.
#' @param sigma1 the bandwidth of the heat kernel for computing similarity of the data vectors.
#' @param sigma2 the bandwidth of the heat kernel for computing similarity of the coordinate vectors.
#'    If this value is small, then relatively larger weights are given to nearby voxels. If is is large, then
#'    spatial weights will be less salient. A relatively large sigma1/sigma2 ratio weights data features more than
#'    spatial features, whereas as large sigma2/sigma1 ration does the opposite.
#' @param iterations the maximum number of cluster iterations
#' @param connectivity
#' @param use_medoid whether to use the medoids to define cluster centroids
#' @param alpha the relative weighting between spatial coherence and feature similarity metrics.
#' alpha = 0, means all feature-weighting, alpha=0 is spatial weighting. Default is .5.
#'
#' @param filter low- and high-pass filter parameters. See details.
#' @export
#' @importFrom neuroim2 NeuroVec NeuroVol
#' @import furrr
#' @import assertthat
#' @return a \code{list} of class \code{supervoxels_cluster_result} with the following elements:
#'
#' \describe{
#'   \item{clustervol}{an instance of type \linkS4class{ClusteredNeuroVol}}
#'   \item{clusters}{ a vector of cluster indices equal to the number of voxels in the mask}
#'   \item{centers}{ a matrix of cluster centers with each column representing the feature vector for a cluster }
#'   \item{coord_centers}{ a matrix of spatial coordinates with each row corresponding to a cluster }
#' }
#'
#' @examples
#'
#' mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
#' bvec <- replicate(10, NeuroVol(array(runif(20*20*20), c(20,20,20)),
#' NeuroSpace(c(20,20,20))), simplify=FALSE)
#' bvec <- do.call(concat, bvec)
#'
#' cres1 <- supervoxels(bvec, mask, K=100, sigma1=1, sigma2=10, sample_frac=.3)
#' cres2 <- supervoxels(bvec, mask, K=100, sigma1=1, sigma2=6, sample_frac=.3)
#' cres3 <- supervoxels(bvec, mask, K=100, sigma1=1, sigma2=4, sample_frac=.3)
#' cres4 <- supervoxels(bvec, mask, K=100, sigma1=1, sigma2=2, sample_frac=.3)
#' cres5 <- supervoxels(bvec, mask, K=100, sigma1=.5, sigma2=1, sample_frac=.3)
#'
#' cres_cons <- merge_clus(cres1, cres2, cres3, cres4, cres5)
#'
#' cres_cons2 <- supervoxels(bvec, mask, K=100, sigma1=c(1,2,3), sigma2=c(3,2,1), sample_frac=.3)
#'
#' ## to access the cluster volume: cres1$clusvol
#'
#' @details
#'
#' Here's a brief overview of the algorithm:
#' 1. It first scales input data (bvec)
#'
#' 2. It initializes the clusters using the gradient of the image and a furthest neighbor approach.
#'
#' 3. It then runs the "supervoxel_fit" function, which iteratively finds spatially constrained clusters using a combination of feature similarity and spatial similarity. The feature similarity is calculated based on the heat kernel with bandwidth sigma1, while the spatial similarity is calculated based on the heat kernel with bandwidth sigma2. The relative weighting between these two similarities is controlled by the alpha parameter.
#'
#' 4. Finally, it returns a list containing the clustered NeuroVol object (clusvol), cluster assignments for each voxel (cluster), the feature vector for each cluster (centers), and the spatial coordinates of each cluster (coord_centers).
#' This algorithm is particularly suitable for brain data analysis, as it takes both the spatial and feature similarities into account, allowing it to identify clusters that are spatially coherent and share similar features.
#'
supervoxel <- function(bvec, mask, K=500, sigma1=1,
                                sigma2=2.5, iterations=50,
                                connectivity=27, use_medoid=FALSE,
                                alpha=.5,
                                filter=c(lp=0, hp=0),
                                filter_method=c("bspline", "locfit"),
                                sample_frac=1, nreps=1) {


  if (length(sigma1) != nreps) {
    sigma1 <- rep(sigma1, length.out=nreps)
  }

  if (length(sigma2) != nreps) {
    sigma2 <- rep(sigma2, length.out=nreps)
  }

  filter_method <- match.arg(filter_method)
  mask.idx <- which(mask > 0)

  ## coordinate grid in mm units
  grid <- index_to_coord(mask, mask.idx)

  if (any(filter > 0)) {
    assert_that(filter$lp <= 1 && filter$hp <= 1)
    message("supervoxel_cluster: pre-filtering time series data with ", filter_method, " smoothing")
    bvec <- filter_vec(bvec, mask, filter$lp, filter$hp, method=filter_method)
  }

  clusinit <- init_cluster(bvec, mask, grid, K)

  feature_mat <- neuroim2::series(bvec, mask.idx)

  res <- purrr::map(1:nreps, function(i) {
    feature_mat <- if (sample_frac < 1) {
      n <- max(sample_frac * nrow(feature_mat),1)
      assert_that(n > 1)
      keep <- sample(1:nrow(feature_mat), n)
      feature_mat[keep,,drop=FALSE]
    } else {
      feature_mat
    }

    ret <- supervoxel_cluster_fit(feature_mat, grid, sigma1=sigma1[i], sigma2=sigma2[i], K=K,
                                iterations=iterations, connectivity=connectivity,
                                use_medoid=use_medoid, alpha=alpha,initclus=clusinit)

    kvol <- ClusteredNeuroVol(as.logical(mask), clusters=ret$clusters)

    structure(
         list(clusvol=kvol,
              cluster=ret$clusters,
              centers=ret$center,
              coord_centers=ret$coord_centers),
         class=c("cluster_result", "list"))
  })

  if (length(res) > 1) {
    message("computing consensus clustering over ", length(res), " partitions")
    kvol <- ClusteredNeuroVol(as.logical(mask), do.call(merge_clus, res))

    cen <- compute_centroids(feature_mat, grid,assignment=kvol@clusters, medoid=use_medoid)
    structure(
      list(clusvol=kvol,
           cluster=kvol@clusters,
           centers=cen$center,
           coord_centers=cen$coord_centers),
      class=c("cluster_result", "list"))
  } else {
    res[[1]]
  }

}

#' @rdname merge_clus
#'
#' @importFrom clue cl_consensus as.cl_hard_partition cl_ensemble
#' @return a \code{linkS4class{ClusteredNeuroVol}} instance
#' @export
merge_clus.cluster_result <- function(x, ...) {
  args <- c(list(x), list(...))
  assert_that(all(map_lgl(args, ~ inherits(., "cluster_result"))))

  ens <- do.call(clue::cl_ensemble, args)
  cons <- clue::cl_consensus(ens, method="SE")
  hpart <- clue::as.cl_hard_partition(cons)
  as.integer(hpart$.Data)
  #ClusteredNeuroVol(x$clusvol@mask, clusters=hpart$.Data)
}


#' @export
merge_clus.cluster_result_time <- function(x, ...) {
  args <- c(list(x), list(...))
  assert_that(all(map_lgl(args, ~ inherits(., "cluster_result"))))

  ens <- do.call(clue::cl_ensemble, args)
  cons <- clue::cl_consensus(ens)
  hpart <- clue::as.cl_hard_partition(cons)
  as.integer(hpart$.Data)
}

#' @export
cl_class_ids.cluster_result <- function(x) {
  x$cluster
}


#' @export
is.cl_partition.cluster_result <- function(x) {
  TRUE
}


#' filter_vec
#'
#' bandpass filter a \code{\linkS4class{NeuroVec}} over the 4th (temporal) dimesion
#'
#' @inheritParams filter_mat
#' @param bvec an instance of type \code{\linkS4class{NeuroVec}} to filter
#' @param mask a binary mask defining the voxels to include in the filter operation
filter_vec <- function(bvec, mask, lp=0, hp=.7, method=c("locfit", "bspline")) {
  feature_mat <- series(bvec, which(mask!=0))
  fmat <- filter_mat(feature_mat, lp,hp,method)
  SparseNeuroVec(fmat, space(bvec), mask=as.logical(mask))
}

#' filter_mat
#'
#' temporally filter a matrix of time-series.
#'
#'
#' @param feature_mat a feature matrix with \code{nrow(feature_mat)} observations and `\code{ncol(feature_mat)} features.
#' @param lp the low-pass bandwidth parameter specified as a fraction
#'    of the number of time-points to include in each window. Must be between [0-1]
#'  @param hp the high-pass bandwidth parameter specified as a fraction
#'    of the number of time-points to include in each window. Must be between [0-1]
#' @importFrom locfit locfit lp
#' @examples
#'
#' feature_mat <- matrix(rnorm(200*10), 200, 10)
#' library(future)
#' plan(multicore)
#' sfeature_mat <- filter_mat(feature_mat, lp=.2, hp=.7)
#' sfeature_mat2 <- filter_mat(feature_mat, lp=.2, hp=.7, method="bspline")
#'
filter_mat <- function(feature_mat, lp=0, hp=.7, method=c("locfit", "bspline")) {
  assert_that(hp >= lp)
  assert_that(hp >= 0 && hp <= 1)
  assert_that(lp >= 0 && lp <= 1)


  method <- match.arg(method)

  if (method == "locfit") {
    nn <- hp
    nn1 <- lp
    time <- seq(1,nrow(feature_mat))

    if (lp ==0 && hp == 0) {
      warning("filter_mat: returning oriignal matrix unchanged, both 'lp' and 'hp' are 0")
      return(feature_mat)
    }

    feature_mat <- do.call(cbind, furrr::future_map(1:ncol(feature_mat), function(i) {
      vals <- feature_mat[,i]
      fvals <- vals - fitted(locfit(vals ~ lp(time, nn=hp)))

      if (lp > 0) {
        fvals <- fitted(locfit(fvals ~ lp(time, nn=lp)))
      }
      fvals
    }))
  } else {
    ## bspline method
    feature_mat <- if (hp > 0) {
      nhp <- ceiling(nrow(feature_mat)/(hp * nrow(feature_mat)))
      message("filter_mat: number of highpass basis functions: ", nhp)
      bmat1 <- if (nhp <= 2) {
        splines::ns(1:nrow(feature_mat), nhp)
      } else {
        splines::bs(1:nrow(feature_mat), nhp)
      }

      resid(lsfit(bmat1, feature_mat))
    } else {
      feature_mat
    }

    feature_mat <- if (lp > 0) {
      nlp <- ceiling(nrow(feature_mat)/(lp * (nrow(feature_mat)*.7)))
      if (nlp > nrow(feature_mat)) {
        nlp <- .8 * nrow(feature_mat)
      }
      message("filter_mat: number of low-pass basis functions: ", nlp)
      bmat1 <- splines::bs(1:nrow(feature_mat), nlp)
      fitted(lm(feature_mat ~ bmat1))
    } else {
      feature_mat
    }

  }
}

#' cluster objects with a temporal constraint
#'
#' @export
#' @inheritParam supervoxel_cluster
#' @examples
#' feature_mat <- matrix(rnorm(100*10), 100, 10)
#' library(future)
#' plan(multicore)
#' cres <- supervoxel_cluster_time(t(feature_mat), K=20)
#'
supervoxel_cluster_time <- function(feature_mat, K=min(nrow(feature_mat), 100), sigma1=1, sigma2=TR*3,
                               iterations=50, TR=2, filter=list(lp=0, hp=0),
                               use_medoid=FALSE, nreps=5) {
  if (filter$lp > 0 || filter$hp > 0) {
    message("supervoxel_cluster: filtering time series")
    feature_mat <- filter_mat(feature_mat, filter$lp, filter$hp)
  }

  nels <- nrow(feature_mat)
  grid <- as.matrix(seq(0, by=TR, length.out=nels))

  fits <- furrr::future_map(1:nreps, function(i) {
    initsamps <- sort(sample(1:nrow(feature_mat), K))
    initcenters <- feature_mat[initsamps,]

    curclus <- FNN::get.knnx(grid[initsamps,,drop=FALSE], grid, k=1)$nn.index[,1]

    ret <- supervoxel_cluster_fit(t(feature_mat), as.matrix(grid),
                           sigma1=sigma1, sigma2=sigma2,
                           K=K,
                           initclus=curclus,
                           iterations=iterations,
                           connectivity=3,
                           use_medoid=use_medoid)

    class(ret) <- c("cluster_result_time", "cluster_result", "list")
    ret
  })

}


#' @export
#' @import neurosurf
supervoxel_cluster_surface <- function(bsurf, K=500, sigma1=1, sigma2=5, iterations=50,
                                  connectivity=6, use_medoid=FALSE, filter=list(lp=-1, hp=-1)) {
  mask.idx <- indices(bsurf)
  grid <- coords(bsurf)[mask.idx,]

  feature_mat <- series(bsurf, mask.idx)

  if (any(filter) > 0) {
    message("supervoxel_cluster: filtering time series")
    feature_mat <- filter_mat(feature_mat, filter$lp, filter$hp)
  }

  ret <- supervoxel_cluster_fit(feature_mat, grid, sigma1=sigma1, sigma2=sigma2, K=K,
                           iterations=iterations, connectivity=connectivity,
                           use_medoid=use_medoid)


  kvol <- BrainSurface(geometry=geometry(bsurf), indices=mask.idx, data=ret$clusters)
  index_sets <- split(mask.idx, ret$clusters)

  ret <- list(clusvol=kvol,
              clusters=ret$clusters,
              centers = ret$centers,
              coord_centers=ret$coord_centers,
              index_sets=index_sets)

  class(ret) <- c("cluster_result_surface", "cluster_result", "list")
}




