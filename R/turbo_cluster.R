


# get centroids for a matrix and set of assignments
#' @keywords internal
#' @importFrom purrr map
compute_centroids <- function(valmat, grid, assignment, medoid=FALSE) {

  csplit <- split(1:length(assignment), assignment)

  if (!medoid) {
    map(csplit, function(id) {
      mat <- valmat[, id]
      coords <- grid[id,,drop=FALSE]
      list(center=rowMeans(mat), centroid=colMeans(coords))
    })
  } else {
    map(csplit, function(id) {
      mat <- valmat[, id, drop=FALSE]
      coords <- grid[id,,drop=FALSE]
      coords_dist <- as.matrix(dist(coords))
      coords_medoid_ind <- which.min(rowSums(coords_dist))
      Dmat <- 1-cor(mat)
      mat_medoid_ind <- which.min(rowSums(Dmat))
      list(center=mat[,mat_medoid_ind], centroid=coords[coords_medoid_ind,])
    })

  }
}

## try multiple kmeans initializations, choose one with best intra-cluster correlation.

#' @import FNN
#' @import assertthat
#' @inheritParams turbo_cluster
turbo_clusterR <- function(mask, bvec, K=500, lambda=.5, iterations=25, connectivity=27, use_medoid=FALSE) {
  assert_that(lambda >= 0 && lambda <= 1)
  assert_that(connectivity > 1 & connectivity <= 27)

  mask.idx <- which(mask > 0)
  grid <- indexToCoord(mask, mask.idx)
  vgrid <- indexToGrid(mask, mask.idx)

  #if (shrink > 0) {
  #  message("running ", shrink, "shrinkage iterations with k = ", 5)
  #  bvec <- shrink_vector(mask, vgrid, bvec, k=5, iter=shrink)
  #}

  kres <- kmeans(grid, K, iter.max=500)
  kvol <- NeuroVol(kres$cluster, space(mask), indices=mask.idx)

  clusid <- sort(unique(kres$cluster))
  neib <- get.knn(grid, k=connectivity)



  iter <- 1
  switches <- 1
  iter.max <- iterations

  valmat <- series(bvec, mask.idx)

  centroids <- compute_centroids(valmat, grid, kres$cluster, medoid=use_medoid)
  sp_centroids <- do.call(rbind, lapply(centroids, "[[", "centroid"))
  num_centroids <- do.call(rbind, lapply(centroids, "[[", "center"))

  denom <- max(get.knn(sp_centroids, k=1)$nn.dist[,1])


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



    centroids <- compute_centroids(bvec, mask.idx, grid, newclus, medoid=use_medoid)
    switches <- sum(newclus != curclus)

    curclus <- newclus
    message("tuboclust: nswitches ", switches)

    iter <- iter + 1
  }

  kvol <- NeuroVol(newclus, space(mask), indices=mask.idx)
  voxsplit <- split(data.frame(vgrid), newclus)

  list(clusvol=kvol,
       clusters=newclus,
       centers = do.call(rbind, lapply(centroids, "[[", "center")),
       coordCenters=do.call(rbind, lapply(centroids, "[[", "centroid")),
       voxelSets=voxsplit)
}


#' @inheritParams turbo_cluster
turbo_cluster_fit <- function(valmat, grid, sigma1=1, sigma2=10, K=min(500, nrow(grid)),
                              iterations=25, connectivity=27, use_medoid=FALSE, initclus=NULL) {

  assert_that(connectivity > 1 & connectivity <= 27)

  if (is.null(initclus)) {
    ## kmeans using coordinates only
    kres <- kmeans(grid, K, iter.max=500)
    clusid <- sort(unique(kres$cluster))
    curclus <- kres$cluster
  } else {
    assert_that(length(initclus) == nrow(grid))
    clusid <- sort(unique(initclus))
    assert_that(length(clusid) == K)
    curclus <- initclus
  }

  ## find k neighbors within 'connectivity' radius
  neib <- FNN::get.knn(grid, k=connectivity)

  ## has to be changed for surface... pass in neighbor_fun?
  dthresh <- median(neib$nn.dist[,connectivity,drop=FALSE])

  centroids <- compute_centroids(valmat, grid, curclus, medoid=use_medoid)

  sp_centroids <- do.call(rbind, lapply(centroids, "[[", "centroid"))
  num_centroids <- do.call(rbind, lapply(centroids, "[[", "center"))

  iter <- 1
  switches <- 1
  iter.max <- iterations

  while (iter < iter.max && switches > 0) {
    candlist <- find_candidates(neib$nn.index-1, neib$nn.dist, curclus, dthresh)

    newclus <- best_candidate(candlist, curclus, t(grid),
                              t(num_centroids), t(sp_centroids),
                              valmat, sigma1, sigma2)

    switches <- attr(newclus, "nswitches")
    if (switches > 0) {
      centroids <- compute_centroids(valmat, grid, newclus, medoid=use_medoid)
      sp_centroids <- do.call(rbind, lapply(centroids, "[[", "centroid"))
      num_centroids <- do.call(rbind, lapply(centroids, "[[", "center"))
      curclus <- newclus
    }

    message("turbo_clust_fit: iter ", iter, " -- num  switches =  ", switches)
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

  valmat <- series(bvec, mask.idx)

  ## find k neighbors within 'connectivity' radius
  neib <- FNN::get.knn(grid, k=connectivity)

  svalmat <- t(do.call(rbind, map(1:nrow(neib$nn.index), function(i) {
    rowMeans(valmat[, c(i, neib$nn.index[i,1:(k-1)])])
  })))

  NeuroVec(valmat, space(bvec), mask=mask)
}

#' turbo_cluster
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
#' @param filter low- and high-pass filter parameters. See details.
#' @export
#' @importFrom neuroim2 NeuroVec NeuroVol
#' @import furrr
#'
#' @return a \code{list} of class \code{turbo_cluster_result} with the following elements:
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
#' cres1 <- turbo_cluster(bvec, mask, K=100, sigma1=1, sigma2=10, sample_frac=.3)
#' cres2 <- turbo_cluster(bvec, mask, K=100, sigma1=1, sigma2=6, sample_frac=.3)
#' cres3 <- turbo_cluster(bvec, mask, K=100, sigma1=1, sigma2=4, sample_frac=.3)
#' cres4 <- turbo_cluster(bvec, mask, K=100, sigma1=1, sigma2=2, sample_frac=.3)
#'
#' cres_cons <- merge_clus(cres1, cres2, cres3, cres4)
#'
#' cres_cons2 <- turbo_cluster(bvec, mask, K=100, sigma1=c(1,2,3), sigma2=c(3,2,1), sample_frac=.3)
#'
#' ## to access the cluster volume: cres1$clusvol
turbo_cluster <- function(bvec, mask, K=500, sigma1=1,
                                sigma2=2.5, iterations=50,
                                connectivity=27, use_medoid=FALSE,
                                filter=c(lp=0, hp=0),
                                sample_frac=1, nreps=1) {
  mask.idx <- which(mask > 0)
  grid <- index_to_coord(mask, mask.idx)

  valmat <- series(bvec, mask.idx)

  if (any(filter > 0)) {
    assert_that(filter$lp <= 1 && filter$hp <= 1)
    message("turbo_cluster: pre-filtering time series data")
    valmat <- filter_mat(valmat, filter$lp, filter$hp)
  }


  if (length(sigma1) != nreps) {
    sigma1 <- rep(sigma1, length.out=nreps)
  }

  if (length(sigma2) != nreps) {
    sigma2 <- rep(sigma2, length.out=nreps)
  }

  res <- furrr::future_map(1:nreps, function(i) {
    if (sample_frac < 1) {
      n <- max(sample_frac * nrow(valmat),1)
      assert_that(n > 1)
      keep <- sample(1:nrow(valmat), n)
      vmat <- valmat[keep,,drop=FALSE]
    }

    ret <- turbo_cluster_fit(vmat, grid, sigma1=sigma1[i], sigma2=sigma2[i], K=K,
                                iterations=iterations, connectivity=connectivity,
                                 use_medoid=use_medoid)
    kvol <- ClusteredNeuroVol(as.logical(mask), clusters=ret$clusters)

    structure(
         list(clusvol=kvol,
              cluster=ret$clusters,
              coord_centers=ret$coord_centers),
         class=c("cluster_result", "list"))
  })

  if (length(res) > 1) {
    message("computing consensus clustering over ", length(res), "partitions")
    kvol <- do.call(merge_clus, res)
    cen <- compute_centroids(valmat, grid,assignment=kvol@clusters, medoid=use_medoid)
    structure(
      list(clusvol=kvol,
           cluster=kvol@clusters,
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
  cons <- clue::cl_consensus(ens)
  hpart <- clue::as.cl_hard_partition(cons)
  hpart$.Data
  #ClusteredNeuroVol(x$clusvol@mask, clusters=hpart$.Data)
}

merge_clus.cluster_result_time <- function(x, ...) {
  args <- c(list(x), list(...))
  assert_that(all(map_lgl(args, ~ inherits(., "cluster_result"))))

  ens <- do.call(clue::cl_ensemble, args)
  cons <- clue::cl_consensus(ens)
  hpart <- clue::as.cl_hard_partition(cons)
  hpart$.Data
}

#' @export
cl_class_ids.cluster_result <- function(x) {
  x$cluster
}

is.cl_partition.cluster_result <- function(x) {
  TRUE
}

#' @param valmat a feature matrix with \code{nrow(valmat)} observations and `\code{ncol(valmat)} features.
#' @param lp the low-pass bandwidth parameter specified as a fraction
#'    of the number of time-points to include in each window. Must be between [0-1]
#'  @param hp the high-pass bandwidth parameter specified as a fraction
#'    of the number of time-points to include in each window. Must be between [0-1]
#' @importFrom locfit locfit lp
#' @examples
#'
#' valmat <- matrix(rnorm(200*10), 200, 10)
#' library(future)
#' plan(multicore)
#' svalmat <- filter_mat(valmat, lp=.2, hp=.7)
#' svalmat2 <- filter_mat(valmat, lp=.2, hp=.7, method="bspline")
#'
filter_mat <- function(valmat, lp=0, hp=.7, method=c("locfit", "bspline")) {
  assert_that(hp >= lp)
  assert_that(hp >= 0 && hp <= 1)
  assert_that(lp >= 0 && lp <= 1)

  method <- match.arg(method)

  if (method == "locfit") {
    nn <- hp
    nn1 <- lp
    time <- seq(1,nrow(valmat))

    if (lp ==0 && hp == 0) {
      warning("filter_mat: returning oriignal matrix unchanged, both 'lp' and 'hp' are 0")
      return(valmat)
    }

    fvalmat <- do.call(cbind, furrr::future_map(1:ncol(valmat), function(i) {
      vals <- valmat[,i]
      fvals <- vals - fitted(locfit(vals ~ lp(time, nn=hp)))

      if (lp > 0) {
        fvals <- fitted(locfit(fvals ~ lp(time, nn=lp)))
      }
      fvals
    }))
  } else {

    fvalmat <- if (hp > 0) {
      nhp <- ceiling(nrow(valmat)/(hp * nrow(valmat)))
      bmat1 <- if (nhp < 2) {
        splines::ns(1:nrow(valmat), nhp)
      } else {
        splines::bs(1:nrow(valmat), nhp)
      }

      resid(lsfit(bmat1, valmat))
    } else {
      valmat
    }

    if (lp > 0) {
      nlp <- nrow(valmat)/(lp * nrow(valmat))
      nlp <- min(nlp, nrow(valmat)-1)

      bmat2 <- splines::bs(1:nrow(valmat), nlp)
      fvalmat <- fitted(lm(fvalmat ~ bmat2))
    }
  }

  fvalmat

}

#' @export
#' @inheritParam turbo_cluster
#' @examples
#' valmat <- matrix(rnorm(100*10), 100, 10)
#' library(future)
#' plan(multicore)
#' cres <- turbo_cluster_time(t(valmat), K=20)
#'
turbo_cluster_time <- function(valmat, K=min(nrow(valmat), 100), sigma1=1, sigma2=TR*3,
                               iterations=50, TR=2, filter=list(lp=0, hp=0),
                               use_medoid=FALSE, nreps=5) {
  if (filter$lp > 0 || filter$hp > 0) {
    message("turbo_cluster: filtering time series")
    valmat <- filter_mat(valmat, filter$lp, filter$hp)
  }

  nels <- nrow(valmat)
  grid <- as.matrix(seq(0, by=TR, length.out=nels))

  fits <- map(1:nreps, function(i) {
    initsamps <- sort(sample(1:nrow(valmat), K))
    initcenters <- valmat[initsamps,]

    curclus <- FNN::get.knnx(grid[initsamps,,drop=FALSE], grid, k=1)$nn.index[,1]

    ret <- turbo_cluster_fit(t(valmat), as.matrix(grid),
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
turbo_cluster_surface <- function(bsurf, K=500, sigma1=1, sigma2=5, iterations=50,
                                  connectivity=6, use_medoid=FALSE, filter=list(lp=-1, hp=-1)) {
  mask.idx <- indices(bsurf)
  grid <- coords(bsurf)[mask.idx,]

  valmat <- series(bsurf, mask.idx)

  if (any(filter) > 0) {
    message("turbo_cluster: filtering time series")
    valmat <- filter_mat(valmat, filter$lp, filter$hp)
  }

  ret <- turbo_cluster_fit(valmat, grid, sigma1=sigma1, sigma2=sigma2, K=K,
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



#' @export
meta_clust.turbo_cluster_surface <- function(x, cuts) {
  D <- 1 - cor(t(x$centers))
  hres <- hclust(as.dist(D), method="ward.D2")

  orig <- x$clusters

  cmat <- do.call(cbind, lapply(cuts, function(i) {
    cind <- cutree(hres, i)
    cind[orig]
  }))

  out <- cbind(cmat, orig)
  BrainSurfaceVector(geometry(x$clusvol), indices=x$clusvol@indices, mat=as.matrix(out))

}

