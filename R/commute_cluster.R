

#' @inheritParams turbo_cluster
#'
#' @import neighborweights
#' @export
commute_cluster <- function(bvec, mask,
                            K=100,
                            ncomp=ceiling(sqrt(K*2)),
                            alpha=.2,
                            sigma1=.73,
                            sigma2=5,
                            connectivity=27,
                            filter=c(lp=0, hp=0),
                            filter_method=c("bspline", "locfit"),
                            weight_mode=c("binary", "heat")) {

  filter_method <- match.arg(filter_method)
  mask.idx <- which(mask > 0)
  grid <- index_to_coord(mask, mask.idx)

  feature_mat <- neuroim2::series(bvec, mask.idx)

  if (any(filter > 0)) {
    assert_that(filter$lp <= 1 && filter$hp <= 1)
    message("commute_cluster: pre-filtering time series data with ", filter_method, " smoothing")
    feature_mat <- filter_mat(feature_mat, filter$lp, filter$hp, method=filter_method)
  }

  feature_mat <- scale(feature_mat)

  message("commute_cluster: computing similarity matrix.")
  W <- neighborweights::weighted_spatial_adjacency(grid, t(feature_mat),
                                                   dthresh=sigma2*2,
                                                   nnk=connectivity,
                                                   wsigma=sigma1,
                                                   sigma=sigma2,
                                                   alpha=alpha,
                                                   weight_mode=weight_mode,
                                                   include_diagonal=FALSE,
                                                   stochastic=TRUE)

  message("commute_cluster: computing commute-time embedding.")
  ct <- neighborweights::commute_time(W, ncomp=ncomp)

  message("commute_cluster: performing k-means with ", K, " clusters.")
  kres <- kmeans(ct$cds, center=K, iter.max=100)
  kvol <- ClusteredNeuroVol(as.logical(mask),kres$cluster)

  message("commute_cluster: computing final centroids")
  centroids <- compute_centroids(feature_mat, grid, kres$cluster, medoid=FALSE)

  ret <- structure(
    list(clusvol=kvol,
         cluster=kres$cluster,
         centers=centroids$center,
         coord_centers=centroids$centroid),
    class=c("commute_time_cluster_result", "cluster_result", "list"))

  ret
}

