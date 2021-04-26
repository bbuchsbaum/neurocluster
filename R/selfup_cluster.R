
hkernel <- function(d, T) {
  exp(-d/T)
}


#' @inheritParams turbo_cluster
#'
#' @import neighborweights
#' @export
selfup_cluster <- function(bvec, mask,
                            alpha=.2,
                            spatial_radius=8,
                            feature_radius=1,
                            sigma1=feature_radius/5,
                            sigma2=spatial_radius/5,
                            heating_rate=.01,
                            sigma2_t0=sigma2/5,
                            sigma1_t0=sigma1,
                            filter=c(lp=0, hp=0),
                            filter_method=c("bspline", "locfit")) {


  filter_method <- match.arg(filter_method)
  mask.idx <- which(mask > 0)
  grid <- index_to_coord(mask, mask.idx)

  feature_mat <- neuroim2::series(bvec, mask.idx)
  csds <- matrixStats::colSds(feature_mat)
  bad <- which(is.na(csds) | csds==0)

  if (length(bad) > 0) {
    warning(paste(length(bad), "have time-courses have NA or 0 stdev. Random vectors"))
    feature_mat[, bad] <- matrix(rnorm(length(bad)*nrow(feature_mat)), nrow(feature_mat), length(bad))
  }

  if (any(filter > 0)) {
    assert_that(filter$lp <= 1 && filter$hp <= 1)
    message("selfup_cluster: pre-filtering time series data with ", filter_method, " smoothing")
    feature_mat <- filter_mat(feature_mat, filter$lp, filter$hp, method=filter_method)
  }

  feature_mat <- scale(feature_mat)

  message("selfup_cluster: computing similarity matrix.")

  converged <- FALSE
  iter=1

  sigma1_cur=sigma1_t0
  sigma2_cur=sigma2_t0


  X <- feature_mat
  cds <- grid

  snn <- rflann::RadiusSearch(cds, cds, spatial_radius^2, max_neighbour=200)

  while (!converged) {

    message("sigma1: ", sigma1_cur)
    message("sigma2: ", sigma2_cur)


    ret <- do.call(rbind, lapply(seq_along(snn$distances), function(i) {

      d1 <- 1 - cor(X[,i], X[,snn$indices[[i]]])
      d1[d1 > feature_radius] <- Inf
      s1<- hkernel(d1, sigma1_cur)
      wts1 <- as.vector(s1/sum(s1))

      #s2 <- hkernel(sqrt(snn$distances[[i]]), sigma2_cur)
      #wts2 <- s2/sum(s2)
      #wts <- (wts1 + wts2)/2
      cbind(i, snn$indices[[i]], wts1)
    }))

    X_old <- X
    #cds_old <- cds

    M <- Matrix::sparseMatrix(i=ret[,1], j=ret[,2], x=ret[,3])
    #cds <- as.matrix(t(t(cds) %*% M))
    X <- t(as.matrix(M %*% t(X)))


    #delta1 <- sum(abs(cds - cds_old))
    delta2 <- sum(abs(X - X_old))
    #message("cds diff:", delta1)
    message("X diff:", delta2)

    sigma1_cur <- min(sigma1_cur + sigma1_cur * heating_rate, sigma1)
    sigma2_cur <- min(sigma2_cur + sigma2_cur * heating_rate, sigma2)

    if (delta2 < 1) {
      break
    }
  }

}
