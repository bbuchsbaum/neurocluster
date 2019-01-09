#' @export
meta_clust.cluster_result <- function(x, cuts=min(as.integer(length(x$centers)/2),2),
                                      algo="hclust", hclust_method="ward.D") {

  orig <- x$cluster

  cvols <- if (algo == "hclust") {
    cen <- do.call(rbind, x$centers)
    #D <- 1 - cor(t(x$centers))

    D <- 1 - cor(t(cen))
    hres <- hclust(as.dist(D), method=hclust_method)
    cmat <- do.call(cbind, lapply(cuts, function(i) {
      cind <- cutree(hres, i)
      cind[orig]
    }))

    cvols <- map(1:ncol(cmat), function(i) {
      ClusteredNeuroVol(x$clusvol@mask, cluster=cmat[,i])
    })

    c(cvols, list(x$clusvol))
  } else {
    stop()
  }
  # } else if (algo == "apcluster") {
  #   S <- cor(t(x$centers))
  #   apres <- apcluster(S, details=TRUE)
  #   out <- integer(nrow(x$centers))
  #   clmap <- do.call(rbind, imap(apres@clusters, function(cl,i) {
  #     cbind(cl,i)
  #   }))
  #
  #   cltab <- as.list(clmap[,2])
  #   names(cltab) <- as.character(clmap[,1])
  #   clusters <- unlist(cltab[as.character(orig)])
  # }

  list(cvols=cvols, cuts=cuts, cutmat=cmat, hclus=hres)

}


# @export
# meta_clust.turbo_cluster_surface <- function(x, cuts) {
#   D <- 1 - cor(t(x$centers))
#   hres <- hclust(as.dist(D), method="ward.D2")
#
#   orig <- x$clusters
#
#   cmat <- do.call(cbind, lapply(cuts, function(i) {
#     cind <- cutree(hres, i)
#     cind[orig]
#   }))
#
#   out <- cbind(cmat, orig)
#   BrainSurfaceVector(geometry(x$clusvol), indices=x$clusvol@indices, mat=as.matrix(out))
#
# }
