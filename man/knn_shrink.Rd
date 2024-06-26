% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supervoxels.R
\name{knn_shrink}
\alias{knn_shrink}
\title{knn_shrink}
\usage{
knn_shrink(bvec, mask, k = 5, connectivity = 27)
}
\arguments{
\item{bvec}{bvec a \code{\linkS4class{NeuroVec}} instance}

\item{mask}{a mask defining the voxels to include in the clustering result.}

\item{k}{the number of nearest neighbors to average over.}

\item{connectivity}{the number of spatial neighbors to include in the search space around each voxel.}
}
\description{
Replace each voxel by the mean of its k nearest neighbors in its local spatial neighborhood
}
\examples{
mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
bvec <- replicate(10, NeuroVol(array(runif(20*20*20), c(20,20,20)),
NeuroSpace(c(20,20,20))), simplify=FALSE)
bvec <- do.call(concat, bvec)

sbvec <- knn_shrink(bvec, mask,k=3)

}
