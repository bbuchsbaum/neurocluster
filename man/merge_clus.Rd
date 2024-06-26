% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allgeneric.R
\name{merge_clus}
\alias{merge_clus}
\title{Merge Clustering Results Using a Consensus Clustering Algorithm}
\usage{
merge_clus(x, method, ...)
}
\arguments{
\item{x}{A clustering result, typically a list or an object of class \code{"cluster_result"}.}

\item{method}{A character string indicating the consensus clustering algorithm to use. Default is "SE".
See \code{\link[clue]{cl_consensus}} for available methods.}

\item{...}{Additional clustering results to be merged.}
}
\value{
A \code{\linkS4class{ClusteredNeuroVol}} instance.
}
\description{
The merge_clus function combines a set of clustering results using a consensus clustering algorithm.
}
\examples{
# Assuming clustering1, clustering2, and clustering3 are objects of class "cluster_result"
merged_clustering <- merge_clus(clustering1, clustering2, clustering3, method="SE")

}
\seealso{
\code{\link[clue]{cl_consensus}}, \code{\link[clue]{as.cl_hard_partition}}, \code{\link[clue]{cl_ensemble}}
}
