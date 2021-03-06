% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_centroids.R
\name{compute_centroids}
\alias{compute_centroids}
\title{compute_centroids}
\usage{
compute_centroids(feature_mat, grid, assignment, medoid = FALSE)
}
\arguments{
\item{feature_mat}{a matrix for features where each column is a feature and each row is an observation}

\item{grid}{a matrix of coordinates, where each row is a coordinates associated with the nth feature.}

\item{assignment}{a vector of integer indices indicating the cluster assignment}

\item{medoid}{compute the medoid rather than the mean of the cluster}
}
\description{
Get centroids for a matrix and set of assignments
}
\keyword{internal}
