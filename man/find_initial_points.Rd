% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snic.R
\name{find_initial_points}
\alias{find_initial_points}
\title{Find Initial Cluster Centers for Supervoxel Algorithm}
\usage{
find_initial_points(cds, grad, K = 100)
}
\arguments{
\item{cds}{A matrix or data frame representing the spatial coordinates of the voxels.}

\item{grad}{A vector representing the gradient values of the voxels.}

\item{K}{The desired number of supervoxels (clusters) in the output (default: 100).}
}
\value{
A list containing two elements:
        selected - a vector of the selected indices corresponding to the initial cluster centers,
        coords - a matrix or data frame with the spatial coordinates of the initial cluster centers.
}
\description{
This function finds the initial cluster centers for a supervoxel algorithm.
Supervoxels are used to partition 3D image data into volumetric regions,
grouping similar voxels together. The initial cluster centers are crucial for
the performance and quality of the final supervoxels.
}
