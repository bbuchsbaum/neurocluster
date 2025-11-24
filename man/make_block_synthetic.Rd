\name{make_block_synthetic}
\alias{make_block_synthetic}
\title{Make a tiny block-structured synthetic volume}
\description{
A scikit-learn style "make\_blobs" for neuroimaging demos. Builds a small
2D/3D mask with three spatial bands (left/mid/right) and clearly separated
time-series patterns plus mild Gaussian noise. Designed to be fast,
deterministic, and visually intuitive in vignettes and tests.
}
\usage{
make_block_synthetic(dims = c(12, 12, 1), ntime = 80, noise = 0.15, seed = NULL)
}
\arguments{
\item{dims}{Integer vector length 3, spatial dimensions. Default \code{c(12,12,1)}.}

\item{ntime}{Number of time points. Default 80.}

\item{noise}{Standard deviation of added Gaussian noise. Default 0.15.}

\item{seed}{Optional integer seed for reproducibility.}
}
\value{
List with elements \code{vec} (\code{NeuroVec}), \code{mask} (\code{NeuroVol}), and integer \code{truth} labels.
}
\examples{
syn <- make_block_synthetic(dims = c(12, 12, 1), ntime = 60, noise = 0.1, seed = 7)
length(unique(syn$truth))
}
