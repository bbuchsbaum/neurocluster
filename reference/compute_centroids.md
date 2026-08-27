# Compute Centroids of Clusters

The compute_centroids function calculates the center and centroid of
each cluster given a feature matrix, grid, and cluster assignments.

## Usage

``` r
compute_centroids(feature_mat, grid, assignment, medoid = FALSE)
```

## Arguments

- feature_mat:

  A matrix of features, where columns represent data points and rows
  represent features.

- grid:

  A matrix representing the spatial grid of data points.

- assignment:

  A vector containing the cluster assignment for each data point.

- medoid:

  A logical value indicating whether to calculate medoids instead of
  means for cluster centers and centroids. Default is FALSE.

## Value

A list containing two elements:

- center:

  A matrix containing the centers of each cluster.

- centroid:

  A matrix containing the centroids of each cluster.

## Examples

``` r
# Simple synthetic example
feature_mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
grid <- cbind(x = runif(10), y = runif(10), z = runif(10))
assignment <- rep(1:2, each = 5)
centroids <- compute_centroids(feature_mat, grid, assignment)
print(names(centroids))
#> [1] "center"   "centroid"

if (FALSE) { # \dontrun{
  # Larger example with real neuroimaging data
  # Assuming `feature_mat`, `grid`, and `assignment` are available
  centroids <- compute_centroids(feature_mat, grid, assignment)
  # To compute medoids instead of means
  medoids <- compute_centroids(feature_mat, grid, assignment, medoid=TRUE)
} # }
```
