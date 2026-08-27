# Initialization Methods for 4D Clustering

Functions for initializing cluster seeds using various strategies.

## Usage

``` r
initialize_clusters(
  coords,
  features,
  n_clusters,
  method = c("gradient", "kmeans", "poisson", "grid", "random"),
  mask = NULL,
  vec = NULL,
  ...
)
```
