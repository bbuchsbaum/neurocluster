# K-means based seed initialization

Uses k-means clustering to find initial seeds.

## Usage

``` r
init_kmeans_seeds(coords, features = NULL, n_clusters, use_features = FALSE)
```

## Arguments

- coords:

  Coordinate matrix

- features:

  Feature matrix (optional)

- n_clusters:

  Number of clusters

- use_features:

  Whether to use features in k-means

## Value

List with seeds, initial_labels, and seed_coords
