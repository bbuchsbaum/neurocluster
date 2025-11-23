# Poisson disk sampling for seed initialization

Places seeds with guaranteed minimum separation distance.

## Usage

``` r
init_poisson_seeds(coords, n_clusters, min_distance = NULL)
```

## Arguments

- coords:

  Coordinate matrix

- n_clusters:

  Approximate number of clusters

- min_distance:

  Minimum distance between seeds (auto-computed if NULL)

## Value

List with seeds, initial_labels, and seed_coords
