# Gradient-based seed initialization

Finds initial cluster centers using spatial gradient information. Seeds
are placed at locations with high gradient and good spatial separation.

## Usage

``` r
init_gradient_seeds(
  coords,
  n_clusters,
  mask = NULL,
  vec = NULL,
  gradient_type = "spatial"
)
```

## Arguments

- coords:

  Coordinate matrix (n_voxels x 3)

- n_clusters:

  Number of seeds to find

- mask:

  NeuroVol mask (required for gradient)

- vec:

  NeuroVec data (required for gradient)

- gradient_type:

  Type of gradient: "spatial", "correlation", "intensity"

## Value

List with seeds, initial_labels, and seed_coords
