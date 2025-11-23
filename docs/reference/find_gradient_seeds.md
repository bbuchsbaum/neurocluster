# Find gradient-based seed points

Internal function that implements the gradient-based seed selection
algorithm (previously find_initial_points).

## Usage

``` r
find_gradient_seeds(coords, grad_vals, K, min_separation_factor = 1.5)
```

## Arguments

- coords:

  Coordinate matrix

- grad_vals:

  Gradient values at each voxel

- K:

  Number of seeds to find

## Value

Vector of selected seed indices
