# Aggregate coordinates by component (mean pooling)

Computes mean coordinate for each component efficiently.

## Usage

``` r
aggregate_coords_cpp(coords, component_labels, n_components)
```

## Arguments

- coords:

  Numeric matrix (voxels x 3)

- component_labels:

  IntegerVector of component labels (0-based)

- n_components:

  Number of unique components

## Value

NumericMatrix of aggregated coordinates (n_components x 3)
