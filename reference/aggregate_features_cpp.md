# Aggregate features by component (mean pooling)

Computes mean feature vector for each component efficiently.

## Usage

``` r
aggregate_features_cpp(feature_mat, component_labels, n_components)
```

## Arguments

- feature_mat:

  Numeric matrix (features x voxels)

- component_labels:

  IntegerVector of component labels (0-based)

- n_components:

  Number of unique components

## Value

NumericMatrix of aggregated features (features x n_components)
