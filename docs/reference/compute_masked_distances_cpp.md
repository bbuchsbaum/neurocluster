# Compute masked distances for ReNA

Computes squared Euclidean distances only for connected pairs in sparse
adjacency.

## Usage

``` r
compute_masked_distances_cpp(feature_mat, adjacency_i, adjacency_j)
```

## Arguments

- feature_mat:

  Numeric matrix (features x voxels)

- adjacency_i:

  Integer vector of row indices for sparse adjacency

- adjacency_j:

  Integer vector of col indices for sparse adjacency

## Value

NumericVector of distances (same length as adjacency_i)
