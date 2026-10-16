# Contract adjacency graph by merging components

Builds new adjacency matrix where nodes are components and edges exist
if any constituent nodes were connected.

## Usage

``` r
contract_graph_cpp(adjacency_i, adjacency_j, component_labels, n_components)
```

## Arguments

- adjacency_i:

  Integer vector of row indices for sparse adjacency

- adjacency_j:

  Integer vector of col indices for sparse adjacency

- component_labels:

  IntegerVector of component labels (0-based)

- n_components:

  Number of unique components

## Value

List with contracted_i and contracted_j vectors
