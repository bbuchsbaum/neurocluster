# Find connected components using Union-Find

Finds weakly connected components in a directed graph defined by 1-NN
edges.

## Usage

``` r
find_connected_components_cpp(n_nodes, nearest_neighbor)
```

## Arguments

- n_nodes:

  Number of nodes

- nearest_neighbor:

  IntegerVector of nearest neighbor for each node (0-based)

## Value

IntegerVector of component labels (0-based, contiguous)
