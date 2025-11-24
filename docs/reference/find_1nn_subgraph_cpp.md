# Find 1-Nearest Neighbor subgraph for ReNA

For each node, finds its single nearest neighbor to form directed 1-NN
graph.

## Usage

``` r
find_1nn_subgraph_cpp(n_nodes, adjacency_i, adjacency_j, distances)
```

## Arguments

- n_nodes:

  Number of nodes

- adjacency_i:

  Integer vector of row indices for sparse adjacency

- adjacency_j:

  Integer vector of col indices for sparse adjacency

- distances:

  Numeric vector of distances for each edge

## Value

IntegerVector of nearest neighbor indices (0-based, -1 if no neighbors)
