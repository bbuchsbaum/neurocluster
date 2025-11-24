# Prune 1-NN edges to achieve target number of components

Sorts 1-NN edges by distance and adds them incrementally until K
components remain. Used for exact-K stopping condition.

## Usage

``` r
prune_edges_for_k_cpp(n_nodes, nearest_neighbor, distances, target_k)
```

## Arguments

- n_nodes:

  Number of nodes

- nearest_neighbor:

  IntegerVector of nearest neighbor for each node (0-based)

- distances:

  NumericVector of distances to nearest neighbors

- target_k:

  Target number of components

## Value

IntegerVector of pruned nearest neighbors (-1 for pruned edges)
