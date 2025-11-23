# Merge Clustering Results Using a Consensus Clustering Algorithm

The merge_clus function combines a set of clustering results using a
consensus clustering algorithm.

## Usage

``` r
merge_clus(x, method, ...)
```

## Arguments

- x:

  A clustering result, typically a list or an object of class
  `"cluster_result"`.

- method:

  A character string indicating the consensus clustering algorithm to
  use. Default is "SE". See
  [`cl_consensus`](https://rdrr.io/pkg/clue/man/cl_consensus.html) for
  available methods.

- ...:

  Additional clustering results to be merged.

## Value

A
[`ClusteredNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.html)
instance.

## See also

[`cl_consensus`](https://rdrr.io/pkg/clue/man/cl_consensus.html),
[`as.cl_hard_partition`](https://rdrr.io/pkg/clue/man/partition.html),
[`cl_ensemble`](https://rdrr.io/pkg/clue/man/cl_ensemble.html)

## Examples

``` r
# Assuming clustering1, clustering2, and clustering3 are objects of class "cluster_result"
merged_clustering <- merge_clus(clustering1, clustering2, clustering3, method="SE")
#> Error in UseMethod("merge_clus", x, method): unused argument (method)
```
