# Merge Clustering Results for ClusteredNeuroVol Objects

This method of merge_clus is specifically designed to merge clustering
results for `ClusteredNeuroVol` objects by performing consensus
clustering across multiple clustering results.

## Usage

``` r
# S3 method for class 'cluster_result'
merge_clus(x, method = "SE", ...)
```

## Arguments

- x:

  A `ClusteredNeuroVol` object or an object of class `"cluster_result"`.

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
instance representing the consensus partition.

## See also

[`cl_consensus`](https://rdrr.io/pkg/clue/man/cl_consensus.html),
[`as.cl_hard_partition`](https://rdrr.io/pkg/clue/man/partition.html),
[`cl_ensemble`](https://rdrr.io/pkg/clue/man/cl_ensemble.html)
