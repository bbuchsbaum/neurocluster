# Test if Object is a Partition

This function tests whether a cluster result object represents a
partition. It is a method for the
[`is.cl_partition`](https://rdrr.io/pkg/clue/man/partition.html) generic
from the `clue` package. For `cluster_result` objects, this always
returns `TRUE` since cluster results represent valid partitions where
each data point belongs to exactly one cluster.

## Usage

``` r
# S3 method for class 'cluster_result'
is.cl_partition(x)
```

## Arguments

- x:

  A `cluster_result` object.

## Value

`TRUE`, indicating that cluster results are always valid partitions.

## See also

[`is.cl_partition`](https://rdrr.io/pkg/clue/man/partition.html) for the
generic function.
