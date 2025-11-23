# Extract Class IDs from Cluster Result

This function extracts the cluster class identifiers from a cluster
result object. It is a method for the
[`cl_class_ids`](https://rdrr.io/pkg/clue/man/n_of_classes.html) generic
from the `clue` package.

## Usage

``` r
# S3 method for class 'cluster_result'
cl_class_ids(x)
```

## Arguments

- x:

  A `cluster_result` object containing clustering information.

## Value

An integer vector of cluster assignments, one for each data point.

## See also

[`cl_class_ids`](https://rdrr.io/pkg/clue/man/n_of_classes.html) for the
generic function.
