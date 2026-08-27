# Common Utilities for 4D Clustering Algorithms

Internal functions shared across all cluster4d methods for consistency
and code reuse.

Internal, fail-closed capability table used to reject common parameters
that a method does not implement. A NULL connectivity set means the
unified connectivity argument is inactive for that method.

## Usage

``` r
ensure_neurovec(vec)

cluster4d_method_contract(method)
```

## Arguments

- method:

  Unified cluster4d method name.

## Value

A list describing supported common parameters.
