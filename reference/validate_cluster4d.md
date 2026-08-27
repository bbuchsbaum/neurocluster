# Validate cluster4d result

Checks validity and quality of clustering results.

## Usage

``` r
validate_cluster4d(result, vec = NULL, mask = NULL, tolerance = 1e-08)
```

## Arguments

- result:

  A cluster4d_result object

- vec:

  Original NeuroVec data, required to verify feature centers.

- mask:

  Original mask, required to verify geometry, voxel order, and physical
  coordinate centers.

- tolerance:

  Relative numeric tolerance used when independently recomputing centers
  and coordinate centers.

## Value

A list with `valid`, `errors`, `warnings`, and a compact `summary`. Any
schema or consistency defect sets `valid = FALSE`.
