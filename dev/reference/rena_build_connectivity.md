# Build connectivity graph for ReNA

Creates exact sparse masked-grid adjacency.

## Usage

``` r
rena_build_connectivity(mask, mask_idx, connectivity)
```

## Arguments

- mask:

  A three-dimensional `NeuroVol` mask.

- mask_idx:

  Included linear indices in canonical mask order.

- connectivity:

  Grid connectivity: 6, 18, or 26.

## Value

A symmetric sparse adjacency matrix in masked-voxel order.
