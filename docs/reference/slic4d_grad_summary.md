# Visualize gradient used for seed relocation

Compute and return the gradient volume that would be used for seed
relocation in slic4d_supervoxels. Useful for debugging and
visualization.

## Usage

``` r
slic4d_grad_summary(
  bvec,
  mask,
  method = c("correlation", "intensity", "spatial")
)
```

## Arguments

- bvec:

  A `NeuroVec` with dims (X, Y, Z, T).

- mask:

  A 3D `NeuroVol` (or logical array) indicating voxels to include.

- method:

  One of "correlation", "intensity", "spatial".

## Value

A 3D array containing the gradient values.
