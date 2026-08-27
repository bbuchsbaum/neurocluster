# Spatial Gradient Calculation

The spatial_gradient function calculates the spatial gradient of a
`NeuroVol` instance within the specified mask.

## Usage

``` r
spatial_gradient(vol, mask, sigma = 0.5)
```

## Arguments

- vol:

  A `NeuroVol` instance for which the spatial gradient should be
  calculated.

- mask:

  A `NeuroVol` mask defining the voxels to include in the spatial
  gradient calculation. If the mask contains `numeric` data, nonzero
  values will define the included voxels. If the mask is a
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.html),
  then `TRUE` will define the set of included voxels.

- sigma:

  A numeric value controlling the spatial weighting function. Default is
  0.5.

## Value

A `NeuroVol` instance containing the spatial gradient values for the
input `vol`.

## See also

[`spatial_laplacian`](https://rdrr.io/pkg/adjoin/man/spatial_laplacian.html),
[`weighted_spatial_adjacency`](https://rdrr.io/pkg/adjoin/man/weighted_spatial_adjacency.html)

## Examples

``` r
if (FALSE) { # \dontrun{
mask <- neuroim2::NeuroVol(array(1, c(20,20,20)), neuroim2::NeuroSpace(c(20,20,20)))
input_vol <- neuroim2::NeuroVol(array(runif(20*20*20), c(20,20,20)),
  neuroim2::NeuroSpace(c(20,20,20)))
gradient_vol <- spatial_gradient(input_vol, mask)
} # }
```
