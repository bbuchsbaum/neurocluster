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

`spatial_laplacian`, `weighted_spatial_adjacency`

## Examples

``` r
mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
#> Error in NeuroVol(array(1, c(20, 20, 20)), NeuroSpace(c(20, 20, 20))): could not find function "NeuroVol"
input_vol <- NeuroVol(array(runif(202020), c(20,20,20)),
NeuroSpace(c(20,20,20)))
#> Error in NeuroVol(array(runif(202020), c(20, 20, 20)), NeuroSpace(c(20,     20, 20))): could not find function "NeuroVol"

gradient_vol <- spatial_gradient(input_vol, mask)
#> Error: object 'mask' not found
```
