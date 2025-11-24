# Make a tiny block-structured synthetic volume

A scikit-learn style "make\\blobs" for neuroimaging demos. Builds a
small 2D/3D mask with three spatial bands (left/mid/right) and clearly
separated time-series patterns plus mild Gaussian noise. Designed to be
fast, deterministic, and visually intuitive in vignettes and tests.

## Usage

``` r
make_block_synthetic(dims = c(12, 12, 1), ntime = 80, noise = 0.15, seed = NULL)
```

## Arguments

- dims:

  Integer vector length 3, spatial dimensions. Default `c(12,12,1)`.

- ntime:

  Number of time points. Default 80.

- noise:

  Standard deviation of added Gaussian noise. Default 0.15.

- seed:

  Optional integer seed for reproducibility.

## Value

List with elements `vec` (`NeuroVec`), `mask` (`NeuroVol`), and integer
`truth` labels.

## Examples

``` r
syn <- make_block_synthetic(dims = c(12, 12, 1), ntime = 60, noise = 0.1, seed = 7)
length(unique(syn$truth))
#> [1] 3
```
