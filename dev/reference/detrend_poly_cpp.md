# Polynomial Detrending (Convenience Function)

Removes polynomial trend up to specified degree from each voxel.

## Usage

``` r
detrend_poly_cpp(data, degree = 1L)
```

## Arguments

- data:

  Numeric matrix with voxels as rows and timepoints as columns

- degree:

  Polynomial degree (0 = mean, 1 = linear, 2 = quadratic, etc.)

## Value

Numeric matrix with polynomial trend removed

## Details

Convenience wrapper that generates polynomial basis and calls
detrend_basis_cpp. For degree=1, this is equivalent to detrend_time_cpp
but slightly slower.
