# Detrend Voxel Timeseries (Fast C++ Implementation)

Removes linear trend from each voxel's timeseries using fast parallel
linear regression.

## Usage

``` r
detrend_time_cpp(data)
```

## Arguments

- data:

  Numeric matrix with voxels as rows and timepoints as columns

## Value

Numeric matrix with linear trend removed from each row (voxel)

## Details

For each voxel, fits a linear model y = intercept + slope \* t and
subtracts the fitted values: residuals = y - fitted.

Uses RcppParallel for fast parallel computation across voxels. The
linear regression is computed efficiently using precomputed time
statistics.
