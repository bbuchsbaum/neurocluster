# Normalize Volumes and Detrend Timeseries (Fast C++ Implementation)

Combined operation: first removes volume mean offsets, then removes
linear trend from each voxel's timeseries. More efficient than calling
both functions separately.

## Usage

``` r
normalize_detrend_cpp(data)
```

## Arguments

- data:

  Numeric matrix with voxels as rows and timepoints as columns

## Value

Numeric matrix with volume offsets and linear trends removed

## Details

Performs two operations in a single parallel pass:

1.  Centers each volume to zero mean

2.  Removes linear trend from each voxel's (now centered) timeseries
