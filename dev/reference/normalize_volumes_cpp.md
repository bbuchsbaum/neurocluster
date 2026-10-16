# Normalize Volumes by Removing Mean Offset (Fast C++ Implementation)

Removes the mean offset from each volume/timepoint in a matrix. This
centers each volume to have zero mean across voxels.

## Usage

``` r
normalize_volumes_cpp(data)
```

## Arguments

- data:

  Numeric matrix with voxels as rows and timepoints as columns

## Value

Numeric matrix with each column (volume) centered to mean zero

## Details

For each timepoint `t`, computes `mean_t = mean(data[, t])` and
subtracts it from all voxels: `output[, t] = data[, t] - mean_t`

Uses RcppParallel for fast parallel computation across voxels.
