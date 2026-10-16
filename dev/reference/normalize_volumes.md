# Normalize Volumes by Removing Mean Offset

Fast removal of volume-wise mean offset from 4D neuroimaging data. Each
volume/timepoint is centered to have zero mean across voxels.

## Usage

``` r
normalize_volumes(vec, mask = NULL)
```

## Arguments

- vec:

  A `NeuroVec` object (4D neuroimaging data)

- mask:

  Optional `NeuroVol` mask. If NULL, uses all non-zero voxels.

## Value

A `NeuroVec` with each volume centered to zero mean

## Details

This function removes global signal offset that can vary across
timepoints due to scanner drift, physiological noise, or other factors.
For each timepoint t:

\$\$output\[,,,t\] = input\[,,,t\] - mean(input\[,,,t\])\$\$

Uses fast C++ parallel implementation via RcppParallel.

## Examples

``` r
if (FALSE) { # \dontrun{
vec <- neuroim2::read_vec("functional.nii.gz")
vec_normalized <- normalize_volumes(vec)
} # }
```
