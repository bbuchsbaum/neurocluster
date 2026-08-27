# Normalize and Detrend Neuroimaging Data

Combined fast preprocessing: removes volume mean offsets and linear
trends in a single efficient operation.

## Usage

``` r
normalize_detrend(vec, mask = NULL)
```

## Arguments

- vec:

  A `NeuroVec` object (4D neuroimaging data)

- mask:

  Optional `NeuroVol` mask. If NULL, uses all non-zero voxels.

## Value

A `NeuroVec` with volume offsets and linear trends removed

## Details

Combines
[`normalize_volumes`](https://bbuchsbaum.github.io/neurocluster/reference/normalize_volumes.md)
and
[`detrend_time`](https://bbuchsbaum.github.io/neurocluster/reference/detrend_time.md)
in a single pass for efficiency. The operations are applied in order:

1.  Remove volume-wise mean offset

2.  Remove linear trend from each voxel's timeseries

Uses fast C++ parallel implementation via RcppParallel.

## Examples

``` r
if (FALSE) { # \dontrun{
vec <- neuroim2::read_vec("functional.nii.gz")
vec_clean <- normalize_detrend(vec)
} # }
```
