# DCT Detrending (Convenience Function)

Removes low-frequency components using discrete cosine transform.

## Usage

``` r
detrend_dct_cpp(data, n_basis = 4L)
```

## Arguments

- data:

  Numeric matrix with voxels as rows and timepoints as columns

- n_basis:

  Number of DCT basis functions to remove (including constant). Default
  4 removes constant, linear, and first two cosine harmonics.

## Value

Numeric matrix with low-frequency components removed

## Details

DCT detrending is commonly used in fMRI preprocessing. The number of
basis functions controls the high-pass filter cutoff:

- n_basis = 1: mean removal only

- n_basis = 2: ~ linear detrend

- n_basis = 4-6: typical for fMRI (removes drift \< ~0.01 Hz for TR=2s)

Higher n_basis = more aggressive high-pass filtering.
