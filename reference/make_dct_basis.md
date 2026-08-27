# Generate DCT-II Basis Matrix

Creates a discrete cosine transform (type II) basis matrix for
detrending.

## Usage

``` r
make_dct_basis(n_time, n_basis)
```

## Arguments

- n_time:

  Number of timepoints

- n_basis:

  Number of DCT basis functions (including constant)

## Value

Numeric matrix (n_time x n_basis) with orthonormal DCT basis

## Details

DCT basis functions are excellent for removing low-frequency drift in
fMRI. The first basis is constant (mean), subsequent bases capture
increasingly higher frequency components.
