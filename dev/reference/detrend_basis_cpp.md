# Flexible Detrending with Basis Functions (Fast C++ Implementation)

Removes trends from each voxel's timeseries by projecting out a basis.
Supports polynomial or DCT (discrete cosine transform) basis functions.

## Usage

``` r
detrend_basis_cpp(data, basis)
```

## Arguments

- data:

  Numeric matrix with voxels as rows and timepoints as columns

- basis:

  Orthonormal basis matrix (n_time x n_basis) to project out. Use
  make_dct_basis() or make_poly_basis() to generate.

## Value

Numeric matrix with basis components removed from each row

## Details

For orthonormal basis B, computes residuals as: y_detrend = y - B \*
(B^T \* y)

This is equivalent to regressing out the basis and keeping residuals.

## Examples

``` r
if (FALSE) { # \dontrun{
# Remove linear + quadratic trend (polynomial degree 2)
basis <- make_poly_basis(n_time = 100, degree = 2)
data_detrend <- detrend_basis_cpp(data, basis)

# Remove low-frequency drift with DCT (first 5 components)
basis <- make_dct_basis(n_time = 100, n_basis = 5)
data_detrend <- detrend_basis_cpp(data, basis)
} # }
```
