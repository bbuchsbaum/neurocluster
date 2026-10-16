# Generate Polynomial Basis Matrix

Creates a polynomial basis matrix for detrending.

## Usage

``` r
make_poly_basis(n_time, degree)
```

## Arguments

- n_time:

  Number of timepoints

- degree:

  Polynomial degree (0 = constant, 1 = linear, 2 = quadratic, etc.)

## Value

Numeric matrix (n_time x (degree+1)) with orthonormalized polynomial
basis

## Details

Uses Gram-Schmidt orthonormalization for numerical stability. Degree 0 =
mean removal, 1 = linear detrend, 2 = quadratic, etc.
