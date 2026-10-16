# Detrend Voxel Timeseries

Fast removal of trends from each voxel's timeseries in 4D neuroimaging
data. Supports linear, polynomial, and DCT (discrete cosine transform)
detrending.

## Usage

``` r
detrend_time(
  vec,
  mask = NULL,
  method = c("linear", "poly", "dct"),
  degree = 2,
  n_basis = 4
)
```

## Arguments

- vec:

  A `NeuroVec` object (4D neuroimaging data)

- mask:

  Optional `NeuroVol` mask. If NULL, uses all non-zero voxels.

- method:

  Detrending method: "linear" (default), "poly", or "dct"

- degree:

  For method="poly", the polynomial degree (1=linear, 2=quadratic, etc.)

- n_basis:

  For method="dct", number of DCT basis functions to remove. Higher
  values = more aggressive high-pass filtering.

## Value

A `NeuroVec` with trend removed from each voxel

## Details

Three detrending methods are available:

**Linear (method="linear"):** Fastest option. Removes linear drift by
fitting y = a + b\*t.

**Polynomial (method="poly"):** Removes polynomial trend up to specified
degree. degree=1 is linear, degree=2 adds quadratic, etc. Uses
orthonormalized basis for stability.

**DCT (method="dct"):** Removes low-frequency components using discrete
cosine transform basis. This is a standard high-pass filter for fMRI
data. The n_basis parameter controls the cutoff:

- n_basis=1: mean removal only

- n_basis=2: approximately linear detrend

- n_basis=4-6: typical for fMRI (removes ~0.01 Hz drift for TR=2s)

All methods use fast C++ parallel implementation via RcppParallel.

## Examples

``` r
if (FALSE) { # \dontrun{
vec <- neuroim2::read_vec("functional.nii.gz")

# Linear detrend (fastest)
vec_lin <- detrend_time(vec)

# Quadratic detrend
vec_quad <- detrend_time(vec, method = "poly", degree = 2)

# DCT high-pass filter (common for fMRI)
vec_dct <- detrend_time(vec, method = "dct", n_basis = 5)
} # }
```
