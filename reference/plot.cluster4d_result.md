# Plot cluster4d result

Creates visualization of clustering results. Shows axial, sagittal, and
coronal slices through the clustered volume.

## Usage

``` r
# S3 method for class 'cluster4d_result'
plot(x, slice = NULL, view = "all", colors = NULL, ...)
```

## Arguments

- x:

  A cluster4d_result object

- slice:

  Slice specification. Can be:

  - NULL (default): Shows middle slices

  - Numeric vector of length 3: c(x, y, z) coordinates

  - "montage": Shows multiple slices

- view:

  Viewing plane: "axial", "sagittal", "coronal", or "all"

- colors:

  Color palette for clusters. Default uses rainbow colors.

- ...:

  Additional arguments passed to plotting functions

## Value

Invisibly returns the plotted data
