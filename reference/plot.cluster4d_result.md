# Plot cluster4d result

Creates visualization of clustering results. Shows axial, sagittal, and
coronal slices through the clustered volume.

## Usage

``` r
# S3 method for class 'cluster4d_result'
plot(x, slice = NULL, view = "all", colors = NULL, zlim = NULL, ...)
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

- zlim:

  Numeric length-two label range shared by every displayed plane. The
  default, `c(0.5, x$n_clusters + 0.5)`, keeps each positive integer
  label mapped to the same categorical color even when a slice omits
  labels.

- ...:

  Additional arguments passed to plotting functions

## Value

Invisibly returns the plotted data
