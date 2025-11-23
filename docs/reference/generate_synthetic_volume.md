# Generate a synthetic 4D neuroimaging volume with labeled clusters

Convenient synthetic data for examples, vignettes, and tests. Each
scenario produces a `NeuroVec` (time series volume), a `NeuroVol` mask,
and the ground-truth cluster labels so you can benchmark algorithms or
build illustrations quickly.

## Usage

``` r
generate_synthetic_volume(
  scenario = c("gaussian_blobs", "z_layers", "checkerboard"),
  dims = c(16, 16, 8),
  n_clusters = 4,
  n_time = 30,
  spread = c(3.5, 3.5, 2),
  noise_sd = 0.05,
  amplitude_range = c(0.8, 1.2),
  spacing_mm = c(3, 3, 3),
  seed = NULL
)
```

## Arguments

- scenario:

  Pattern to embed. One of `"gaussian_blobs"`, `"z_layers"`, or
  `"checkerboard"`.

- dims:

  Integer vector of length 3 giving the spatial grid dimensions.

- n_clusters:

  Number of latent clusters to embed.

- n_time:

  Number of time points.

- spread:

  Characteristic spread of each Gaussian blob (scalar or length 3).

- noise_sd:

  Standard deviation of additive Gaussian noise.

- amplitude_range:

  Range used when scaling each voxel's amplitude.

- spacing_mm:

  Physical voxel size supplied to the generated `NeuroSpace`.

- seed:

  Optional seed for reproducibility.

## Value

A list with elements `vec`, `mask`, `truth`, `coords`, `patterns`,
`weights`, `dims`, `n_clusters`, and `scenario`.

## Examples

``` r
syn <- generate_synthetic_volume(
  scenario = "gaussian_blobs",
  dims = c(12, 12, 6),
  n_clusters = 4,
  seed = 1
)
str(syn)
#> List of 9
#>  $ vec       :Formal class 'DenseNeuroVec' [package "neuroim2"] with 3 slots
#>   .. ..@ .Data: num [1:12, 1:12, 1:6, 1:30] -0.35 -0.311 -0.344 -0.515 -0.328 ...
#>   .. ..@ label: chr "none"
#>   .. ..@ space:Formal class 'NeuroSpace' [package "neuroim2"] with 6 slots
#>   .. .. .. ..@ dim    : int [1:4] 12 12 6 30
#>   .. .. .. ..@ origin : num [1:3] 0 0 0
#>   .. .. .. ..@ spacing: num [1:3] 3 3 3
#>   .. .. .. ..@ axes   :Formal class 'AxisSet3D' [package "neuroim2"] with 4 slots
#>   .. .. .. .. .. ..@ k   :Formal class 'NamedAxis' [package "neuroim2"] with 2 slots
#>   .. .. .. .. .. .. .. ..@ axis     : chr "Inferior-to-Superior"
#>   .. .. .. .. .. .. .. ..@ direction: num [1:3] 0 0 1
#>   .. .. .. .. .. ..@ j   :Formal class 'NamedAxis' [package "neuroim2"] with 2 slots
#>   .. .. .. .. .. .. .. ..@ axis     : chr "Posterior-to-Anterior"
#>   .. .. .. .. .. .. .. ..@ direction: num [1:3] 0 1 0
#>   .. .. .. .. .. ..@ i   :Formal class 'NamedAxis' [package "neuroim2"] with 2 slots
#>   .. .. .. .. .. .. .. ..@ axis     : chr "Left-to-Right"
#>   .. .. .. .. .. .. .. ..@ direction: num [1:3] 1 0 0
#>   .. .. .. .. .. ..@ ndim: int 3
#>   .. .. .. ..@ trans  : num [1:4, 1:4] 3 0 0 0 0 3 0 0 0 0 ...
#>   .. .. .. ..@ inverse: num [1:4, 1:4] 0.333 0 0 0 0 ...
#>   .. ..$ dim: int [1:4] 12 12 6 30
#>  $ mask      :Formal class 'DenseNeuroVol' [package "neuroim2"] with 2 slots
#>   .. ..@ .Data: logi [1:12, 1:12, 1:6] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>   .. ..@ space:Formal class 'NeuroSpace' [package "neuroim2"] with 6 slots
#>   .. .. .. ..@ dim    : int [1:3] 12 12 6
#>   .. .. .. ..@ origin : num [1:3] 0 0 0
#>   .. .. .. ..@ spacing: num [1:3] 3 3 3
#>   .. .. .. ..@ axes   :Formal class 'AxisSet3D' [package "neuroim2"] with 4 slots
#>   .. .. .. .. .. ..@ k   :Formal class 'NamedAxis' [package "neuroim2"] with 2 slots
#>   .. .. .. .. .. .. .. ..@ axis     : chr "Inferior-to-Superior"
#>   .. .. .. .. .. .. .. ..@ direction: num [1:3] 0 0 1
#>   .. .. .. .. .. ..@ j   :Formal class 'NamedAxis' [package "neuroim2"] with 2 slots
#>   .. .. .. .. .. .. .. ..@ axis     : chr "Posterior-to-Anterior"
#>   .. .. .. .. .. .. .. ..@ direction: num [1:3] 0 1 0
#>   .. .. .. .. .. ..@ i   :Formal class 'NamedAxis' [package "neuroim2"] with 2 slots
#>   .. .. .. .. .. .. .. ..@ axis     : chr "Left-to-Right"
#>   .. .. .. .. .. .. .. ..@ direction: num [1:3] 1 0 0
#>   .. .. .. .. .. ..@ ndim: int 3
#>   .. .. .. ..@ trans  : num [1:4, 1:4] 3 0 0 0 0 3 0 0 0 0 ...
#>   .. .. .. ..@ inverse: num [1:4, 1:4] 0.333 0 0 0 0 ...
#>   .. ..$ dim: int [1:3] 12 12 6
#>  $ truth     : int [1:864] 4 4 4 4 4 4 4 4 4 4 ...
#>  $ coords    : int [1:864, 1:3] 1 2 3 4 5 6 7 8 9 10 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "x" "y" "z"
#>  $ patterns  : num [1:4, 1:30] 0.667 0.167 0.972 -0.99 0.435 ...
#>  $ weights   : num [1:864] 0.248 0.312 0.362 0.386 0.38 ...
#>  $ dims      : num [1:3] 12 12 6
#>  $ n_clusters: num 4
#>  $ scenario  : chr "gaussian_blobs"
truth_vol <- array(syn$truth, syn$dims)
clusvol <- neuroim2::ClusteredNeuroVol(syn$mask > 0, clusters = syn$truth)
plot(clusvol, slice = c(6, 6, 3), view = "axial")

#> Warning: "slice" is not a graphical parameter
#> Warning: "view" is not a graphical parameter
#> Error in plot.xy(xy, type, ...): 'x' and 'y' lengths differ in plot.xy()
```
