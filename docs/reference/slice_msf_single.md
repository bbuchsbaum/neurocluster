# Single Run SLiCE-MSF Segmentation

Lower-level function that performs a single run of SLiCE-MSF
segmentation without consensus fusion. This function is useful for
testing parameters or when you need direct access to the DCT sketch and
reliability weights. Most users should use `slice_msf` instead for
better stability through consensus.

## Usage

``` r
slice_msf_single(
  vec,
  mask,
  r = 12,
  k = 0.32,
  min_size = 80,
  nbhd = 8,
  stitch_z = TRUE,
  theta_link = 0.85,
  min_contact = 1,
  gamma = 1.5,
  z_mult = 0
)
```

## Arguments

- vec:

  A `NeuroVec` instance supplying the time series data.

- mask:

  A `NeuroVol` mask defining the voxels to include.

- r:

  DCT sketch rank (number of basis functions). Higher values preserve
  more temporal detail but increase computation. Range: 4-20. Default is
  12.

- k:

  Scale parameter for Felzenszwalb-Huttenlocher segmentation (0-2).
  Controls the trade-off between under- and over-segmentation. Smaller
  values (0.1-0.3) create more clusters, larger values (0.5-1.0) create
  fewer, larger clusters. Default is 0.32.

- min_size:

  Minimum cluster size in voxels. Clusters smaller than this are merged
  with neighbors. Default is 80.

- nbhd:

  Neighborhood connectivity for within-slice edges (4 or 8). Note that 6
  is automatically mapped to 8. Default is 8.

- stitch_z:

  Enable cross-slice stitching to create 3D clusters. When FALSE, each
  slice is treated independently. Default is TRUE.

- theta_link:

  Centroid correlation threshold for z-stitching (0-1). Lower values
  allow more aggressive merging across slices. Only used if
  stitch_z=TRUE. Default is 0.85.

- min_contact:

  Minimum number of vertically adjacent voxels required to consider
  stitching two clusters. Only used if stitch_z=TRUE. Default is 1.

- gamma:

  Reliability weighting exponent based on split-half correlations.
  Higher values give more weight to reliable voxels. Default is 1.5.

- z_mult:

  Z-smoothing factor (0-1). Values \> 0 blend slice sketches prior to
  clustering to reduce boundary artifacts. Default is 0.0.

## Value

A list with the following components:

- labels:

  Integer vector of cluster labels for each voxel

- weights:

  Numeric vector of reliability weights for each voxel

- sketch:

  Matrix of DCT coefficients (r × n_voxels)

- params:

  List of parameters used in the segmentation

## Details

This function exposes the core SLiCE-MSF algorithm without consensus
fusion. The returned sketch matrix contains the DCT coefficients that
summarize each voxel's time series, and the weights indicate the
reliability of each voxel based on split-half correlations. These can be
useful for diagnostic purposes or custom post-processing.

The k parameter is particularly important: it directly controls the
scale of segmentation. For typical fMRI data, values between 0.2-0.5
work well. Lower values are needed for fine-grained parcellation, while
higher values produce coarser segmentation.
