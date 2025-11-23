# Consensus Fusion for SLiCE-MSF

Combines multiple SLiCE-MSF segmentation runs using consensus clustering
to improve stability and reduce variability. This function can help
reduce z-plane artifacts by averaging cluster assignments across
multiple independent runs.

## Usage

``` r
slice_msf_consensus(
  run_results,
  mask,
  nbhd = 8,
  k_fuse = 0.3,
  min_size_fuse = 80,
  use_features = FALSE,
  lambda = 0.7
)
```

## Arguments

- run_results:

  List of results from `slice_msf_single`, each containing labels,
  weights, and sketch matrices.

- mask:

  A `NeuroVol` mask used in the original segmentation.

- nbhd:

  Neighborhood connectivity (4, 6, or 8). Should match the value used in
  the original segmentation. Default is 8.

- k_fuse:

  Scale parameter for consensus fusion graph (0-2). Lower values create
  more refined clusters during fusion. If too low, may fragment
  consensus clusters. Default is 0.30.

- min_size_fuse:

  Minimum cluster size during fusion. Can be set lower than the original
  min_size to preserve small but consistent clusters. Default is 80.

- use_features:

  Include feature similarity in consensus fusion. When TRUE, uses both
  label co-occurrence and temporal similarity from DCT sketches. This
  can improve cross-slice consistency and is recommended when targeting
  exact K. Default is FALSE.

- lambda:

  Mixing parameter (0-1) controlling the balance between label agreement
  and feature similarity when use_features=TRUE. Higher values (0.7-0.9)
  emphasize label agreement, lower values (0.3-0.6) emphasize features.
  Default is 0.7.

## Value

A list with the following components:

- labels:

  Integer vector of consensus cluster labels

- params:

  List of parameters used in consensus fusion

## Details

Consensus fusion works by building a co-association matrix that counts
how often each pair of voxels is assigned to the same cluster across
runs. This matrix is then used as edge weights in a new graph
segmentation. When use_features=TRUE, the edge weights also incorporate
the similarity of DCT sketches, providing additional stability.

For reducing z-artifacts, using use_features=TRUE with lambda around
0.5-0.6 can help ensure that clusters are consistent both in terms of
their boundaries and their temporal characteristics.
