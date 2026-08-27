# Explicit cluster-quality estimands

Computes partition, within-cluster temporal, and physical spatial
estimands from final voxel labels. Temporal coherence is the mean
Pearson correlation over all unordered within-cluster voxel pairs.
Spatial dispersion is the root mean squared Euclidean distance from each
voxel coordinate to its final cluster centroid; lower values are more
compact.

## Usage

``` r
cluster_metrics(result, feature_mat = NULL, coords = NULL, truth = NULL)
```

## Arguments

- result:

  A `cluster_result` with one final label per voxel.

- feature_mat:

  Optional finite numeric matrix with voxels in rows and at least two
  time/features in columns. A transposed matrix is accepted when only
  its columns match the number of labels. Every row must be non-constant
  for Pearson correlation.

- coords:

  Optional finite numeric matrix with voxels in rows and spatial
  dimensions in columns. When omitted for a canonical
  `cluster4d_result`, physical millimetre coordinates are reconstructed
  from its mask and affine provenance.

- truth:

  Optional partition labels over exactly the same voxels.

## Value

A named list. Always contains `n_clusters` and `size_summary`. With
`truth`, it adds adjusted Rand, variation of information in bits, and
pairwise Dice. With `feature_mat`, it adds
`temporal_pairwise_correlation`. With coordinates, it adds
`spatial_rms_distance_mm` (or `spatial_rms_distance` when explicit
coordinates are not declared to be millimetres by a canonical result).
