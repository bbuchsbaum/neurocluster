# API overview

## Interfaces

- Unified entry: [`cluster4d()`](../reference/cluster4d.md)
- Wrappers: `cluster4d_supervoxels`, `cluster4d_snic`, `cluster4d_slic`,
  `cluster4d_slice_msf`, `cluster4d_flash3d`

## Algorithms

- Legacy cores: `supervoxels`, `snic`, `slic4d_supervoxels`,
  `slice_msf`, `supervoxels_flash3d`

## Utilities

- Validation & comparison: `validate_cluster4d`, `compare_cluster4d`
- Helpers: `compute_centroids`, `knn_shrink`, `spatial_gradient`,
  `tesselate`

## Variants

- Surfaces: `supervoxel_cluster_surface`
- Time: `supervoxel_cluster_time` includes: in_header: \|-
