# Find boundary voxels (voxels with neighbors having different labels)

Find boundary voxels (voxels with neighbors having different labels)

## Usage

``` r
find_boundary_voxels_cpp(voxel_labels, neighbor_indices)
```

## Arguments

- voxel_labels:

  Integer vector of cluster assignments

- neighbor_indices:

  Integer matrix (voxels x K) of nearest neighbor indices (1-based)

## Value

Integer vector of boundary voxel indices (1-based for R)
