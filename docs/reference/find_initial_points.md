# Find Initial Cluster Centers for Supervoxel Algorithm (Deprecated)

**Deprecated:** This function is now available as
[`find_gradient_seeds()`](find_gradient_seeds.md) in `cluster4d_init.R`.
This version is maintained for backward compatibility.

This function finds the initial cluster centers for a supervoxel
algorithm. Supervoxels are used to partition 3D image data into
volumetric regions, grouping similar voxels together. The initial
cluster centers are crucial for the performance and quality of the final
supervoxels.

## Usage

``` r
find_initial_points(cds, grad, K = 100)
```

## Arguments

- cds:

  A matrix or data frame representing the spatial coordinates of the
  voxels.

- grad:

  A vector representing the gradient values of the voxels.

- K:

  The desired number of supervoxels (clusters) in the output (default:
  100).

## Value

A list containing two elements: selected - a vector of the selected
indices corresponding to the initial cluster centers, coords - a matrix
or data frame with the spatial coordinates of the initial cluster
centers.

## See also

[`find_gradient_seeds`](find_gradient_seeds.md),
[`initialize_clusters`](cluster4d_init.md)
