# Supervoxel Clustering on a Surface

Cluster feature data on a cortical surface or mesh using a
supervoxel-like approach.

## Usage

``` r
supervoxel_cluster_surface(
  bsurf,
  K = 500,
  sigma1 = 1,
  sigma2 = 5,
  iterations = 50,
  connectivity = 6,
  use_medoid = FALSE
)
```

## Arguments

- bsurf:

  A `NeuroSurface` or similar object with geometry, coords, and data.

- K:

  Number of clusters.

- sigma1:

  Heat kernel bandwidth for feature similarity (data vectors).

- sigma2:

  Heat kernel bandwidth for spatial similarity (coordinate vectors).

- iterations:

  Max iterations.

- connectivity:

  Neighborhood size on the surface (e.g., \# of nearest mesh neighbors).

- use_medoid:

  Whether to use medoids for cluster centers.

## Value

A `list` with:

- clusvol:

  A `NeuroSurface` storing the final clustering result.

- clusters:

  Integer vector of cluster assignments (one per vertex).

- centers:

  Matrix of cluster centers.

- coord_centers:

  Matrix of spatial centroid coordinates.

- index_sets:

  List of vertex indices for each cluster.
