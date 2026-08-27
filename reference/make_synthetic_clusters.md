# Create Synthetic Ground Truth Clustering Data

Generates a volume with perfectly separable clusters for testing
clustering algorithm correctness. Each cluster is a spatially contiguous
cubic sub-region with identical timeseries within the cluster and
orthogonal signals between clusters.

## Usage

``` r
make_synthetic_clusters(
  n_time = 100,
  noise_sd = 0,
  n_clusters = 27,
  dim_per_cluster = 3,
  seed = NULL
)
```

## Arguments

- n_time:

  Number of timepoints. Default 100.

- noise_sd:

  Standard deviation of Gaussian noise to add. Default 0 (no noise =
  perfect separability).

- n_clusters:

  Number of clusters. Must be a perfect cube (8, 27, 64, etc.) for cubic
  cluster geometry. Default 27.

- dim_per_cluster:

  Voxels per dimension per cluster. Default 3.

- seed:

  Random seed for reproducibility. Default NULL.

## Value

A list with:

- vec:

  NeuroVec object

- mask:

  NeuroVol mask (all ones)

- true_clusters:

  NeuroVol with ground truth cluster labels

- true_labels:

  Integer vector of true cluster assignments

- n_clusters:

  Number of clusters

- cluster_size:

  Size of each cluster in voxels

- signals:

  Matrix of cluster signals (n_time x n_clusters)

## Details

For the default n_clusters=27, creates a 9x9x9 volume divided into 27
cubic 3x3x3 clusters arranged in a 3x3x3 grid:


    z=1-3:           z=4-6:           z=7-9:
      y                y                y
      ^                ^                ^
      | 7 | 8 | 9 |    | 16| 17| 18|    | 25| 26| 27|
      | 4 | 5 | 6 |    | 13| 14| 15|    | 22| 23| 24|
      | 1 | 2 | 3 |    | 10| 11| 12|    | 19| 20| 21|
      +-----------> x  +-----------> x  +-----------> x

Each cluster has a unique orthogonal signal generated from QR
decomposition of a random matrix, ensuring zero correlation between
clusters.

This creates an "easy" test case where any reasonable clustering
algorithm should achieve perfect or near-perfect recovery with
noise_sd=0.

## Examples

``` r
# Default: 27 clusters, perfect separability
data <- make_synthetic_clusters()
table(data$true_labels)  # 27 voxels per cluster
#> 
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
#> 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 27 
#> 27 
#> 27 

# 8 clusters (2x2x2 grid of 3x3x3 cubes = 6x6x6 volume)
data <- make_synthetic_clusters(n_clusters = 8)

# With noise for more realistic testing
data <- make_synthetic_clusters(noise_sd = 0.5)
```
