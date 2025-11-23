# Adaptive Correlation Superclustering (ACSC)

Clusters fMRI voxels into spatially-coherent groups based on temporal
correlation and spatial proximity. Includes optional refinement for
boundary corrections. The algorithm supports parallel processing via the
future package for improved performance.

## Usage

``` r
acsc(
  bvec,
  mask,
  block_size = 2,
  ann_k = 10,
  alpha = 0.5,
  correlation_metric = c("pearson", "spearman", "robust"),
  spatial_weighting = c("gaussian", "binary"),
  refine = TRUE,
  max_refine_iter = 5,
  K = NULL
)
```

## Arguments

- bvec:

  A NeuroVec-like object containing 4D fMRI data.

- mask:

  A NeuroVol-like object (logical or numeric mask).

- block_size:

  Approximate side length of blocks (e.g., 2 or 3). Must be \> 0.

- ann_k:

  Number of approximate (or exact) nearest neighbors per block. Must be
  \>= 1.

- alpha:

  Weighting for correlation vs. spatial proximity (0 \<= alpha \<= 1).

- correlation_metric:

  Correlation metric ("pearson", "spearman", "robust").

- spatial_weighting:

  Spatial adjacency weighting ("gaussian", "binary").

- refine:

  Logical; whether to refine boundaries.

- max_refine_iter:

  Maximum iterations for boundary refinement. Must be \>= 0.

- K:

  (Optional) Desired number of clusters.

## Value

A list with elements:

- cluster_map:

  3D array with cluster labels per voxel.

- graph:

  An `igraph` object used for clustering.

- init_block_label:

  Initial coarse partition (3D array) matching `mask` dimensions.

## Details

### C++ Acceleration

ACSC now includes **C++ acceleration for boundary refinement** using
RcppParallel, providing 3-6x speedup for typical datasets. The C++
implementation:

- Uses optimized correlation via normalized dot products (10-15x faster
  than R's [`cor()`](https://rdrr.io/r/stats/cor.html))

- Processes boundary voxels in parallel using multiple CPU cores

- Automatically falls back to R implementation if C++ fails

#### Performance Gains:

- **Small datasets** (\<1,000 voxels): 1.5-2x overall speedup

- **Medium datasets** (1,000-5,000 voxels): 2-4x overall speedup

- **Large datasets** (\>5,000 voxels): 3-6x overall speedup

- **Boundary refinement phase**: 6-8x faster than pure R implementation

#### C++ Implementation Details:

The C++ acceleration normalizes feature vectors to unit length, enabling
fast correlation computation via dot products. This is mathematically
equivalent to Pearson correlation for centered data and provides
significant performance benefits.

### Parallelization Strategy

ACSC uses **dual-layer parallelization**:

1.  **R-level parallelization** (via future package):

    - Data preprocessing (detrending)

    - Block summary computation

    - Graph edge construction

    - User-configurable across platforms

2.  **C++ thread parallelization** (via RcppParallel):

    - Boundary voxel refinement

    - Automatic multi-core utilization

    - No configuration needed

#### Configuring R-level Parallelization:

    library(future)

    # Sequential (default)
    plan(sequential)
    result <- acsc(bvec, mask, K = 100)

    # Parallel on local machine (uses all cores)
    plan(multisession)
    result <- acsc(bvec, mask, K = 100)

    # Parallel with specific number of workers
    plan(multisession, workers = 4)
    result <- acsc(bvec, mask, K = 100)

    # On a cluster
    plan(cluster, workers = c("node1", "node2", "node3"))
    result <- acsc(bvec, mask, K = 100)

    # Reset to sequential
    plan(sequential)

#### Performance Characteristics:

- **Best speedup**: Boundary refinement (6-8x) and graph construction
  (2-3x)

- **Overhead**: Small for C++ calls, moderate for future parallelization

- **Memory**: Each future worker needs data copy; C++ threads share
  memory

- **Optimal workers**: Usually matches physical cores (not threads)

- **Break-even point**: C++ benefits all dataset sizes; future benefits
  \>5,000 voxels

#### Performance Tips:

- C++ acceleration is enabled by default and recommended for all use
  cases

- For small datasets (\<1,000 voxels), sequential future plan may be
  faster

- Use `plan(multisession)` on Windows/macOS for stability

- Use `plan(multicore)` on Linux for lower memory overhead

- C++ and future parallelization work together without conflicts

- Monitor memory usage with many workers on large datasets
