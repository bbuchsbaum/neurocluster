# neurocluster

R package for spatially constrained clustering of neuroimaging data.

## Overview

The `neurocluster` package provides multiple algorithms for clustering
4D neuroimaging data (3D space + time) while respecting spatial
structure. All algorithms balance between feature similarity (e.g., time
series correlation) and spatial proximity to produce coherent,
interpretable clusters.

## Installation

``` r
# Install from local source
devtools::install()

# Or build and install
R CMD build .
R CMD INSTALL neurocluster_*.tar.gz
```

## Quick Start

The package provides a unified [`cluster4d()`](reference/cluster4d.md)
interface for all clustering methods:

``` r
library(neurocluster)
library(neuroim2)

# Load your data
mask <- read_vol("brain_mask.nii")
vec <- read_vec("fmri_data.nii.gz")

# Cluster with default method (supervoxels)
result <- cluster4d(vec, mask, n_clusters = 100)

# Try different methods
result_snic <- cluster4d(vec, mask, n_clusters = 100, method = "snic")
result_msf <- cluster4d(vec, mask, n_clusters = 100, method = "slice_msf")
```

## Available Methods

### 1. Supervoxels (`method = "supervoxels"`)

- **Best for**: General purpose, high-quality clusters
- **Characteristics**: Iterative refinement using heat kernels
- **Key parameters**: `sigma1` (feature bandwidth), `sigma2` (spatial
  bandwidth)
- **Parallelized**: Yes

### 2. SNIC - Simple Non-Iterative Clustering (`method = "snic"`)

- **Best for**: Fast results, moderate-sized data
- **Characteristics**: Priority queue-based, single pass
- **Key parameters**: `compactness`
- **Parallelized**: No (inherently sequential)

### 3. SLIC - Simple Linear Iterative Clustering (`method = "slic"`)

- **Best for**: Consistent cluster sizes, preservation of K
- **Characteristics**: Iterative with local search windows
- **Key parameters**: `compactness`, `preserve_k`
- **Parallelized**: Yes
- **Note**: Requires slic4d_core C++ implementation

### 4. Slice-MSF (`method = "slice_msf"`)

- **Best for**: Large datasets, parallel processing
- **Characteristics**: Slice-wise processing with consensus
- **Key parameters**: `num_runs`, `consensus`, `stitch_z`
- **Parallelized**: Yes (across slices)

### 5. FLASH-3D (`method = "flash3d"`)

- **Best for**: Fast processing, temporal coherence
- **Characteristics**: DCT-based compression
- **Key parameters**: `lambda_t`, `bits`, `dctM`
- **Parallelized**: Yes

## Common Parameters

All methods accept these standardized parameters:

- `n_clusters`: Target number of clusters
- `spatial_weight`: Balance between spatial and feature similarity (0-1)
- `max_iterations`: Maximum iterations for convergence
- `connectivity`: Neighborhood structure (6, 26, or 27)
- `parallel`: Enable parallel processing (where supported)
- `verbose`: Print progress information

## Examples

### Basic Usage

``` r
# Load test data
mask <- NeuroVol(array(1, c(64,64,32)), NeuroSpace(c(64,64,32)))
vec <- replicate(100, NeuroVol(array(rnorm(64*64*32), c(64,64,32)),
                 NeuroSpace(c(64,64,32))), simplify=FALSE)
vec <- do.call(concat, vec)

# Cluster with balanced spatial/feature weighting
result <- cluster4d(vec, mask, 
                   n_clusters = 50,
                   spatial_weight = 0.5)

# Access results
clusters <- result$cluster          # Cluster assignments
centers <- result$centers            # Feature centers
spatial_centers <- result$coord_centers  # Spatial centroids
```

### Choosing Parameters

``` r
# Get intelligent parameter suggestions
n_voxels <- sum(mask > 0)
n_timepoints <- dim(vec)[4]

# Get balanced recommendations
params <- suggest_cluster4d_params(n_voxels, n_timepoints, 
                                   priority = "balanced")

# Use recommended settings
result <- cluster4d(vec, mask, 
                   n_clusters = params$n_clusters,
                   method = params$recommended_method)
```

### Method-Specific Features

``` r
# Supervoxels with custom kernel bandwidths
result <- cluster4d(vec, mask, n_clusters = 100,
                   method = "supervoxels",
                   sigma1 = 2.0,    # Feature kernel
                   sigma2 = 3.0,    # Spatial kernel
                   use_gradient = TRUE)

# Slice-MSF with consensus clustering
result <- cluster4d(vec, mask, n_clusters = 100,
                   method = "slice_msf",
                   num_runs = 5,      # Multiple runs
                   consensus = TRUE,  # Consensus clustering
                   stitch_z = TRUE)   # Reduce z-plane artifacts

# FLASH-3D with temporal smoothing
result <- cluster4d(vec, mask, n_clusters = 100,
                   method = "flash3d",
                   lambda_t = 1.5,    # Temporal regularization
                   bits = 128)        # Quantization precision
```

### Comparing Methods

``` r
# Run multiple methods
methods <- c("supervoxels", "snic", "slice_msf")
results <- lapply(methods, function(m) {
  cluster4d(vec, mask, n_clusters = 50, method = m)
})
names(results) <- methods

# Compare results
comparison <- compare_cluster4d(results)
print(comparison)
```

## Performance Considerations

### Small Datasets (\< 10,000 voxels)

- Use `supervoxels` for quality
- Use `snic` for speed
- Parallel processing overhead may not be worth it

### Medium Datasets (10,000 - 100,000 voxels)

- `supervoxels` with `parallel = TRUE`
- `slice_msf` with `num_runs = 1` for speed
- `flash3d` for balanced speed/quality

### Large Datasets (\> 100,000 voxels)

- `slice_msf` with parallel processing
- `flash3d` for fastest results
- Consider reducing `n_clusters` for tractability

### Memory Constraints

- Use `snic` (non-iterative, low memory)
- Use `slice_msf` (processes slices independently)
- Reduce time series length if possible

## Advanced Usage

### Custom Initialization

``` r
# Use gradient-based initialization for better seeds
result <- cluster4d(vec, mask, n_clusters = 100,
                   method = "supervoxels",
                   use_gradient = TRUE)
```

### Working with Results

``` r
# Extract ClusteredNeuroVol
clusvol <- result$clusvol

# Get cluster sizes
cluster_sizes <- table(result$cluster)

# Find cluster for specific voxel
voxel_coord <- c(32, 32, 16)
voxel_idx <- coord_to_index(mask, voxel_coord)
voxel_cluster <- result$cluster[voxel_idx]

# Extract time series for a cluster
cluster_id <- 5
cluster_mask <- result$cluster == cluster_id
cluster_timeseries <- series(vec, which(cluster_mask))
```

### Visualization

``` r
# Plot cluster centers
plot(result)

# Visualize spatial distribution
library(neuroim2)
writeVol(result$clusvol, "clusters.nii.gz")

# View specific cluster
cluster_vol <- NeuroVol(array(0, dim(mask)), space(mask))
cluster_vol[result$cluster == 5] <- 1
writeVol(cluster_vol, "cluster_5.nii.gz")
```

## Backward Compatibility

The original function interfaces are still available:

``` r
# Original supervoxels syntax
result <- supervoxels(vec, mask, K = 100, alpha = 0.3)

# Original SNIC syntax  
result <- snic(vec, mask, K = 100, compactness = 5)

# Original slice_msf syntax
result <- slice_msf(vec, mask, target_k_global = 100)
```

## Citation

If you use this package in your research, please cite:

    @software{neurocluster,
      title = {neurocluster: Spatially Constrained Clustering for Neuroimaging},
      author = {Your Name},
      year = {2024},
      url = {https://github.com/yourusername/neurocluster}
    }

## Contributing

Contributions are welcome! Please submit issues and pull requests on
GitHub.

## License

This package is licensed under [LICENSE](#license).

## See Also

- [neuroim2](https://github.com/bbuchsbaum/neuroim2): Neuroimaging data
  structures
- [neurosurf](https://github.com/bbuchsbaum/neurosurf): Surface-based
  clustering
- [neighborweights](https://github.com/bbuchsbaum/neighborweights):
  Spatial weight matrices

## Albers theme

This package uses the albersdown theme. Vignettes are styled with
`vignettes/albers.css` and a local `vignettes/albers.js`; the palette
family is provided via `params$family` (default ‘red’). The pkgdown site
uses `template: { package: albersdown }`.

## Albers theme

This package uses the albersdown theme. Vignettes are styled with
`vignettes/albers.css` and a local `vignettes/albers.js`; the palette
family is provided via `params$family` (default ‘red’). The pkgdown site
uses `template: { package: albersdown }`.

## Albers theme

This package uses the albersdown theme. Vignettes are styled with
`vignettes/albers.css` and a local `vignettes/albers.js`; the palette
family is provided via `params$family` (default ‘red’). The pkgdown site
uses `template: { package: albersdown }`.
