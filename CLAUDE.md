# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Commands

### R Package Development
- **Build package**: `R CMD build .` or use RStudio build tools
- **Install package**: `R CMD INSTALL .` or `devtools::install()`
- **Check package**: `R CMD check .` or `devtools::check()`
- **Run tests**: `devtools::test()` or `testthat::test_dir("tests/testthat")`
- **Generate documentation**: `devtools::document()` (uses roxygen2)

### Testing
- Tests are in `tests/testthat/` directory
- Use `testthat::test_file("tests/testthat/test_turboclust.R")` to run single test
- Test data is available in `testdata/` directory with example neuroimaging files

### C++ Compilation
- Package uses Rcpp for C++ integration
- Source files in `src/` are automatically compiled during package build
- Use `devtools::clean_dll()` if experiencing C++ compilation issues

## Architecture Overview

### Core Functionality
This is an R package for spatially constrained clustering of neuroimaging data. The package provides multiple clustering algorithms designed to respect spatial structure in 3D brain volumes.

### Unified Interface: cluster4d Framework
The package provides a unified `cluster4d()` interface (`R/cluster4d.R`) that standardizes access to all clustering methods:
- Single entry point with `method` parameter selection
- Standardized parameters across all methods
- Consistent result structure (`cluster4d_result` class)
- Backward compatibility with original function interfaces

### Key Clustering Algorithms

1. **Supervoxels** (`method = "supervoxels"` in cluster4d, or `supervoxels()` in `R/supervoxels.R`)
   - Primary algorithm for spatially constrained clustering
   - Iterative algorithm balancing feature similarity and spatial proximity
   - Uses heat kernels with bandwidths `sigma1` (features) and `sigma2` (coordinates)
   - Parallelized for performance

2. **SNIC** (`method = "snic"` in cluster4d, or `snic()` in `R/snic.R`) 
   - Simple Non-Iterative Clustering algorithm
   - Uses priority queue approach with compactness parameter
   - Sequential processing (not parallelized due to algorithm design)

3. **SLIC** (`method = "slic"` in cluster4d, or via `R/slic4d.R`)
   - Simple Linear Iterative Clustering with 4D support
   - Preserves cluster count with local search windows
   - Parallelized implementation

4. **Slice-MSF** (`method = "slice_msf"` in cluster4d, or `slice_msf()` in `R/slice_cluster.R`)
   - Processes 2D slices independently then merges
   - Consensus clustering option for stability
   - Parallelized across slices

5. **FLASH-3D** (`method = "flash3d"` in cluster4d, or `supervoxels_flash3d()` in `R/flash3d.R`)
   - Fast clustering using DCT compression
   - Temporal coherence preservation
   - Parallelized implementation

6. **Meta-clustering** (`meta_clust()` in `R/meta_clust.R`)
   - Hierarchical clustering of cluster results
   - Consensus clustering via `merge_clus()` functions

### Data Types and Dependencies
- Built on `neuroim2` package for neuroimaging data structures
- Primary data types: `NeuroVec` (4D data), `NeuroVol` (3D volumes), `ClusteredNeuroVol` (clustering results)
- Uses `neighborweights` package for spatial adjacency and Laplacian calculations
- Surface clustering support via `neurosurf` package

### C++ Implementation
- Performance-critical functions in `src/` directory:
  - `find_best.cpp`: Heat kernel computations and cluster assignment optimization
  - `snic.cpp`: SNIC algorithm implementation  
  - `correlation_gradient.cpp`: Gradient calculations
- All C++ functions exported via `RcppExports.cpp`

### Cluster Result Structure
All clustering functions return standardized `cluster4d_result` objects (inheriting from `cluster_result`) with:
- `clusvol`: `ClusteredNeuroVol` with spatial cluster assignments
- `cluster`: Integer vector of cluster IDs
- `centers`: Feature space centroids
- `coord_centers`: Spatial coordinate centroids
- `n_clusters`: Actual number of clusters found
- `method`: Method used for clustering
- `parameters`: List of parameters used

### Experimental Features
`experimental/` directory contains developmental algorithms:
- `patch_cluster.R`: Patch-based clustering
- `selfup_cluster.R`: Self-updating clustering
- `turbo_cluster.R`: Fast clustering variants

### Spatial Operations
- Coordinate conversion via `index_to_coord()` and `index_to_grid()` (from neuroim2)
- Spatial gradient computation in `R/gradient.R` using neighbor weights
- Tessellation via k-means on spatial coordinates (`tesselate()`)

## Working with Test Data
- `testdata/S1betas.rds`: Example beta coefficients
- `testdata/mask.nii`: Brain mask in NIfTI format  
- `testdata/rscan01.nii.gz`: Example 4D neuroimaging scan
- Use `neuroim2::read_vol()` and `neuroim2::read_vec()` to load test data