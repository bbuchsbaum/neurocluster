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

### Key Clustering Algorithms

1. **Supervoxels** (`supervoxels()` in `R/supervoxels.R`)
   - Primary algorithm for spatially constrained clustering
   - Iterative algorithm balancing feature similarity and spatial proximity
   - Uses heat kernels with bandwidths `sigma1` (features) and `sigma2` (coordinates)
   - Alpha parameter weights data vs spatial similarity

2. **SNIC** (`snic()` in `R/snic.R`) 
   - Simple Non-Iterative Clustering algorithm
   - Uses priority queue approach with compactness parameter
   - Gradient-based initialization via `find_initial_points()`

3. **Meta-clustering** (`meta_clust()` in `R/meta_clust.R`)
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
All clustering functions return `cluster_result` objects with:
- `clusvol`: `ClusteredNeuroVol` with spatial cluster assignments
- `cluster`: Integer vector of cluster IDs
- `centers`: Feature space centroids
- `coord_centers`: Spatial coordinate centroids

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