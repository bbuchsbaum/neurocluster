# Clustering Methods Audit (Correctness + Performance)

This note captures a code-level audit of the main clustering algorithms in **neurocluster** and the highest-impact correctness/soundness and speed opportunities.

## Build/Runtime Soundness

- **RcppParallel linkage**: `src/Makevars` and `src/Makevars.win` previously omitted `RcppParallel::RcppParallelLibs()`. On systems where `RcppParallel.so` is not already loaded globally, this can manifest as a dynamic loading failure (missing `RcppParallel::tbbParallelFor` symbol) when the package DLL is loaded.
  - Fix: link against the RcppParallel-provided TBB libraries via `RcppParallel::RcppParallelLibs()`.

## Algorithm-Specific Findings

### SNIC (`R/snic.R`, `src/snic.cpp`)

**Correctness/soundness issues**

- **Coordinate scaling bug**: `index_to_grid()` returns voxel/grid coordinates; these were previously **divided** by `spacing(mask)` to form `norm_coords`. Since `spacing()` is **mm/voxel**, correct conversion to physical coordinates is to use `index_to_coord()` (or multiply by spacing), not divide.
- **Spatial scale bug**: the expected seed spacing for a 3D volume should scale like `(volume / K)^(1/3)` (or `(N/K)^(1/3)` in voxel-units). The prior `sqrt(N/K)` heuristic is a 2D spacing rule and can mis-weight the spatial term substantially in 3D.
- **Result coordinate units**: coordinate centers should be consistently in physical coordinates (mm) across methods; SNIC previously reported centers in voxel units via `index_to_grid()`.

**Performance issues**

- **Dead computation**: a sample-pairwise `dist()` statistic (`sf`) was computed but never used.
- **Extra conversions**: initializing `L <- array(0, dim(mask))` forces coercion to integer inside Rcpp; using `0L` avoids avoidable copies.

**Implemented changes**

- Use `index_to_coord(mask, mask.idx)` for spatial distances and for computing/reporting `coord_centers`.
- Compute `S` from the physical bounding box and use `s = (1.1 * S)^2` because the C++ SNIC distance uses squared spatial distances divided by `s`.
- Remove unused `dist()` statistic.
- Initialize `L` as an integer array (`0L`).

### G3S (`src/g3s_core.cpp`, `R/g3s.R`)

**Correctness/soundness**

- `g3s_propagate_cpp()` assumed seed indices are always valid and unique and indexed `centroids[k]` based on the seed loop index. If invalid/duplicate seeds ever appear, this can cause **out-of-bounds access** (undefined behavior).

**Implemented change**

- Make the seed initialization robust by assigning labels contiguously as centroids are appended and always using `centroids.back()` when seeding neighbor pushes; skip duplicates defensively.

### SVD Compression (`R/g3s_compression.R`)

**Performance**

- `U %*% diag(d)` allocates a diagonal matrix. This is avoidable.

**Implemented change**

- Replace `U %*% diag(d)` with `sweep(U, 2, d, "*")`.

### Center Computation (`R/cluster4d_common.R`)

**Performance**

- `compute_cluster_centers(..., method="mean")` previously looped over clusters and called `colMeans()` repeatedly, which becomes a bottleneck for large `N` and/or large `K`.

**Implemented change**

- For `"mean"`, compute group sums via `rowsum()` (C-accelerated) and divide by counts.

### SLIC 4D (`src/slic4d.cpp`, `R/slic4d.R`)

**Correctness/soundness**

- **Thread-safety bug causing segfaults**: the center-grid index stored bbox minima in an `Rcpp::NumericVector` and accessed it from within `RcppParallel::parallelFor()`. Rcpp objects are not safe to access concurrently, and this can manifest as intermittent crashes.
- **Incorrect distance upper-bound reuse**: the assignment kernel initialized the per-voxel best distance from the prior iteration (`dists[i]`). Center locations/features change between iterations, so those distances are not comparable and can prevent label updates (and can leave `dists` inconsistent with the current centers).

**Implemented change**

- Store bbox minima as plain doubles in the grid index, and add stricter input-shape validation at the start of `slic4d_core()`.
- Recompute best distances per iteration (ignore prior `dists[i]`), and compute feature distances via cached norms + dot products (`||x-c||^2 = ||x||^2 + ||c||^2 - 2 x·c`) with a per-thread cached voxel feature row.

### Supervoxels / Fused Assignment (`src/fused_assignment_binned.cpp`, `R/supervoxels.R`)

**Performance**

- The fused-assignment inner loop previously computed `sum((x-c)^2)` directly, which is still O(D) but does extra scalar work per feature per candidate.

**Implemented change**

- Compute feature distances via cached norms + dot products with cached centroid norms (`||c||^2`), reducing arithmetic per feature per candidate while preserving exact results.

### ACSC (`src/acsc_boundary.cpp`, `R/acsc.R`)

**Correctness/soundness**

- **Strided-memory bug** (boundary refinement): the C++ boundary refinement treated each voxel time series as contiguous in memory, but R matrices are column-major. This made the dot-products (and centroid accumulation) incorrect, reducing refinement effectiveness.
- **Neighbor definition**: boundary refinement used `FNN::get.knn(coords, k=6)` as a proxy for 6-connected adjacency; this is slower and can connect across gaps in sparse masks.

**Performance**

- **kNN hotspot**: `FNN::get.knn()` on block summaries is a dominant cost for large numbers of blocks and long time-series.

**Implemented change**

- Fix dot-products and centroid accumulation to respect R’s column-major layout (strided access) and replace map-based centroids with a dense centroid table for speed.
- Use grid-based 6-connected neighbor indices in R, avoiding `get.knn()` for boundary refinement.
- Add an optional block-summary projection (`knn_proj_dim`, `knn_proj_method`) before kNN graph construction (DCT basis or random projection), reducing kNN dimension when timepoints are large.

## Cross-Method Opportunities (Still Open)

### ReNA (`R/rena.R`, `src/rena.cpp`)

**Performance**

- Extracting sparse edges via `Matrix::which(current_adjacency > 0, arr.ind=TRUE)` forces work proportional to the dense matrix size.

**Implemented change**

- Use `Matrix::summary(current_adjacency)` to pull nonzeros directly (sparse), then filter to upper-triangular edges.

- **Avoid `as.array(NeuroVec)` for gradients**: `slic4d_supervoxels()` and some initialization paths materialize full 4D arrays just to compute mean/gradients. For fMRI-sized data this can dominate both time and memory; prefer masked operations via `series()` and local neighborhoods when possible.
- **Empty-cluster handling in iterative methods**: `supervoxel_cluster_fit()` uses a fixed `K` and can produce empty clusters; empty centroids can silently become attractors (e.g., if left at zeros). A robust strategy is to track cluster sizes and skip/bin only active centroids or re-seed empties.
- **Coordinate-unit consistency**: some methods report centers in voxel units while others use mm coordinates. Standardizing on mm (and documenting it) simplifies downstream comparison and plotting.
- **kNN cost**: repeated `FNN::get.knn()` calls can dominate runtime; cache neighbors when possible and avoid recomputation across refinement stages.

## Files Changed

- `src/Makevars`
- `src/Makevars.win`
- `R/snic.R`
- `src/g3s_core.cpp`
- `R/g3s_compression.R`
- `R/cluster4d_common.R`
- `src/slic4d.cpp`
- `src/fused_assignment_binned.cpp`
- `R/acsc.R`
- `src/acsc_boundary.cpp`
- `R/rena.R`
- `tests/testthat/test_acsc_boundary_cpp.R`
- `tests/testthat/test_acsc_knn_projection.R`
- `tests/testthat/test_cluster_centers.R`
- `tests/testthat/test_fused_assignment_binned.R`
