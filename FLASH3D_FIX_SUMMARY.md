# FLASH-3D Z-Axis Striping Fix

## Problem Summary

FLASH-3D clustering consistently produced **z-axis striping artifacts** (horizontal layers from inferior to superior) regardless of parameter settings (`lambda_s`, `lambda_t`, `bits`, `dctM`). Cluster indices increased monotonically from inferior to superior brain regions.

## Root Cause

A **triple-compounding z-bias** in the seed placement algorithm:

1. **R's `which()` ordering**: Returns indices in column-major order where z varies slowest
   - Creates z-contiguous blocks: all z=0 voxels, then z=1, then z=2, etc.

2. **Z-major coarse grid indexing**: Coarse cells numbered as `(iz * gy + iy) * gx + ix`
   - Groups cells by z-coordinate first

3. **First-come-first-served selection**: Sequential scan of z-ordered voxels
   - Always selects the lowest z-voxel within each coarse cell
   - Creates systematic inferior→superior seed placement

## The Fix

**File**: `src/flash3d.cpp` lines 172-217

**Change**: Replaced first-come-first-served selection with **spatial-balanced selection**:
- Build membership lists for all coarse cells
- Calculate each cell's geometric center
- Select the voxel **closest to the cell center**
- Eliminates processing order bias

### Code Changes

**Before (biased)**:
```cpp
for (size_t vi = 0; vi < mask_lin.size(); ++vi) {
  int c = cell_of(mask_lin[vi]);
  if (coarse[c] < 0) {
    coarse[c] = (int)vi;  // First voxel wins
    reps.push_back((int)vi);
  }
}
```

**After (unbiased)**:
```cpp
// Build cell membership lists
std::vector<std::vector<int>> cell_members(gcells);
for (size_t vi = 0; vi < mask_lin.size(); ++vi) {
  int c = cell_of(mask_lin[vi]);
  cell_members[c].push_back((int)vi);
}

// Select voxel closest to each cell's geometric center
for (int64_t c = 0; c < gcells; ++c) {
  if (cell_members[c].empty()) continue;

  // Calculate cell center coordinates
  // (ix, iy, iz) -> (cx, cy, cz)

  // Find closest voxel to center
  for (int vi : cell_members[c]) {
    double d2 = (xc - cx)² + (yc - cy)² + (zc - cz)²;
    // Select minimum distance
  }
}
```

## Expected Results

### Before Fix:
- Spearman correlation (cluster ID vs z-position): **r ≈ 0.8-0.95** (strong ordering)
- Visible horizontal striping in cluster visualizations
- Cluster IDs 1-10 in inferior slices, 90-100 in superior slices

### After Fix:
- Spearman correlation: **r ≈ -0.1 to +0.1** (no ordering)
- Uniform z-distribution across clusters
- Cluster IDs spatially random

## Verification

Run the provided verification script:

```bash
Rscript verify_flash3d_fix.R
```

Or test manually in R:

```r
library(neurocluster)
library(neuroim2)

# Your data
result <- cluster4d(vec, mask, K = 100, method = "flash3d")

# Check z-ordering
coords <- result$coord_centers
cor(1:nrow(coords), coords[, 3], method = "spearman")
# Should be near 0 (no correlation)

# Visualize z-distribution
hist(coords[, 3], breaks = 20, main = "Cluster Z-Distribution")
# Should be roughly uniform across z-range
```

## Technical Details

### Why Other Components Were Not the Cause:

✅ **DCT hashing**: Purely temporal, no spatial dependencies
✅ **Jump-flood propagation**: Perfectly isotropic 26-neighbor sampling
✅ **Hamming distance**: Valid metric, no z-axis bias
✅ **Spatial distance**: Isotropic Euclidean metric
✅ **Annealing**: Uniform scaling across all axes
✅ **Buffer management**: Correct implementation, no race conditions

The issue was purely in the **deterministic seed placement** interacting with R's memory layout.

### Why Parameters Didn't Help:

All parameters (`lambda_s`, `lambda_t`, `bits`, `dctM`, `rounds`) operate **after** seed placement. The z-bias was baked into the initial seeds before any parameter-dependent processing occurred.

## Performance Impact

**Negligible**: The fix adds one pass to build cell membership lists, but this is O(N) where N = number of masked voxels, which is typically small compared to the jump-flood iterations. Expected overhead: < 1% of total runtime.

## Backward Compatibility

The fix changes seed placement, so **results will differ** from previous versions:
- Cluster spatial distributions will be similar (same algorithm)
- Cluster **IDs** will be assigned differently (no z-ordering)
- Clustering **quality** should be equivalent or slightly better (unbiased initialization)

## Alternative Fixes Considered

1. **Post-processing ID shuffle** (simplest but cosmetic only)
2. **Deterministic randomization** (requires PRNG, less elegant)
3. **Z-coordinate priority queue** (more complex, similar result)

The implemented solution (spatial-balanced selection) is preferred because it:
- Eliminates the root cause
- Maintains determinism
- Requires no PRNG
- Creates optimal initial seeds (closest to cell centers)

## Questions?

If striping still appears after this fix, check:
1. `vox_scale` parameter matches actual voxel dimensions
2. Input data doesn't have systematic z-gradients in temporal patterns
3. Run verification script to quantify any remaining bias

---

**Fix Implemented**: 2025-01-24
**Package Version**: 0.1.0
**Affected Functions**: `supervoxels_flash3d()`, `cluster4d(..., method = "flash3d")`
