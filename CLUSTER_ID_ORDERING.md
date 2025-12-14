# Cluster ID Ordering in Neurocluster

## The Problem

Many clustering algorithms in `neurocluster` assign cluster IDs sequentially as clusters are encountered during processing. When voxels are processed in **memory order** (which is typically **z-major** in 3D arrays due to R's column-major storage), this creates systematic spatial ordering:

**Cluster ID 1, 2, 3, ... → Inferior to Superior brain regions**

### Affected Algorithms

| Method | Affected? | Spearman r_z | Mechanism |
|--------|-----------|--------------|-----------|
| `flash3d` | ✅ **YES** | **0.825** | Coarse cell iteration in z-order + sequential ID assignment |
| `slice_msf` | ✅ **YES** | **0.885** | Final label assignment in z-major linear order |
| `supervoxels` | ✅ **NO** | 0.013 | Iterative refinement breaks initial z-ordering |
| `snic` | ✅ **NO** | -0.193 | Gradient-based seed selection not z-ordered |
| `slic4d` | ⚠️ **Unknown** | — | Not yet tested |
| `acsc` | ⚠️ **Unknown** | — | Not yet tested |

**Note**: Correlations measured empirically on test data (96×96×26 voxels, K=50). |r_z| > 0.3 indicates ordering.

## Why This Matters

### 1. Visualization Artifacts

When viewing clusters with **continuous color scales** (e.g., viridis, rainbow), ordered cluster IDs create artificial gradients that suggest spatial structure where none exists:

```r
# With ordered IDs (BAD)
result <- cluster4d(vec, mask, K = 100, method = "slice_msf")
plot(result$clusvol)
# Shows smooth color gradient from blue (inferior) → yellow (superior)
# This is ARTIFICIAL - not real brain structure!
```

### 2. Interpretation Bias

Ordered IDs may mislead users into thinking cluster numbering has anatomical significance:
- "Cluster 1 is always inferior ventral regions"
- "Cluster 100 is always superior dorsal regions"

This is **incorrect** - the ordering is a computational artifact.

### 3. Statistical Issues

Some analyses may inadvertently treat cluster ID as a continuous variable, which is meaningless if IDs are ordered.

---

## Solutions

### Option 1: Use `randomize_cluster_ids()` (Recommended)

Apply cluster ID randomization to any clustering result:

```r
library(neurocluster)

# Standard clustering
result <- cluster4d(vec, mask, K = 100, method = "slice_msf")

# Randomize cluster IDs
result_random <- randomize_cluster_ids(result, seed = 42)

# Now visualize with continuous color scale
plot(result_random$clusvol)  # Distinct color patches, no gradient
```

**Advantages**:
- Works with any clustering method
- Preserves all cluster properties (centroids, sizes, etc.)
- Reproducible with seed parameter
- Can be applied post-hoc to saved results

### Option 2: Check Before Randomizing

Use `check_cluster_ordering()` to detect if randomization is needed:

```r
# Check for spatial ordering
ordering <- check_cluster_ordering(result)
print(ordering)
#> $cor_z
#> [1] 0.87  # Strong correlation with z-axis!
#>
#> $is_ordered
#> [1] TRUE
#>
#> $interpretation
#> [1] "Strong spatial ordering detected - randomization strongly recommended"

# Randomize if needed
if (ordering$is_ordered) {
  result <- randomize_cluster_ids(result)
}
```

### Option 3: Fix at Algorithm Level

For algorithms where ordering is undesirable, modify the source code to randomize seed placement or label assignment order.

**Example**: FLASH-3D now uses spatial-balanced seed selection instead of first-come-first-served (see `FLASH3D_FIX_SUMMARY.md`).

---

## Best Practices

### For Visualization

**Always use discrete color scales** or randomize IDs:

```r
# GOOD: Discrete color palette
library(neuroim2)
library(RColorBrewer)
discrete_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_clusters)
plot(result$clusvol, col = discrete_colors)

# GOOD: Randomized IDs with continuous scale
result <- randomize_cluster_ids(result)
plot(result$clusvol)

# BAD: Ordered IDs with continuous scale
plot(result_ordered$clusvol)  # Creates artificial gradient!
```

### For Publications

1. **Always mention** whether cluster IDs were randomized
2. **Report** the correlation between cluster ID and spatial position
3. **Use** discrete color schemes for cluster visualizations
4. **Include** the randomization seed for reproducibility

```r
# Example for methods section:
result <- cluster4d(vec, mask, K = 100, method = "flash3d")
result <- randomize_cluster_ids(result, seed = 42)
# "Cluster IDs were randomized to avoid visualization artifacts
#  (seed = 42 for reproducibility)"
```

### For Analysis Pipelines

Add randomization as a standard post-processing step:

```r
# Standard pipeline
cluster_and_process <- function(vec, mask, K, method = "flash3d") {
  # Cluster
  result <- cluster4d(vec, mask, K, method = method)

  # Always randomize
  result <- randomize_cluster_ids(result)

  # Continue with analysis
  metrics <- compute_cluster_metrics(result)

  list(result = result, metrics = metrics)
}
```

---

## Technical Details

### Why Z-Major Ordering Occurs

1. **R's memory layout**: 3D arrays stored in column-major order
   - `array[x, y, z]` stored as: all x for (y=1, z=1), then all x for (y=2, z=1), ...
   - **x varies fastest, z varies slowest**

2. **Linear indexing**: `which(mask > 0)` returns indices in storage order
   - Index 1: x=1, y=1, z=1
   - Index 2: x=2, y=1, z=1
   - ...
   - Index N: x=max, y=max, z=max

3. **Sequential processing**: Most algorithms iterate voxels in index order
   - First voxels encountered → lower z-slices
   - Later voxels encountered → higher z-slices

4. **Label assignment**: First unique cluster encountered gets ID 1
   - Lower z-slices processed first → get lower cluster IDs
   - Higher z-slices processed later → get higher cluster IDs

### Spearman Correlation Thresholds

| |r| | Interpretation | Action |
|------|----------------|---------|
| > 0.7 | Strong ordering | **Must randomize** |
| 0.3-0.7 | Moderate ordering | **Should randomize** |
| < 0.3 | Weak ordering | Randomization optional |

---

## Algorithm-Specific Notes

### FLASH-3D

**Status**: ⚠️ **PARTIALLY FIXED** - spatial balance improved but cluster ID ordering remains

**Current behavior**:
- Seeds are spatially balanced (no bias toward any z-slice)
- However, seeds are selected in z-order (coarse cells iterated z-first)
- Cluster IDs assigned sequentially: seeds[0] → ID 0, seeds[1] → ID 1, etc.
- **Result**: Strong z-correlation (r_z = 0.825)

**Root cause**: Lines 183-707 in `src/flash3d.cpp`:
1. Coarse cells numbered z-first: `iz = c / (gx*gy)`
2. Loop iterates: `for (c = 0; c < gcells; ++c)`
3. Seeds accumulated in z-order
4. Labels assigned: `label[g] = k` for k = 0..K-1

**Solution**: Always use `randomize_cluster_ids()` after clustering:
```r
result <- supervoxels_flash3d(vec, mask, K = 100)
result <- randomize_cluster_ids(result)
```

### slice_msf

**Status**: ⚠️ **NOT FIXED** - produces very strong z-ordering (r_z = 0.885)

**Mechanism**: Final label assignment (`src/slice_cluster.cpp:791-797`) iterates voxels in z-major linear order:
```cpp
for (int v=0; v<N; ++v) {  // v in linear (z-major) order
  if (root2lab.find(rG) == root2lab.end()) {
    root2lab.emplace(rG, lab);
    labels[v] = lab;
    lab++;  // Sequential assignment
  }
}
```
First cluster encountered → label 1, next → label 2, etc.
Since low-z voxels processed first, they get low cluster IDs.

**Solution**: Always use `randomize_cluster_ids()` after clustering:
```r
result <- slice_msf(vec, mask, target_k_global = 100)
result <- randomize_cluster_ids(result)
```

**Alternative**: Could randomize voxel iteration order in C++, but post-processing is simpler and safer.

---

## FAQ

**Q: Does randomization change my clustering results?**
A: No. It only changes the numerical IDs assigned to clusters. All voxel assignments, centroids, and cluster properties remain identical.

**Q: Should I randomize before or after computing cluster metrics?**
A: Doesn't matter - metrics are invariant to ID labels. But randomize before visualization.

**Q: Can I use the same randomization across subjects?**
A: No - each subject's clustering is independent. Use the same seed for reproducibility within a subject, but each subject will have different cluster-to-ID mappings.

**Q: Will this be fixed in all algorithms?**
A: FLASH-3D is fixed at the algorithm level. For others, `randomize_cluster_ids()` provides a simple, universal solution that works across all methods.

**Q: Does this affect cluster correspondence across subjects?**
A: Cluster correspondence across subjects requires alignment methods (e.g., template matching, functional alignment) regardless of ID ordering. Randomization doesn't affect this.

---

## References

- **FLASH-3D fix**: See `FLASH3D_FIX_SUMMARY.md` for details on spatial-balanced seed selection
- **Memory layout**: R Internals manual, Section 1.1.4 "Attributes"
- **Column-major order**: Golub & Van Loan (2013), Matrix Computations, Section 1.2.6

### supervoxels

**Status**: ✅ **NO Z-ORDERING** (r_z = 0.013)

**Why not affected**: Although initialization uses evenly-spaced z-ordered seed points, iterative refinement scrambles the initial ordering:
- K-means initialization on spatial coordinates creates z-ordered clusters
- However, 10+ iterations of heat kernel-based reassignment break this pattern
- Final clusters show no systematic z-correlation

**Verification**: Empirically tested - no randomization needed for supervoxels.

### SNIC

**Status**: ✅ **NO Z-ORDERING** (r_z = -0.193)

**Why not affected**: Gradient-based seed selection does not produce z-ordered seeds:
- Seeds selected at high-gradient locations (edges, boundaries)
- Gradient magnitude uncorrelated with z-position
- Priority queue processing maintains seed order but seeds themselves not z-ordered

**Verification**: Empirically tested - no randomization needed for SNIC.

---

**Last updated**: 2025-01-24
**Package version**: 0.1.0
**Verification**: Empirical testing on 96×96×26 test data (56,542 voxels, K=50)
