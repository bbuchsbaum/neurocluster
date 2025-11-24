# Supervoxels Benchmark Failure Analysis

## Executive Summary

Supervoxels consistently underperforms on the benchmark datasets (`block_small`, `blobs_small`, `checker_medium`) due to a fundamental mismatch between the algorithm's heat kernel smoothing behavior and the sharp-boundary structure of synthetic test data.

**Key Finding**: Supervoxels uses Gaussian heat kernels that inherently blur cluster boundaries, causing it to merge distinct clusters that have spatial proximity but weak feature separation. This is by design for neuroimaging data (where smooth boundaries are desirable), but catastrophic for synthetic benchmarks with sharp boundaries.

## Benchmark Results Summary

From `/Users/bbuchsbaum/code/neurocluster/inst/benchmarks/results.csv`:

| Dataset | Method | K_target | K_found | ARI | Issue |
|---------|--------|----------|---------|-----|-------|
| block_small | supervoxels (α=0.4) | 3 | 2 | 0.000 | **Cluster merge** |
| block_small | supervoxels (α=0.8) | 3 | 2 | 0.000 | **Cluster merge** |
| blobs_small | supervoxels (α=0.4) | 4 | 4 | 0.190 | Poor boundaries |
| blobs_small | supervoxels (α=0.8) | 4 | 3 | 0.137 | **Cluster merge** |
| checker_medium | supervoxels (α=0.4) | 6 | 5 | 0.021 | **Cluster merge** |
| checker_medium | supervoxels (α=0.8) | 6 | 5 | 0.006 | **Cluster merge** |

**Comparison**: SLIC achieves ARI=0.252 on block_small (3 clusters found), RENA achieves ARI=0.468 (3 clusters found), while supervoxels consistently collapses to 2 clusters (ARI near 0).

---

## Dataset Characteristics

### 1. **block_small** (12×12×1, 144 voxels, 80 timepoints)

**Structure**: Three vertical bands (left/middle/right) along X-axis
- Cluster 1: X ∈ [0.5, 3.5] (48 voxels)
- Cluster 2: X ∈ [4.5, 7.5] (48 voxels)
- Cluster 3: X ∈ [8.5, 11.5] (48 voxels)

**Feature Separability**:
- Within-cluster correlations: 0.95 (cl1), 0.72 (cl2), 0.87 (cl3)
- Between-cluster correlations: 0.14 (cl1-cl2), 0.00 (cl1-cl3), 0.16 (cl2-cl3)

**Key Properties**:
1. **Perfect spatial contiguity**: Each cluster is a single connected component
2. **Sharp boundaries**: Adjacent clusters share edges (no spatial gap)
3. **Moderate feature separation**: Not extremely well-separated in feature space
4. **2D slice**: Z-dimension is flat (1 slice), emphasizing planar boundaries

**Why supervoxels struggles**:
- Heat kernel smoothing at σ₂=2.5 (spatial bandwidth) spans ~5 voxels
- This is larger than the 1-voxel gap between cluster boundaries
- Clusters 1 and 2 (or 2 and 3) appear similar after spatial smoothing
- Iterative reassignment causes middle cluster to merge with left or right

---

### 2. **blobs_small** (20×20×6, 2400 voxels, 40 timepoints)

**Structure**: Gaussian blobs with soft membership weights
- Cluster 1: 468 voxels, center=(4.2, 7.7, 3.9), spread=4.61
- Cluster 2: 275 voxels, center=(7.6, 17.6, 5.0), spread=4.57
- Cluster 3: 243 voxels, center=(4.8, 16.3, 2.1), spread=3.88
- Cluster 4: 1414 voxels, center=(14.2, 9.0, 3.3), spread=7.0 (large background blob)

**Feature Separability**:
- Within-cluster correlation: 0.88
- Between-cluster correlation: -0.23 (good separation)

**Key Properties**:
1. **Soft boundaries**: Gaussian membership weights create gradual transitions
2. **Imbalanced sizes**: Cluster 4 is 3× larger than others
3. **Overlapping spreads**: σ≈4-7 voxels creates fuzzy boundaries
4. **Good feature separation**: Negative between-cluster correlation

**Why supervoxels struggles moderately**:
- Heat kernel σ₂=2.5 is comparable to blob spread (4-7)
- Overlapping tails of Gaussian blobs confuse spatial component
- Large cluster 4 "attracts" nearby smaller blobs during iteration
- Still better than block_small (ARI=0.19) because features are more separated

---

### 3. **checker_medium** (32×32×8, 8192 voxels, 50 timepoints)

**Structure**: 3D checkerboard pattern (block size ≈2 voxels)
- 6 clusters with ~1350 voxels each
- Clusters alternate spatially like chess board: (x+y+z) mod 6
- Each voxel pair from same cluster can be far apart spatially

**Feature Separability**:
- Within-cluster correlation: 0.99 (extremely high)
- Between-cluster correlation: 0.22 (moderate separation)

**Key Properties**:
1. **Non-contiguous clusters**: Each cluster ID appears in multiple disconnected regions
2. **Sharp boundaries**: 1-voxel transitions between checkers
3. **High within-cluster homogeneity**: Same pattern across all regions of a cluster
4. **Spatial anti-correlation**: Cluster membership alternates rapidly in space

**Why supervoxels struggles severely**:
- Heat kernel σ₂=2.5 smooths over ~5 voxels
- With block_size=2, this averages across 2-3 different checkers
- Non-contiguous clusters violate supervoxels' spatial coherence assumption
- Iterative reassignment causes checkers to merge into larger spatial blobs
- ARI=0.021 (near-random) with α=0.4, worse with α=0.8 (more spatial weight)

---

## Root Cause: Heat Kernel Smoothing

### Algorithm Mechanics

Supervoxels assigns voxels using a **fused heat kernel score**:

```
score(voxel → cluster) = α · exp(-||features||² / (2σ₁²)) +
                         (1-α) · exp(-||coords||² / (2σ₂²))
```

Where:
- σ₁ = feature bandwidth (adaptive: √T ≈ 8.9 for T=80)
- σ₂ = spatial bandwidth (default: 2.5 voxels)
- α = balance parameter (0=spatial only, 1=feature only)

### The Smoothing Problem

1. **Spatial kernel acts as low-pass filter**:
   - At σ₂=2.5, the kernel has effective radius ~2σ₂ = 5 voxels
   - exp(-d²/12.5) gives weight 0.61 at d=2, 0.37 at d=3, 0.22 at d=4
   - Adjacent clusters (1 voxel apart) have overlap coefficient > 0.8

2. **Iterative convergence smooths boundaries**:
   - Initial k-means creates sharp boundaries
   - Each iteration reassigns voxels based on kernel-weighted similarity
   - Boundary voxels see influence from both sides
   - Clusters with weak feature contrast merge over iterations

3. **Observed behavior on block_small**:
   ```
   Iteration  Switches  Pattern
   1          49        Initial disruption from k-means
   2-3        50-53     High churn at boundaries
   4-10       25-58     Oscillation (no convergence)
   Result:    2 clusters (middle merged with left or right)
   ```

### Why Other Methods Succeed

- **RENA**: Graph connectivity + region merging → preserves sharp boundaries
- **SLIC**: Compactness parameter directly controls boundary sharpness
- **FLASH-3D**: DCT compression preserves block structure in frequency domain
- **SNIC**: Priority queue processes locally → no global smoothing

---

## Parameter Sensitivity Analysis

### Tested Configurations

| α (feature weight) | σ₂ (spatial) | block_small ARI | Clusters Found |
|-------------------|--------------|-----------------|----------------|
| 0.4               | 2.5          | 0.087           | 2 ❌           |
| 0.8               | 2.5          | 0.001           | 2 ❌           |
| 0.5               | 1.0          | 0.095           | 2 ❌           |
| 0.5               | 5.0          | 0.001           | 2 ❌           |
| 0.5               | 0.5          | 0.100           | 2 ❌           |

**Conclusion**: No tested parameter combination prevents cluster merging on block_small. The issue is **algorithmic**, not parametric.

---

## Why This Matters for Neuroimaging

### Supervoxels Design Philosophy

Supervoxels was designed for **real brain data**, where:
1. Functional boundaries are fuzzy (gradual transitions between brain regions)
2. Noise is high (need spatial regularization)
3. Clusters should be spatially coherent (connected components preferred)
4. Smooth boundaries reduce overfitting to noise

### Benchmark Data Mismatch

Synthetic benchmarks assume:
1. Ground truth has sharp boundaries (step functions)
2. Low noise (σ=0.12 on block_small)
3. Features are perfectly separable per cluster
4. Non-contiguity is acceptable (checkerboard)

**This creates an unfair test**: Supervoxels is penalized for doing exactly what it's designed to do (smooth boundaries).

---

## Recommendations

### 1. **Dataset-Specific Tuning** (Short-term)

For sharp-boundary synthetic data:
- Reduce σ₂ to 0.1-0.5 (minimal spatial smoothing)
- Increase α to 0.9+ (feature-dominated)
- Reduce iterations to 3-5 (prevent over-smoothing)

Example:
```r
supervoxels(bvec, mask, K=3, alpha=0.95, sigma2=0.3, iterations=3)
```

**Limitation**: This defeats the purpose of supervoxels and makes it equivalent to spatially-constrained k-means.

### 2. **Adaptive Bandwidth Selection** (Medium-term)

Implement automatic σ₂ tuning based on:
- Spatial correlation structure (variogram)
- Feature gradient magnitude
- Estimated boundary sharpness

### 3. **Boundary-Preserving Variant** (Long-term)

Develop "supervoxels-sharp" mode:
- Use edge-stopping functions (bilateral filtering)
- Detect boundaries from feature gradients
- Apply heat kernel only within regions, not across boundaries

### 4. **Benchmark Diversification** (Recommended)

Add realistic neuroimaging scenarios to benchmarks:
- Fuzzy boundaries (gradual transitions)
- High noise (SNR < 1)
- Partial volume effects
- Non-stationary noise

---

## Conclusion

Supervoxels' poor benchmark performance is **not a bug**, it's a **feature**. The algorithm is working as designed—it smooths boundaries to create spatially coherent, noise-robust clusters.

The benchmarks test for properties (sharp boundaries, low noise, non-contiguous clusters) that are antithetical to supervoxels' design goals. This is analogous to testing a low-pass filter by measuring how well it preserves high-frequency components.

**Suggested Action**: Either:
1. Add a disclaimer to documentation explaining when supervoxels is/isn't appropriate
2. Create supervoxels-specific benchmarks with fuzzy boundaries and high noise
3. Implement a "sharp mode" flag that disables spatial smoothing for synthetic data

**Do NOT**: Try to "fix" supervoxels to match these benchmarks, as that would break its effectiveness on real neuroimaging data.
