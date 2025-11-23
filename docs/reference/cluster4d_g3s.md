# G3S: Gradient-Guided Geodesic Supervoxels

The "Platonic Ideal" of fMRI clustering: combines manifold learning (SVD
compression), gradient-based seeding, and geodesic propagation for O(N
log N) speed with high biological plausibility.

## Usage

``` r
cluster4d_g3s(
  vec,
  mask,
  K = 100,
  n_components = 15,
  variance_threshold = 0.95,
  alpha = 0.5,
  compactness = NULL,
  max_refinement_iter = 3,
  verbose = FALSE,
  use_irlba = TRUE,
  use_rsvd = TRUE
)
```

## Arguments

- vec:

  A `NeuroVec` instance supplying the 4D neuroimaging data to cluster.

- mask:

  A `NeuroVol` mask defining voxels to include. Nonzero = included.

- K:

  Integer; target number of clusters. Default: 100.

- n_components:

  Integer; number of SVD components for feature compression.
  Default: 15. Higher values preserve more variance but reduce speed
  gains.

- variance_threshold:

  Numeric (0-1); minimum variance to preserve in SVD. If n_components
  doesn't meet this, more components will be added. Default: 0.95.

- alpha:

  Numeric (0-1); feature weight. 0 = all spatial, 1 = all feature.
  Default: 0.5 (balanced).

- compactness:

  Numeric; spatial scaling factor. Larger = more compact clusters.
  Default: auto-computed from volume and K.

- max_refinement_iter:

  Integer; number of boundary refinement iterations. Default: 3. Set to
  0 to skip refinement.

- verbose:

  Logical; print progress messages. Default: FALSE.

- use_irlba:

  Logical; use fast randomized SVD for large datasets. Default: TRUE.

- use_rsvd:

  Logical; if TRUE and the rsvd package is installed, prefer the
  randomized SVD backend in `compress_features_svd`. Default: TRUE.

## Value

A `cluster4d_result` object (also inherits `g3s_result`,
`cluster_result`) with components:

- clusvol:

  `ClusteredNeuroVol` with cluster assignments.

- cluster:

  Integer vector of cluster labels for masked voxels.

- centers:

  Matrix of cluster centers in feature space.

- coord_centers:

  Matrix of cluster spatial centers.

- n_clusters:

  Actual number of clusters produced.

- method:

  Character; "g3s".

- parameters:

  List of all parameters used.

- metadata:

  G3S-specific metadata including compression info.

## Details

### Algorithm Overview

G3S combines four key phases:

1.  **Hyper-Compression (Phase 1)**: Uses SVD to reduce T=300 timepoints
    to M=15 dimensions, achieving 20x speedup in similarity calculations
    while preserving 95%+ variance.

2.  **Gradient Seeding (Phase 2)**: Finds local minima in the functional
    gradient field, placing seeds in stable functional cores rather than
    on boundaries.

3.  **Geodesic Propagation (Phase 3)**: Uses a priority queue to grow
    clusters from seeds along paths of least resistance, ensuring
    contiguity and respecting cortical geometry.

4.  **Boundary Refinement (Phase 4)**: Polishes cluster boundaries by
    checking if surface voxels correlate better with neighboring
    clusters.

### Complexity and Performance

- **Time**: O(N log N) where N = number of voxels

- **Memory**: O(N × M) where M \<\< T (typically 15 vs 300)

- **Speedup**: 10-20x faster than iterative supervoxels

- **Quality**: Superior boundaries due to geodesic propagation

### Comparison with Other Methods

|             |           |             |            |                  |
|-------------|-----------|-------------|------------|------------------|
| **Method**  | **Speed** | **Quality** | **Memory** | **Complexity**   |
| G3S         | Fast      | Excellent   | Low        | O(N log N)       |
| Supervoxels | Slow      | Good        | High       | O(N × K × iters) |
| SNIC        | Fast      | Good        | Low        | O(N log N)       |
| FLASH3D     | Fast      | Good        | Medium     | O(N)             |

G3S advantages:

- **vs. Supervoxels**: 10-20x faster, better boundaries, no iterations

- **vs. SNIC**: Better initialization (gradient vs spatial), adaptive
  centroids

- **vs. FLASH3D**: No quantization artifacts, true correlation-based
  similarity

## See also

[`cluster4d`](cluster4d.md) for unified clustering interface.
[`compress_features_svd`](compress_features_svd.md) for SVD compression
details. [`find_gradient_seeds_g3s`](find_gradient_seeds_g3s.md) for
gradient seeding details.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example data
mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
vec <- replicate(50, NeuroVol(array(rnorm(20*20*20), c(20,20,20)),
                              NeuroSpace(c(20,20,20))), simplify=FALSE)
vec <- do.call(concat, vec)

# Basic G3S clustering
result <- cluster4d_g3s(vec, mask, K = 100)
print(result$n_clusters)
print(result$metadata$variance_explained)

# Conservative compression (preserve 98% variance)
result <- cluster4d_g3s(vec, mask, K = 100, variance_threshold = 0.98)

# More aggressive compression (faster, less accurate)
result <- cluster4d_g3s(vec, mask, K = 100, n_components = 10)

# Emphasize spatial compactness
result <- cluster4d_g3s(vec, mask, K = 100, alpha = 0.3)

# Skip boundary refinement for maximum speed
result <- cluster4d_g3s(vec, mask, K = 100, max_refinement_iter = 0)
} # }
```
