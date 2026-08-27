# Compress High-Dimensional Features via Randomized SVD

Reduces the dimensionality of fMRI time series data using Singular Value
Decomposition (SVD) while preserving the majority of variance. This is a
key component of the G3S (Gradient-Guided Geodesic Supervoxels)
algorithm, enabling 10-20x speedup in similarity computations.

## Usage

``` r
compress_features_svd(
  feature_mat,
  n_components = 15,
  variance_threshold = 0.95,
  use_irlba = TRUE,
  use_rsvd = TRUE
)
```

## Arguments

- feature_mat:

  A numeric matrix with voxels in rows and timepoints in columns (N x
  T). The function will center and scale the data automatically.

- n_components:

  Target number of dimensions (default: 15). Higher values preserve more
  variance but reduce speed gains.

- variance_threshold:

  Minimum proportion of variance to preserve (0-1). If the requested
  n_components doesn't meet this threshold, more components will be
  retained. Default: 0.95 (95% variance).

- use_irlba:

  Logical; if TRUE and the irlba package is available, use fast
  randomized SVD for large datasets (\>10,000 voxels). Default: TRUE.

- use_rsvd:

  Logical; if TRUE and the rsvd package is available, prefer the
  randomized SVD implementation from rsvd for large datasets (\>10,000
  voxels). This can outperform irlba on tall-and-skinny matrices.
  Default: TRUE.

## Value

A list with components:

- features:

  Compressed feature matrix (N x M), starting from the requested
  `n_components` and increasing M when needed to satisfy
  `variance_threshold`. Each row is normalized to unit length for cosine
  similarity calculations.

- variance_explained:

  Proportion of total variance explained by the retained components
  (0-1).

- n_components:

  Actual number of components retained (may be \> n_components input if
  needed to meet variance_threshold).

- rotation:

  The right singular vectors (V matrix) for transforming new data.

- singular_values:

  The singular values (d vector) for each component.

- center:

  Column means used for centering (for transforming new data).

- scale:

  Column standard deviations used for scaling (for transforming new
  data).

## Details

### Algorithm

The compression follows these steps:

1.  **Center and scale**: Each timepoint is centered (mean=0) and scaled
    (sd=1) to ensure equal contribution from all timepoints.

2.  **SVD decomposition**: Computes `X = U * D * t(V)` where:

    - U (N x k): Left singular vectors (voxel loadings)

    - D (k x k): Diagonal matrix of singular values

    - V (T x k): Right singular vectors (time loadings)

3.  **Compression**: Retains only the first M components where M is
    chosen to balance speed (fewer components) and accuracy (more
    variance explained).

4.  **Normalization**: Each compressed feature vector is normalized to
    unit length, enabling fast cosine similarity via dot products:
    `cor(x, y) approx x * y`

### Performance Characteristics

- **Memory**: Compressed data uses M/T of original size (e.g., 15/300 =
  5%)

- **Speed**: Similarity calculations are M/T times faster (e.g., 20x for
  15/300)

- **Accuracy**: Typically preserves 95%+ of signal with M=10-20 for fMRI
  data

### Choosing n_components

- **Small (M=10)**: Maximum speed, good for noisy or heavily smoothed
  data

- **Medium (M=15)**: Balanced, recommended default for most fMRI data

- **Large (M=25)**: High accuracy, use for high SNR data or critical
  applications

The function will automatically increase n_components if needed to meet
variance_threshold. If that requires the full available rank, the exact
full-rank decomposition is used.

## Reproducibility

Both randomized backends draw random projections. They are run under a
fixed internal seed and the caller's RNG state is restored afterwards,
so repeated calls on the same data return the same decomposition and the
function never perturbs the surrounding stream. Without that, two
identical
[`cluster4d_g3s`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_g3s.md)
runs could disagree on half their voxel labels.

## See also

[`cluster4d_g3s`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_g3s.md)
for the full G3S clustering algorithm.
[`transform_new_data_svd`](https://bbuchsbaum.github.io/neurocluster/reference/transform_new_data_svd.md)
for applying the compression to new data.

## Examples

``` r
# Simulate fMRI time series data
n_voxels <- 1000
n_timepoints <- 300
latent <- matrix(rnorm(n_voxels * 8), n_voxels, 8)
loadings <- matrix(rnorm(n_timepoints * 8), n_timepoints, 8)
feature_mat <- latent %*% t(loadings) +
  matrix(rnorm(n_voxels * n_timepoints, sd = 0.05), n_voxels, n_timepoints)

# Compress to 15 dimensions (typical for G3S)
compressed <- compress_features_svd(feature_mat, n_components = 15)
print(compressed$variance_explained)  # At least 0.95
#> [1] 0.9996279
print(ncol(compressed$features))      # 15 here; may increase when necessary
#> [1] 15

# More aggressive compression
compressed_fast <- compress_features_svd(feature_mat, n_components = 10)

# Conservative compression (preserve 98% variance)
compressed_accurate <- compress_features_svd(
  feature_mat,
  n_components = 15,
  variance_threshold = 0.98
)
```
