# Generate synthetic temporal patterns for cluster signals

Creates distinct temporal patterns for each cluster. Used internally by
[`generate_synthetic_volume`](https://bbuchsbaum.github.io/neurocluster/reference/generate_synthetic_volume.md)
but exported for custom synthetic data generation.

## Usage

``` r
synthetic_time_patterns(
  n_clusters,
  n_time,
  pattern_types = c("orthogonal", "harmonic", "random")
)
```

## Arguments

- n_clusters:

  Number of clusters/patterns to generate.

- n_time:

  Number of time points.

- pattern_types:

  Type of patterns: `"orthogonal"` (QR-decomposition for maximally
  distinct signals), `"harmonic"` (phase/frequency-shifted sinusoids),
  or `"random"` (random linear combinations of harmonics).

## Value

A matrix of dimension `n_clusters x n_time` with one temporal pattern
per row.
