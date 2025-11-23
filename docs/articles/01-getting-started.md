# Getting started with neurocluster

## Goal

Cluster a small 4D volume (3D space × time) into a few spatially
coherent parcels, inspect the result, and export a NIfTI file.

## TL;DR

``` r
# Create a small but structured synthetic dataset
syn <- generate_synthetic_volume(
  scenario = "gaussian_blobs",
  dims = c(10, 10, 6),
  n_clusters = 4,
  n_time = 24,
  seed = 42
)

# Minimal clustering: SNIC is fast on small data
result <- cluster4d(syn$vec, syn$mask, n_clusters = 4, method = "snic", max_iterations = 1)
```

## Setup

``` r
syn <- generate_synthetic_volume(
  scenario = "gaussian_blobs",
  dims = c(10, 10, 6),
  n_clusters = 4,
  n_time = 24,
  seed = 42
)
brain_mask <- syn$mask
vec  <- syn$vec
```

## Peek at the ground truth

``` r
# Display the ground truth cluster labels
truth_array <- array(syn$truth, dim = syn$dims)
image(truth_array[,,3], col = rainbow(syn$n_clusters),
      main = "Ground Truth: Slice 3", axes = FALSE)
```

![Axial slice showing four distinct
clusters.](01-getting-started_files/figure-html/truth-peek-1.png)

Ground-truth Gaussian blobs (axial slice).

## Walkthrough

``` r
result <- cluster4d(vec, mask, n_clusters = 4)
print(result)
summary(result)
plot(result)
```

## Quick visual (toy blocks)

A tiny, seeded toy illustrates a clear 2D grid of contiguous clusters
(K=6) on a single axial slice. Voxels within each grid cell share a
similar time course (high correlation), so clustering should recover the
underlying block structure. Colors are arbitrary cluster IDs.

``` r
toy <- generate_synthetic_volume(
  scenario = "checkerboard",
  dims = c(12, 12, 1),
  n_clusters = 6,
  n_time = 16,
  seed = 123
)
toy_res <- cluster4d(
  toy$vec, toy$mask, n_clusters = 6,
  method = "slic", spatial_weight = 0.5, connectivity = 26,
  max_iterations = 10
)
plot(toy_res, slice = c(1, 1, 1), view = "axial")
```

![Axial slice showing several contiguous colored
regions.](01-getting-started_files/figure-html/quick-visual-1.png)

Toy axial slice clustered with K=6; contiguous blocks roughly follow a
2×3 grid (colors are arbitrary cluster IDs).

## Export

``` r
writeVol(result$clusvol, "clusters.nii.gz")
```

## See also

- Compare methods: articles/compare-methods.html#run-methods
- Choose parameters: articles/choose-parameters.html#heuristics
- Visualize & export: articles/visualize-export.html#plot
