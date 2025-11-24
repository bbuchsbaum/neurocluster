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

Here’s a scikit-learn-style “blobs” example: three vertical bands with
distinct time courses plus light noise. It’s tiny, deterministic, and
easy to see what the algorithm is doing.

``` r
toy <- make_block_synthetic(dims = c(12, 12, 1), ntime = 60, noise = 0.1, seed = 7)
toy_res <- cluster4d(
  toy$vec, toy$mask, n_clusters = 3,
  method = "snic", max_iterations = 1
)
plot(toy_res, slice = c(1, 1, 1), view = "axial")
```

![Axial slice showing three contiguous colored
bands.](01-getting-started_files/figure-html/quick-visual-1.png)

Blocky synthetic with three temporal signatures. SNIC recovers the
spatial bands cleanly.

## Export

``` r
writeVol(result$clusvol, "clusters.nii.gz")
```

## See also

- Compare methods: articles/compare-methods.html#run-methods
- Choose parameters: articles/choose-parameters.html#heuristics
- Visualize & export: articles/visualize-export.html#plot
