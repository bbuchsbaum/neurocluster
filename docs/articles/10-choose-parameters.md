# Choose parameters for your data

## Heuristics

The package’s helper
[`suggest_cluster4d_params()`](../reference/suggest_cluster4d_params.md)
derives a baseline K as `n_voxels / 250` (see source in
`cluster4d_common.R`) and adjusts by priority. This is a starting point;
actual K depends on analysis goals and mask size.

## Suggestions

``` r
n_vox <- 50000
n_time <- 200
params <- suggest_cluster4d_params(n_vox, n_time, priority = "balanced")
params$recommended_method
params$n_clusters
```

## Tradeoffs

- Quality: increase iterations; methods with iterative refinement (e.g.,
  [`supervoxels()`](../reference/supervoxels.md),
  [`cluster4d_slic()`](../reference/cluster4d_slic.md)); consider higher
  K.
- Speed: `slice_msf` (slice-wise) or `flash3d` (hash/DCT) with fewer
  iterations and lower K.
- Memory: [`snic()`](../reference/snic.md) (single pass), `slice_msf`
  with fewer coefficients.

## Recipes

- Whole-brain (2–3 mm): start from `n_voxels/250`, `spatial_weight`
  around 0.4–0.6, connectivity 26.
- ROI: smaller K (e.g., 10–100), lower `spatial_weight` (0.2–0.4),
  connectivity 6.

## Toy data for illustrations

To make parameter effects concrete, we use
[`make_block_synthetic()`](../reference/make_block_synthetic.md)—three
vertical bands with distinct time courses, light noise, and a single
slice. It’s fast, deterministic, and spatially local.

## Effect of K (number of clusters)

Larger K produces finer partitions. Below we fix the method (`snic`) and
compactness, varying only K.

``` r
par(mfrow = c(1, 3))
res_k2 <- snic(toy$vec, toy$mask, K = 2, compactness = 3)
plot(res_k2, slice = c(1, 1, 1), view = "axial"); title("K = 2")
```

![Three panels showing axial slices clustered with K=2, K=3, and
K=5.](10-choose-parameters_files/figure-html/fig-k-1.png)

Effect of K on a block synthetic (snic).

``` r
res_k3 <- snic(toy$vec, toy$mask, K = 3, compactness = 3)
plot(res_k3, slice = c(1, 1, 1), view = "axial"); title("K = 3")
```

![Three panels showing axial slices clustered with K=2, K=3, and
K=5.](10-choose-parameters_files/figure-html/fig-k-2.png)

Effect of K on a block synthetic (snic).

``` r
res_k5 <- snic(toy$vec, toy$mask, K = 5, compactness = 3)
plot(res_k5, slice = c(1, 1, 1), view = "axial"); title("K = 5")
```

![Three panels showing axial slices clustered with K=2, K=3, and
K=5.](10-choose-parameters_files/figure-html/fig-k-3.png)

Effect of K on a block synthetic (snic).

``` r
par(mfrow = c(1,1))
```

## Effect of compactness (SNIC)

Higher `compactness` makes clusters more spatially tight; lower values
follow feature patterns more closely.

``` r
par(mfrow = c(1, 3))
res_c2 <- snic(toy$vec, toy$mask, K = 3, compactness = 2)
plot(res_c2, slice = c(1, 1, 1), view = "axial"); title("compactness = 2")
```

![Three panels with compactness 2, 4,
6.](10-choose-parameters_files/figure-html/fig-compactness-1.png)

Effect of compactness on the block synthetic (K=3, snic).

``` r
res_c4 <- snic(toy$vec, toy$mask, K = 3, compactness = 4)
plot(res_c4, slice = c(1, 1, 1), view = "axial"); title("compactness = 4")
```

![Three panels with compactness 2, 4,
6.](10-choose-parameters_files/figure-html/fig-compactness-2.png)

Effect of compactness on the block synthetic (K=3, snic).

``` r
res_c6 <- snic(toy$vec, toy$mask, K = 3, compactness = 6)
plot(res_c6, slice = c(1, 1, 1), view = "axial"); title("compactness = 6")
```

![Three panels with compactness 2, 4,
6.](10-choose-parameters_files/figure-html/fig-compactness-3.png)

Effect of compactness on the block synthetic (K=3, snic).

``` r
par(mfrow = c(1,1))
```

includes: in_header: \|-
