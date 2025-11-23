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

To make parameter effects concrete, we call
[`generate_synthetic_volume()`](../reference/generate_synthetic_volume.md)
to create a tiny 4D volume where the true clusters form a clear 2D 2×3
grid on a single axial slice. Voxels within each grid cell share a
similar time course (high correlation), so a well‑tuned clustering
should recover these blocks. This runs quickly and deterministically.

## Effect of K (number of clusters)

Larger K produces finer partitions. The toy volume below is clustered
with the same method and weights, varying only K.

``` r
par(mfrow = c(1, 3))
res_k3 <- cluster4d(toy$vec, toy$mask, n_clusters = 3)
plot(res_k3, slice = c(1, 1, 1), view = "axial"); title("K = 3")
```

![Three panels showing axial slices clustered with K=3, K=6, and
K=9.](10-choose-parameters_files/figure-html/fig-k-1.png)

Effect of K on a toy axial view (default method).

``` r
res_k6 <- cluster4d(toy$vec, toy$mask, n_clusters = 6)
plot(res_k6, slice = c(1, 1, 1), view = "axial"); title("K = 6")
```

![Three panels showing axial slices clustered with K=3, K=6, and
K=9.](10-choose-parameters_files/figure-html/fig-k-2.png)

Effect of K on a toy axial view (default method).

``` r
res_k9 <- cluster4d(toy$vec, toy$mask, n_clusters = 9)
plot(res_k9, slice = c(1, 1, 1), view = "axial"); title("K = 9")
```

![Three panels showing axial slices clustered with K=3, K=6, and
K=9.](10-choose-parameters_files/figure-html/fig-k-3.png)

Effect of K on a toy axial view (default method).

``` r
par(mfrow = c(1,1))
```

## Effect of spatial_weight

Higher `spatial_weight` emphasizes spatial compactness relative to
feature similarity.

``` r
par(mfrow = c(1, 3))
res_sw02 <- cluster4d(toy$vec, toy$mask, n_clusters = 6, spatial_weight = 0.2)
plot(res_sw02, slice = c(1, 1, 1), view = "axial"); title("spatial_weight = 0.2")
```

![Three panels with spatial_weight 0.2, 0.5,
0.8.](10-choose-parameters_files/figure-html/fig-sw-1.png)

Effect of spatial_weight on a toy axial view (K=6).

``` r
res_sw05 <- cluster4d(toy$vec, toy$mask, n_clusters = 6, spatial_weight = 0.5)
plot(res_sw05, slice = c(1, 1, 1), view = "axial"); title("spatial_weight = 0.5")
```

![Three panels with spatial_weight 0.2, 0.5,
0.8.](10-choose-parameters_files/figure-html/fig-sw-2.png)

Effect of spatial_weight on a toy axial view (K=6).

``` r
res_sw08 <- cluster4d(toy$vec, toy$mask, n_clusters = 6, spatial_weight = 0.8)
plot(res_sw08, slice = c(1, 1, 1), view = "axial"); title("spatial_weight = 0.8")
```

![Three panels with spatial_weight 0.2, 0.5,
0.8.](10-choose-parameters_files/figure-html/fig-sw-3.png)

Effect of spatial_weight on a toy axial view (K=6).

``` r
par(mfrow = c(1,1))
```

## Effect of connectivity

Connectivity controls the local neighborhood (6 = faces; 26 =
faces/edges/corners). Higher connectivity typically yields smoother
boundaries.

``` r
par(mfrow = c(1, 2))
res_c6  <- cluster4d(toy$vec, toy$mask, n_clusters = 6, spatial_weight = 0.5, connectivity = 6)
plot(res_c6, slice = c(1, 1, 1), view = "axial"); title("connectivity = 6")
```

![Two panels: connectivity 6 and connectivity
26.](10-choose-parameters_files/figure-html/fig-connectivity-1.png)

Effect of connectivity on a toy axial view (K=6, spatial_weight=0.5).

``` r
res_c26 <- cluster4d(toy$vec, toy$mask, n_clusters = 6, spatial_weight = 0.5, connectivity = 26)
plot(res_c26, slice = c(1, 1, 1), view = "axial"); title("connectivity = 26")
```

![Two panels: connectivity 6 and connectivity
26.](10-choose-parameters_files/figure-html/fig-connectivity-2.png)

Effect of connectivity on a toy axial view (K=6, spatial_weight=0.5).

``` r
par(mfrow = c(1,1))
```

includes: in_header: \|-
