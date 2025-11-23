# Consensus and slice stitching

## Consensus

``` r
res_cons <- cluster4d(toy$vec, toy$mask, n_clusters = 6, method = "slice_msf",
                      spatial_weight = 0.5, connectivity = 26,
                      num_runs = 3, consensus = TRUE, stitch_z = FALSE)
```

## Stitch across z

``` r
res_stitch <- cluster4d(toy$vec, toy$mask, n_clusters = 6, method = "slice_msf",
                        spatial_weight = 0.5, connectivity = 26,
                        num_runs = 1, consensus = FALSE, stitch_z = TRUE)
```

## Visual comparison

``` r
ok <- exists("res_cons", inherits = TRUE) && exists("res_stitch", inherits = TRUE)
if (ok) {
  par(mfrow = c(1,2))
  plot(res_cons,   view = "axial"); title("consensus (no stitch)")
  plot(res_stitch, view = "axial"); title("stitch_z = TRUE")
  par(mfrow = c(1,1))
} else {
  cat("Consensus/stitch demo unavailable on this platform or dataset size.")
}
#> Consensus/stitch demo unavailable on this platform or dataset size.
```

## Post-merge

``` r
# Example: merge clusters if too small (toy illustration)
val <- validate_cluster4d(res_stitch)
val$summary
```

## Notes & limitations

- Slice-wise clustering operates independently per slice; boundaries can
  be discontinuous across z. `stitch_z = TRUE` merges adjacent slices
  using contact and similarity thresholds.
- Consensus (`num_runs > 1, consensus = TRUE`) aggregates multiple runs;
  improves stability at additional runtime cost. includes: in_header:
  \|-
