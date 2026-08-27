# Consensus and slice stitching

## Consensus

``` r

res_cons <- slice_msf(
  toy$vec, toy$mask, target_k_global = -1,
  compactness = 5, nbhd = 8,
  num_runs = 3, consensus = TRUE, stitch_z = FALSE
)
```

## Stitch across z

``` r

res_stitch <- slice_msf(
  toy$vec, toy$mask, target_k_global = -1,
  compactness = 5, nbhd = 8,
  num_runs = 1, consensus = FALSE, stitch_z = TRUE
)
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
val <- validate_cluster4d(res_stitch, toy$vec, toy$mask)
val$summary
```

## Notes & limitations

- With `stitch_z = FALSE`, each axial slice has an independent graph and
  a global exact-K target is unavailable. With `stitch_z = TRUE`, the
  graph adds axial edges and exact `target_k_global` becomes available.
- Consensus (`num_runs > 1, consensus = TRUE`) aggregates multiple runs;
  improves stability at additional runtime cost. includes: in_header:
  \|-
