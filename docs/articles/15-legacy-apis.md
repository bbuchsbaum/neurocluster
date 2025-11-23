# Legacy APIs and backward compatibility

## Map old to new

- `supervoxels(vec, mask, K, alpha, ...)` →
  `cluster4d(vec, mask, n_clusters=K, method="supervoxels", spatial_weight=1-alpha, ...)`
- `snic(vec, mask, K, compactness, ...)` →
  `cluster4d(..., method="snic", spatial_weight=compactness/10, ...)`
- `slic4d_supervoxels(...)` → `cluster4d(..., method="slic", ...)`
- `slice_msf(...)` → `cluster4d(..., method="slice_msf", ...)`
- `supervoxels_flash3d(...)` → `cluster4d(..., method="flash3d", ...)`

## Examples

``` r
# Old
# res <- supervoxels(vec, mask, K = 100, alpha = 0.3)

# New
# res <- cluster4d(vec, mask, n_clusters = 100, method = "supervoxels", spatial_weight = 0.7)
```

## Deprecations

- Prefer [`cluster4d()`](../reference/cluster4d.md) for consistency and
  shared validation. includes: in_header: \|-
