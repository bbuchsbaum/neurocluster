# Spatially constrained clustering

## Intuition

Spatially constrained clustering balances data similarity (time series)
with spatial proximity to produce contiguous parcels.

## Alpha vs. compactness

- Supervoxels: `alpha` controls the relative weight of feature
  similarity; in
  [`cluster4d_supervoxels()`](../reference/cluster4d_supervoxels.md),
  `spatial_weight` is mapped as `alpha = 1 - spatial_weight`.
- SNIC/SLIC: use a `compactness` term in their native interfaces; via
  [`cluster4d()`](../reference/cluster4d.md), set `spatial_weight` which
  the wrappers map to the appropriate compactness scale.

## Connectivity

- Use 6 for face-adjacency; 26 adds edges/corners. In
  [`cluster4d()`](../reference/cluster4d.md), `connectivity` is passed
  through to the chosen method (where applicable).

## Artifacts

- Slice-wise methods can exhibit z‑plane seams; stitch or increase
  spatial weighting. includes: in_header: \|-
