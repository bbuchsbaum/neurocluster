# Fast 4D SLIC supervoxels (mask-aware, gradient relocation, preserve-K)

Cluster a 4D `NeuroVec` (x,y,z,time) into compact, spatially contiguous
3D supervoxels using an enhanced SLIC-style algorithm with mask-aware
seeding, gradient-based seed relocation, and exact K preservation.

## Usage

``` r
slic4d_supervoxels(
  bvec,
  mask,
  K,
  compactness = 10,
  max_iter = 10,
  n_threads = 0,
  step_mm = NULL,
  n_components = 0,
  feature_norm = c("zscale", "l2", "none"),
  seed_method = c("mask_poisson", "mask_grid", "grid", "farthest"),
  seed_relocate = c("none", "correlation", "intensity", "spatial"),
  seed_relocate_radius = 1L,
  connectivity = c(26L, 6L),
  strict_connectivity = TRUE,
  enforce_connectivity = TRUE,
  preserve_k = FALSE,
  topup_iters = 2L,
  min_size = 0L,
  verbose = FALSE
)
```

## Arguments

- bvec:

  A `NeuroVec` with dims (X, Y, Z, T).

- mask:

  A 3D `NeuroVol` (or logical array) indicating voxels to include.

- K:

  Target number of supervoxels.

- compactness:

  Spatial vs feature tradeoff (like SLIC 'm'). Larger = more compact.

- max_iter:

  Maximum iterations (default 10).

- n_threads:

  Number of CPU threads to use (0 = auto).

- step_mm:

  Optional approximate spacing between seeds in millimeters; if NULL,
  computed as cubic root of bounding-box volume / K.

- n_components:

  If \> 0, random-project each voxel's time series to this dimension for
  speed (Johnson-Lindenstrauss style). 0 = use full time series.

- feature_norm:

  One of "zscale", "l2", "none".

- seed_method:

  One of "mask_poisson" (Poisson disk in mask), "mask_grid" (grid in
  mask), "grid" (regular grid), "farthest" (farthest point sampling).

- seed_relocate:

  One of "correlation" (correlation gradient), "intensity" (mean
  intensity gradient), "spatial" (spatial gradient using
  neighborweights), "none" (no relocation).

- seed_relocate_radius:

  Search radius in voxels for gradient-based seed relocation (default
  1).

- connectivity:

  Neighborhood connectivity: 6 (face neighbors) or 26 (all neighbors).

- strict_connectivity:

  Enforce exactly one connected component per label (default TRUE).

- enforce_connectivity:

  Alias for strict_connectivity (backward compatibility).

- preserve_k:

  Ensure exactly K non-empty labels by splitting largest clusters
  (default FALSE).

- topup_iters:

  Number of refinement iterations after splitting for preserve_k
  (default 2).

- min_size:

  Minimum component size (voxels) to keep before relabel (default 0).

- verbose:

  Logical.

## Value

A list of class `cluster_result` with elements:

- `clusvol`: `ClusteredNeuroVol` with the final labels

- `cluster`: integer vector of length = \#masked voxels

- `centers`: matrix (K x d_feat) center features

- `coord_centers`: matrix (K x 3) spatial centers in mm

## Note

Consider using [`cluster4d`](cluster4d.md) with `method = "slic"` for a
standardized interface across all clustering methods.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(neuroim2)
  # Basic usage
  res <- slic4d_supervoxels(bvec, mask, K = 1000, compactness = 15)
  
  # With mask-aware seeding and gradient relocation
  res <- slic4d_supervoxels(bvec, mask, K = 1000, 
                           seed_method = "mask_poisson",
                           seed_relocate = "correlation",
                           preserve_k = TRUE)
} # }
```
