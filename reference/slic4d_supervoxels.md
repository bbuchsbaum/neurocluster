# Fast 4D SLIC supervoxels (mask-aware, gradient relocation, preserve-K)

Cluster a 4D `NeuroVec` (x,y,z,time) into compact, spatially contiguous
3D supervoxels using an enhanced SLIC-style algorithm with mask-aware
seeding, gradient-based seed relocation, and topology-aware K
preservation.

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
  verbose = FALSE,
  .input_contract = NULL
)
```

## Arguments

- bvec:

  A `NeuroVec` with dims (X, Y, Z, T). Can also be a 3D
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.html)
  for structural image segmentation, which will be automatically
  converted to a single-timepoint NeuroVec internally.

- mask:

  A 3D `NeuroVol` indicating voxels to include. Finite values greater
  than zero are included; zero and negative values are excluded.

- K:

  Finite integer target number of supervoxels, between one and the
  number of included voxels.

- compactness:

  Finite non-negative spatial vs feature tradeoff (like SLIC 'm').
  Larger values produce more compact clusters.

- max_iter:

  Positive finite integer iteration limit (default 10).

- n_threads:

  Non-negative finite integer number of CPU threads. Zero requests
  RcppParallel's automatic thread count, one forces serial assignment,
  and values above one request that many threads. Assignment is
  intentionally serial for masks with fewer than 2,000 included voxels.

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
  intensity gradient), "spatial" (spatial gradient using adjoin), "none"
  (no relocation).

- seed_relocate_radius:

  Search radius in voxels for gradient-based seed relocation (default
  1).

- connectivity:

  Neighborhood connectivity: 6 (face neighbors) or 26 (all neighbors).

- strict_connectivity:

  Enforce exactly one connected component per label (default TRUE).

- enforce_connectivity:

  Deprecated alias for strict connectivity. When omitted, it follows
  `strict_connectivity`; if either supplied value is TRUE, connectivity
  is enforced.

- preserve_k:

  Return exactly K non-empty labels when compatible with strict
  connectivity (default FALSE). See Details.

- topup_iters:

  Non-negative number of native refinement iterations before the final
  connectivity and exact-K repair (default 2).

- min_size:

  Minimum component size (voxels) to keep before relabel (default 0).

- verbose:

  Logical.

- .input_contract:

  Internal prevalidated input receipt used by
  [`cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md);
  direct callers should leave this as `NULL`.

## Value

A `cluster4d_result` (also inheriting from `cluster_result`) with
elements:

- `clusvol`: `ClusteredNeuroVol` with the final labels

- `cluster`: integer vector of length = \#masked voxels

- `centers`: matrix (actual_k x d_feat) center features

- `coord_centers`: matrix (actual_k x 3) spatial centers in mm

## Details

Connectivity is the final partition constraint. With
`strict_connectivity = TRUE`, every returned label has exactly one
connected component under the requested neighborhood. If
`preserve_k = TRUE`, an adjacency-preserving repair reaches exactly K
whenever K is at least the number of connected mask components. When K
is below that topological minimum, strict connectivity takes precedence,
the minimum feasible number of labels is returned, and a warning is
emitted. Returned feature and coordinate centers are recomputed from the
final public labels in the original feature space and physical
coordinate system.

## Note

Consider using
[`cluster4d`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md)
with `method = "slic"` for a standardized interface across all
clustering methods.

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
