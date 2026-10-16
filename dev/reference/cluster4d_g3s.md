# Gradient-guided geodesic clustering on a masked grid

Compresses voxel features, selects low-gradient seeds, and assigns
voxels by deterministic multi-source Dijkstra propagation. Edges are
exactly the requested masked-grid neighbors; excluded gaps are never
bridged. The edge objective is
`alpha * feature_distance + (1 - alpha) * physical_distance / compactness`.

## Usage

``` r
cluster4d_g3s(
  vec,
  mask,
  K = 100,
  n_components = 15,
  variance_threshold = 0.95,
  alpha = 0.5,
  compactness = NULL,
  max_refinement_iter = 3,
  verbose = FALSE,
  use_irlba = TRUE,
  use_rsvd = TRUE,
  connectivity = 26
)
```

## Arguments

- vec:

  A `NeuroVec` or `SparseNeuroVec` containing the input features.

- mask:

  A `NeuroVol`; exactly finite values greater than zero are included.

- K:

  Positive requested cluster count.

- n_components:

  Requested SVD feature count.

- variance_threshold:

  Minimum retained variance in `[0, 1]`.

- alpha:

  Feature-distance weight in `[0, 1]`; zero is spatial-only and one is
  feature-only.

- compactness:

  Positive finite physical length scale. Smaller values impose a
  stronger spatial penalty. The default uses masked voxel measure, K,
  and the effective dimension, so single-slice masks remain
  well-defined.

- max_refinement_iter:

  Non-negative number of boundary passes. Refinement uses the same
  feature/spatial mixture and final labels are repaired through the
  shared adjacency-preserving exact-K engine.

- verbose:

  Whether to report phase progress.

- use_irlba:

  Whether large SVDs may use `irlba`.

- use_rsvd:

  Whether randomized SVD may be used when available.

- connectivity:

  Exact grid connectivity: 6, 18, or 26.

## Value

A typed `g3s_result` and `cluster4d_result` with K by T centers,
physical coordinate centers, graph provenance, and physical-scale
metadata.

## References

Dijkstra, E. W. (1959). A note on two problems in connexion with graphs.
Numerische Mathematik, 1, 269-271.
