# Volumetric DCT minimum-spanning-forest clustering

Compresses each voxel time series into a non-DC DCT sketch, uses cosine
dissimilarity on a masked 3-D neighborhood graph, and applies
Felzenszwalb-Huttenlocher segmentation. Multi-run fits use seeded DCT
subspace and mode-weight resampling followed by uniformly weighted
consensus. Nonzero sketches are normalized to unit length. A
zero-information sketch remains zero and therefore has edge
dissimilarity one to every neighbor.

## Usage

``` r
slice_msf(
  vec,
  mask,
  target_k_global = -1,
  r = 12,
  compactness = 5,
  min_size = 80,
  num_runs = 3,
  consensus = TRUE,
  stitch_z = TRUE,
  nbhd = 8,
  gamma = 0,
  k_fuse = NULL,
  min_size_fuse = NULL,
  use_features = FALSE,
  lambda = 0.7,
  z_mult = 0,
  seed = 1L,
  ensemble_fraction = 0.8
)
```

## Arguments

- vec:

  A `NeuroVec` or `SparseNeuroVec` containing the time series.

- mask:

  A `NeuroVol`; exactly finite values greater than zero are included.

- target_k_global:

  `-1` for the natural segmentation or a positive exact cluster count.
  Exact K uses the shared adjacency-preserving split/merge engine.

- r:

  Requested number of non-DC DCT modes. It is capped at `T - 1`.

- compactness:

  Finite value in `[0, 10]` controlling FH scale. Larger values use a
  smaller scale and generally produce finer components.

- min_size:

  Positive minimum component size for each run.

- num_runs:

  Positive number of runs. Values above one require consensus.

- consensus:

  Whether to fuse a multi-run ensemble.

- stitch_z:

  Whether graph edges may cross axial slices. When `FALSE`, slices are
  independent and a global exact-K target is unavailable.

- nbhd:

  `4`, `6`, or `8`. Values `4` and `6` select axial neighbors; `8`
  selects the full diagonal neighborhood.

- gamma:

  Reserved compatibility parameter. It must be zero; positive
  reliability weighting is rejected because it can collapse
  non-identical feature vectors to zero-distance edges.

- k_fuse:

  Positive FH scale for consensus; defaults to the run scale.

- min_size_fuse:

  Positive minimum consensus component size; defaults to `min_size`.

- use_features:

  Include sketch similarity as well as label agreement in multi-run
  consensus.

- lambda:

  Consensus mixture in `[0, 1]`; one uses label agreement only, zero
  uses feature similarity only. Active only when `use_features = TRUE`.

- z_mult:

  DCT-sketch smoothing fraction in `[0, 1]` across axial slices.

- seed:

  Integer seed used to construct multi-run DCT subspaces.

- ensemble_fraction:

  Fraction in `(0, 1]` controlling the candidate DCT frequency pool for
  a multi-run ensemble.

## Value

A `slice_msf_cluster_result` and `cluster4d_result`. `cluster` is in
mask order, `centers` is K by T, and `coord_centers` is K by 3.
Multi-run results also contain native `runs`. `metadata$exact_k_repair`
records the natural and requested K, pre/post sizes, whether and how
repair ran, connectivity, minimum size, and effective gamma, z
smoothing, seed, run count, and consensus controls for both natural and
exact-K calls.

## Experimental status

`slice_msf()` is available for explicit evaluation but is not a
recommended production method. A previous reliability-weighting path
collapsed distinct feature vectors to zero-distance edges, and the
former exact-K repair could create singleton-dominated partitions. Those
paths are now rejected or replaced and covered by regression tests, but
release-level recertification is still pending. See
[`vignette("slice-msf")`](https://bbuchsbaum.github.io/neurocluster/articles/slice-msf.md)
for the evaluation contract.

## References

Felzenszwalb, P. F., & Huttenlocher, D. P. (2004). Efficient graph-based
image segmentation. IJCV, 59(2), 167-181.
