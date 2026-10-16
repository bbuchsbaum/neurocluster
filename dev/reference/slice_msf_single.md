# Run one Slice-MSF segmentation

Low-level diagnostic interface returning native full-volume labels, a
signed adjacent-pair temporal-smoothness diagnostic, and an `r_used` by
voxel DCT sketch. Temporal smoothness does not affect segmentation.

## Usage

``` r
slice_msf_single(
  vec,
  mask,
  r = 12,
  k = 0.32,
  min_size = 80,
  nbhd = 8,
  stitch_z = TRUE,
  gamma = 0,
  z_mult = 0,
  dct_frequencies = NULL,
  dct_weights = NULL
)
```

## Arguments

- vec:

  A `NeuroVec` or `SparseNeuroVec` containing the time series.

- mask:

  A `NeuroVol`; exactly finite values greater than zero are included.

- r:

  Requested number of non-DC DCT modes. It is capped at `T - 1`.

- k:

  Positive FH segmentation scale.

- min_size:

  Positive minimum component size for each run.

- nbhd:

  `4`, `6`, or `8`. Values `4` and `6` select axial neighbors; `8`
  selects the full diagonal neighborhood.

- stitch_z:

  Whether graph edges may cross axial slices. When `FALSE`, slices are
  independent and a global exact-K target is unavailable.

- gamma:

  Reserved compatibility parameter. It must be zero; positive
  reliability weighting is rejected because it can collapse
  non-identical feature vectors to zero-distance edges.

- z_mult:

  DCT-sketch smoothing fraction in `[0, 1]` across axial slices.

- dct_frequencies:

  Optional unique non-DC DCT frequencies in `[1, T - 1]`.

- dct_weights:

  Optional positive weights paired with `dct_frequencies`.

## Value

A list with full-volume `labels`, `temporal_smoothness`, `sketch`, and
`params`.
