# Fuse Slice-MSF runs on a masked graph

Builds neighbor-edge co-association scores from at least two run results
and optionally mixes those scores with DCT-sketch similarity.

## Usage

``` r
slice_msf_consensus(
  run_results,
  mask,
  nbhd = 8,
  k_fuse = 0.3,
  min_size_fuse = 80,
  use_features = FALSE,
  lambda = 0.7,
  stitch_z = TRUE,
  gamma = 0
)
```

## Arguments

- run_results:

  At least two results from
  [`slice_msf_single()`](https://bbuchsbaum.github.io/neurocluster/reference/slice_msf_single.md).

- mask:

  The original `NeuroVol` mask. Outside-mask labels are discarded.

- nbhd:

  `4`, `6`, or `8`; see
  [`slice_msf()`](https://bbuchsbaum.github.io/neurocluster/reference/slice_msf.md).

- k_fuse:

  Positive FH scale for the consensus graph.

- min_size_fuse:

  Positive minimum consensus component size.

- use_features:

  Whether DCT-sketch similarity contributes to edges.

- lambda:

  Mixture in `[0, 1]`; active when `use_features = TRUE`.

- stitch_z:

  Whether consensus edges may cross axial slices.

- gamma:

  Reserved compatibility parameter. It must be zero.

## Value

A list with full-volume consensus `labels` and complete `params`.
