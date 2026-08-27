# Compare multiple cluster4d results

Compares clustering results over exactly the same included mask voxels
and physical coordinate system. Ambiguous center-based scores are not
used.

## Usage

``` r
compare_cluster4d(
  ...,
  metrics = c("summary", "spatial_dispersion"),
  feature_mat = NULL
)
```

## Arguments

- ...:

  cluster4d_result objects to compare

- metrics:

  Character vector selecting explicit estimands:

  - `"summary"`: cluster-count and cluster-size summaries.

  - `"spatial_dispersion"`: root mean squared physical distance in
    millimetres from voxels to their final cluster centroid (lower is
    more compact).

  - `"temporal_coherence"`: mean Pearson correlation across all
    unordered within-cluster voxel pairs, requiring `feature_mat`.

  - `"partition_agreement"`: adjusted Rand index, variation of
    information in bits, and pairwise Dice for exactly two results.

- feature_mat:

  A shared finite numeric matrix with voxels in rows and at least two
  time/features in columns. Required only for `"temporal_coherence"`.
  Rows must be non-constant.

## Value

A data frame with one row per result and explicitly named metrics.

## Details

Every result must carry canonical mask and affine provenance. The exact
included voxel indices and physical geometry are checked before any
cross-result comparison. Partition metrics use only observed contingency
cells and never construct an N-by-N co-membership matrix. Pairwise
agreement columns are repeated on the two result rows so the return
value remains a data frame.
