# Read, validate, and save clustering results

A clustering is useful only if its labels can be tied back to the mask,
feature space, and physical geometry that produced them.
`cluster4d_result` keeps those pieces together. Use its public fields
and validate the object against the original inputs before
interpretation or export.

## What does a valid result look like?

``` r

fit <- cluster4d(
  syn$vec, syn$mask, n_clusters = 4,
  method = "slic", preserve_k = TRUE,
  max_iterations = 6, parallel = FALSE
)
receipt <- validate_cluster4d(fit, syn$vec, syn$mask)
receipt$summary
#> $n_clusters
#> [1] 4
#> 
#> $n_voxels
#> [1] 576
#> 
#> $cluster_size_range
#> [1] 132 155
#> 
#> $has_centers
#> [1] TRUE
#> 
#> $has_spatial_centers
#> [1] TRUE
```

The validation call independently checks the label schema, centers,
physical centroids, mask ordering, and geometry. Calling it without
`vec` and `mask` cannot certify all of those relationships.

## How are labels related to centers?

| Field | Meaning |
|----|----|
| `labels` | Contiguous positive IDs in included-mask voxel order |
| `cluster` | An identical compatibility alias of `labels` |
| `label_ids` | The label represented by each row of both center matrices |
| `centers` | Means in the original feature space; rows are labels |
| `coord_centers` | Mean physical coordinates in millimetres; rows are labels |
| `actual_k` | Number of final labels; may differ from a requested target when topology makes the target infeasible |
| `clusvol` | A [`neuroim2::ClusteredNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.html) carrying labels and mask geometry |

Never index a center matrix by an undocumented label ordering. The
explicit `label_ids` mapping is the contract.

## What should you notice in a slice?

![Axial grid divided into four colored regions; connectivity through
other slices cannot be inferred from this panel
alone.](objects-results_files/figure-html/result-plot-1.png)

One axial slice of the validated SLIC fit. Colours identify labels, not
an ordering or effect magnitude. A slice is a visual screen, not proof
of 3-D connectivity.

**What to notice.** The figure screens label coverage and gross spatial
shape, but it cannot establish that each label has one connected
component in 3-D. Colour numbers are arbitrary; compare partitions with
defined metrics, not by matching hues.

| label | voxels | components_26 |
|------:|-------:|--------------:|
|     1 |    134 |             1 |
|     2 |    155 |             1 |
|     3 |    155 |             1 |
|     4 |    132 |             1 |

Observed 26-neighbor component count for every final label. {.table}

The fitted result has one 26-neighbor component per label on this
fixture. That statement comes from an explicit whole-volume traversal
and hidden assertion, not from the middle slice. Use the connectivity
definition declared by the method when reproducing this diagnostic.

## What belongs in a saved receipt?

Save the fitted object together with the original mask identity and
analysis metadata. At minimum retain `method`, `parameters`, `actual_k`,
`provenance`, and the result of
[`validate_cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/validate_cluster4d.md).
For comparison and export workflows, continue to [Validate and compare
clusterings](https://bbuchsbaum.github.io/neurocluster/articles/validate-compare.md)
and [Visualize and export
results](https://bbuchsbaum.github.io/neurocluster/articles/visualize-export.md).
