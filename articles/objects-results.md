# Objects and results

## cluster4d_result

Canonical fields: `labels` (with identical backward-compatible alias
`cluster`), `clusvol`, `centers`, `coord_centers`, `actual_k` (with
alias `n_clusters`), `label_ids`, `method`, `parameters`, typed
`provenance`, and method-specific `metadata`.

## Print/plot/summary

- `print(x)`: method, counts, cluster size stats
- `summary(x)`: parameters, cluster size distribution, ranges
- `plot(x)`: axial/sagittal/coronal slices

## neuroim2 objects

- `NeuroVol`, `NeuroVec`, `ClusteredNeuroVol`, helpers (`series`,
  `index_to_coord`, `spacing`)
