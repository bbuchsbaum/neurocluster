# Compare multiple cluster4d results

Compares clustering results from different methods or parameters.

## Usage

``` r
compare_cluster4d(
  ...,
  metrics = c("summary", "spatial_coherence", "temporal_coherence")
)
```

## Arguments

- ...:

  cluster4d_result objects to compare

- metrics:

  Comparison metrics to compute. Options:

  - "summary": Basic statistics

  - "spatial_coherence": Spatial compactness measure

  - "temporal_coherence": Feature similarity within clusters

  - "overlap": Dice coefficient between methods (requires exactly 2
    results)

## Value

A comparison data frame
