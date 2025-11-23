# Validate and compare clusterings

## Checks

``` r
val <- validate_cluster4d(res1)
val
```

## Sizes

``` r
table(res1$cluster)
table(res2$cluster)
```

## Compare

``` r
compare_cluster4d(res1, res2)
```

## Report

- Flag tiny clusters for potential merging.
- Note between-cluster correlations for interpretability.
