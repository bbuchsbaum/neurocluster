# Check for Spatial Ordering in Cluster IDs

Diagnostic function to detect whether cluster IDs are spatially ordered
along any axis. Useful for identifying clustering algorithms that may
produce ordered IDs.

## Usage

``` r
check_cluster_ordering(cluster_result)
```

## Arguments

- cluster_result:

  A clustering result object with `coord_centers` component

## Value

A list with correlation statistics:

- `cor_x`: Spearman correlation between cluster ID and x-coordinate

- `cor_y`: Spearman correlation between cluster ID and y-coordinate

- `cor_z`: Spearman correlation between cluster ID and z-coordinate

- `max_abs_cor`: Maximum absolute correlation across all axes

- `is_ordered`: Logical indicating if any axis shows strong ordering
  (\|r\| \> 0.3)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- cluster4d(vec, mask, K = 100, method = "flash3d")

# Check for ordering
ordering <- check_cluster_ordering(result)
print(ordering)

# If strongly ordered, randomize
if (ordering$is_ordered) {
  result <- randomize_cluster_ids(result)
}
} # }
```
