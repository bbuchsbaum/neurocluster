# Supervoxel Clustering in Time

Cluster feature matrix (rows = time points) in a "supervoxel" style but
over temporal dimension.

## Usage

``` r
supervoxel_cluster_time(
  feature_mat,
  K = min(nrow(feature_mat), 100),
  sigma1 = 1,
  sigma2 = 3,
  iterations = 50,
  TR = 2,
  filter = list(lp = 0, hp = 0),
  use_medoid = FALSE,
  nreps = 5
)
```

## Arguments

- feature_mat:

  A matrix (nrows = time points, ncols = features) or vice versa.

- K:

  Number of clusters.

- sigma1:

  Heat kernel bandwidth for feature similarity (data vectors).

- sigma2:

  Heat kernel bandwidth for spatial similarity (coordinate vectors).

- iterations:

  Maximum number of cluster iterations.

- TR:

  Repetition time (seconds).

- filter:

  List specifying optional frequency filters, e.g.,
  `list(lp=0.1, hp=0)`.

- use_medoid:

  Whether to use medoids for cluster centers.

- nreps:

  Number of repeated initializations.

## Value

A list of cluster results (one per repetition), each of which has the
same structure as `supervoxel_cluster_fit()`.

## Examples

``` r
feature_mat <- matrix(rnorm(100 * 10), 100, 10)
library(future)
#> Warning: package ‘future’ was built under R version 4.5.2
plan(multicore)
cres <- supervoxel_cluster_time(t(feature_mat), K=20)
#> Error in (function (.x, .f, ..., .progress = FALSE) {    map_("list", .x, .f, ..., .progress = .progress)})(.x = 1L, .f = function (...) {    NULL    NULL    ...furrr_out <- ...furrr_fn(...)    ...furrr_out}): ℹ In index: 1.
#> Caused by error in `sample.int()`:
#> ! cannot take a sample larger than the population when 'replace = FALSE'
```
