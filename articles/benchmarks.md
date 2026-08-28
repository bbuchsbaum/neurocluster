# Benchmark with a reproducible receipt

A benchmark is evidence only for its measured workload. This article
does not read a checked-in `results.csv` or present old machine timings
as current. It runs a deliberately small comparison during rendering and
shows the receipt needed to reproduce it. Use the pattern—not the tiny
timings—to design a study for your data.

## What exactly is measured here?

The estimand is elapsed wall-clock seconds for a complete
[`cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md)
call on one generated 400-voxel, 18-time-point volume, after one warm-up
call. We keep three repetitions for each method and require every fit to
be valid and to produce the requested K. These are smoke-scale timings,
not whole-brain performance claims.

``` r

methods <- c("snic", "slic")
benchmark_order <- c("snic", "slic", "slic", "snic", "snic", "slic")
repetitions <- as.integer(ave(
  seq_along(benchmark_order), benchmark_order, FUN = seq_along
))
data.frame(run = seq_along(benchmark_order), method = benchmark_order,
           repetition = repetitions)
#>   run method repetition
#> 1   1   snic          1
#> 2   2   slic          1
#> 3   3   slic          2
#> 4   4   snic          2
#> 5   5   snic          3
#> 6   6   slic          3

timings <- do.call(rbind, lapply(seq_along(benchmark_order), function(run) {
  method <- benchmark_order[run]
  elapsed <- system.time(fit <- fit_one(method))[["elapsed"]]
  validation <- validate_cluster4d(fit, syn$vec, syn$mask)
  data.frame(
    run = run, method = method, repetition = repetitions[run],
    elapsed_seconds = unname(elapsed),
    actual_k = fit$actual_k, valid = validation$valid
  )
}))
timings
#>   run method repetition elapsed_seconds actual_k valid
#> 1   1   snic          1           0.045        4  TRUE
#> 2   2   slic          1           0.008        4  TRUE
#> 3   3   slic          2           0.008        4  TRUE
#> 4   4   snic          2           0.047        4  TRUE
#> 5   5   snic          3           0.047        4  TRUE
#> 6   6   slic          3           0.007        4  TRUE
```

![Jittered timing points and median diamonds for SNIC and serial SLIC on
a 400-voxel
example.](benchmarks_files/figure-html/benchmark-figure-1.png)

Three elapsed-time repetitions per method on the vignette’s tiny
generated workload; diamonds are medians.

**What to notice.** Repeated timings vary even on a tiny fixed input.
The plot supports only this render’s workload and environment; it does
not rank methods for a different volume, K, parameter set, or machine.

## What belongs in the receipt?

``` r

list(
  package_version = as.character(packageVersion("neurocluster")),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  platform = R.version$platform,
  dimensions = syn$dims,
  included_voxels = length(syn$truth),
  time_points = ncol(syn$patterns),
  requested_k = 4L,
  repetitions = 3L,
  run_order = benchmark_order
)
#> $package_version
#> [1] "0.1.0"
#> 
#> $r_version
#> [1] "4.6.1"
#> 
#> $platform
#> [1] "x86_64-pc-linux-gnu"
#> 
#> $dimensions
#> [1] 10 10  4
#> 
#> $included_voxels
#> [1] 400
#> 
#> $time_points
#> [1] 18
#> 
#> $requested_k
#> [1] 4
#> 
#> $repetitions
#> [1] 3
#> 
#> $run_order
#> [1] "snic" "slic" "slic" "snic" "snic" "slic"
```

For a real benchmark, also record the commit, CPU, RAM, operating
system, thread backend and count, full method calls, warm-up policy,
allocation or peak RSS measurement method, and every repetition.
Predefine quality estimands when comparing approximation settings;
elapsed time alone cannot establish a useful clustering. See
[Performance and
memory](https://bbuchsbaum.github.io/neurocluster/articles/performance-memory.md)
and [Validate and compare
clusterings](https://bbuchsbaum.github.io/neurocluster/articles/validate-compare.md).
