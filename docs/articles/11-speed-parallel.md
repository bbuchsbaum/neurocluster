# Speed up and parallelize

## Threads

``` r
RcppParallel::setThreadOptions(numThreads = 2)
```

## Method choices

- Large N: `slice_msf` (slice-wise MSF) or `flash3d` (hash/DCT); both
  avoid global iteration over all voxels.
- Preserve K: [`cluster4d_slic()`](../reference/cluster4d_slic.md) with
  `preserve_k = TRUE`.
- Iterative refinement: [`supervoxels()`](../reference/supervoxels.md)
  performs multiple reassignment iterations; uses RcppParallel for
  updates.

## Grain size

- For `supervoxels`, a practical start is
  `grain_size = max(100, nvox / (threads * 10))`. Smaller values improve
  balancing at some overhead.

## Memory tips

- [`snic()`](../reference/snic.md) is single-pass with low memory
  overhead.
- Reduce time points or compress features if memory-bound.

## Toy runtime illustration (small N)

Illustrative timing on a small, deterministic dataset (3-band
synthetic). Absolute times are not representative of large datasets but
show relative ordering on the same volume.

``` r
toy <- make_block_synthetic(dims = c(14, 14, 2), ntime = 40, noise = 0.1, seed = 321)

run_and_time <- function(tag, fn) {
  t <- system.time(fn())
  unname(t["elapsed"])
}

elapsed <- c(
  snic = run_and_time("snic", function() snic(toy$vec, toy$mask, K = 3, compactness = 3, max_iter = 3)),
  slice_msf = run_and_time("slice_msf", function() slice_msf(toy$vec, toy$mask,
                                                            target_k_global = 3,
                                                            r = 8, min_size = 6,
                                                            compactness = 3,
                                                            num_runs = 1,
                                                            stitch_z = FALSE)),
  supervoxels = run_and_time("supervoxels", function() supervoxels(toy$vec, toy$mask,
                                                                  K = 3,
                                                                  alpha = 0.6,
                                                                  connectivity = 6,
                                                                  iterations = 5)),
  slic = run_and_time("slic", function() cluster4d(toy$vec, toy$mask,
                                                  n_clusters = 3,
                                                  method = "slic",
                                                  spatial_weight = 0.5,
                                                  connectivity = 26,
                                                  max_iterations = 5,
                                                  preserve_k = TRUE))
)

barplot(elapsed, ylab = "seconds", main = "Toy runtimes (small N)")
```

![Bar plot of elapsed seconds for several methods on tiny
data.](11-speed-parallel_files/figure-html/fig-runtime-1.png)

Illustrative runtimes on a 3-band synthetic (small N).

includes: in_header: \|-
