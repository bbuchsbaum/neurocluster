# Use parallel execution responsibly

Parallel execution is a property of a particular method path, input
size, and backend. More threads can lose to scheduling overhead on small
masks, and a method that accepts `parallel = TRUE` may still choose a
serial path for a tiny job. Treat thread settings as part of the
benchmark receipt, not as a universal speed switch.

## Which unified methods accept `parallel`?

In the current
[`cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md)
contract, `parallel` is supported for `method = "supervoxels"` and
`method = "slic"`. An explicit `parallel = TRUE` is rejected for methods
without a parallel unified path. This fail-closed behaviour prevents a
script from appearing parallel while silently ignoring the request.

``` r

small_requested <- cluster4d(
  small_syn$vec, small_syn$mask, n_clusters = 4,
  method = "slic", parallel = TRUE,
  preserve_k = TRUE, max_iterations = 3
)
large_serial <- cluster4d(
  large_syn$vec, large_syn$mask, n_clusters = 4,
  method = "slic", parallel = FALSE,
  preserve_k = TRUE, max_iterations = 3
)
large_parallel <- cluster4d(
  large_syn$vec, large_syn$mask, n_clusters = 4,
  method = "slic", parallel = TRUE,
  preserve_k = TRUE, max_iterations = 3
)
```

| workload | masked_voxels | requested | used  | native_n_threads |
|:---------|--------------:|:----------|:------|-----------------:|
| small    |           588 | TRUE      | FALSE |                1 |
| large    |          2592 | FALSE     | FALSE |                1 |
| large    |          2592 | TRUE      | TRUE  |                0 |

Requested and effective SLIC assignment modes in this render. {.table}

The small 588-voxel call requests parallel assignment but deliberately
uses the serial kernel; the 2,592-voxel call exercises the parallel
path. The serial and parallel large-workload calls return identical
labels in this render. These are execution-contract checks, not a
speedup claim. Inspect the recorded metadata rather than inferring
thread use from elapsed time alone.

## How do direct APIs control threads?

- `supervoxels(..., parallel = TRUE, num_threads = n)` scopes an
  optional RcppParallel thread override to that call and restores the
  previous setting.
- `slic4d_supervoxels(..., n_threads = n)` passes the requested thread
  count to its native assignment operations. `0` means backend default
  and `1` requests serial assignment.
- `cluster4d(..., method = "slic", parallel = FALSE)` maps to one native
  thread. `parallel = TRUE` requests the backend default, but inputs
  below 2,000 masked voxels use the serial assignment kernel and report
  that fallback in metadata.

Prefer these per-call controls. Changing `RcppParallel` global options
in an analysis or vignette can affect unrelated work in the same R
session and makes the execution receipt harder to interpret.

## How should a speedup claim be tested?

Run serial and parallel calls on the same installed package, input, K,
method arguments, and machine. Warm both paths, repeat them in
interleaved order, and validate every returned fit. Report elapsed-time
distributions and whether labels or quality estimands changed. Record
CPU, backend, requested threads, and method metadata.

Do not tune `grain_size` from a formula presented as universally
optimal. It is a scheduling control for
[`supervoxels()`](https://bbuchsbaum.github.io/neurocluster/reference/supervoxels.md)
and should be benchmarked over a small prespecified grid on the actual
workload. The executable receipt pattern is in [Benchmark with a
reproducible
receipt](https://bbuchsbaum.github.io/neurocluster/articles/benchmarks.md);
memory boundaries are in [Performance and
memory](https://bbuchsbaum.github.io/neurocluster/articles/performance-memory.md).
