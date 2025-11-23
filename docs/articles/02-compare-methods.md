# Compare clustering methods

## Goal

Run several methods with the same `n_clusters` and compare basic
characteristics.

Throughout this vignette we use
[`generate_synthetic_volume()`](../reference/generate_synthetic_volume.md)
to build a tiny checkerboard dataset (single axial slice, 2×3
structure). The generator provides both the noisy time series volume and
the latent labels, so we can see how close each method gets to the known
answer.

## Ground truth grid (for reference)

``` r
# Display the ground truth as a simple image
truth_array <- array(toy$truth, dim = toy$dims)
image(truth_array[,,1], col = rainbow(toy$n_clusters),
      main = "Ground Truth: 2x3 Grid", axes = FALSE)
```

![Axial slice showing a clean 2×3 grid of colored
blocks.](02-compare-methods_files/figure-html/fig-toy-ground-truth-1.png)

Ground truth 2×3 grid of clusters on the toy axial slice (labels 1–6).

## Run methods (same K)

``` r
methods <- c("supervoxels", "snic", "slic")
run_one <- function(m) {
  # Method-specific tuning for this tiny toy so that each method plausibly recovers the 2×3 grid.
  if (m == "supervoxels") {
    args <- list(
      vec = toy$vec, mask = toy$mask,
      n_clusters = 6,
      method = "supervoxels",
      spatial_weight = 0.3,          # more feature-driven on this toy
      connectivity = 6,
      max_iterations = 15,
      sigma2 = 1,                    # tighter spatial kernel
      use_gradient = FALSE
    )
  } else if (m == "snic") {
    args <- list(
      vec = toy$vec, mask = toy$mask,
      n_clusters = 6,
      method = "snic",
      spatial_weight = 0.25,         # lower compactness to respect patterns
      max_iterations = 50
    )
  } else if (m == "slic") {
    args <- list(
      vec = toy$vec, mask = toy$mask,
      n_clusters = 6,
      method = "slic",
      spatial_weight = 0.5,
      connectivity = 26,
      max_iterations = 10,
      preserve_k = TRUE,
      seed_method = "grid"
    )
  } else {
    stop("Unknown method")
  }
  out <- try(do.call(cluster4d, args), silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}
results <- setNames(lapply(methods, run_one), methods)
results_ok <- Filter(Negate(is.null), results)
```

## Axial slices by method

``` r
# Fallback in case prior chunk failed in a different environment
if (!exists("results_ok", inherits = TRUE)) {
  methods <- c("supervoxels", "snic", "slic")
  run_one <- function(m) {
    # Simple fallback with default parameters
    out <- try(cluster4d(toy$vec, toy$mask,
                         n_clusters = 6, method = m,
                         max_iterations = 5),
               silent = TRUE)
    if (inherits(out, "try-error")) NULL else out
  }
  results <- setNames(lapply(methods, run_one), methods)
  results_ok <- Filter(Negate(is.null), results)
}
n <- length(results_ok); if (n == 0) n <- 1
par(mfrow = c(1, n))
for (nm in names(results_ok)) {
  plot(results_ok[[nm]], slice = c(1, 1, 1), view = "axial")
  title(nm)
}
```

![Panels showing axial slices for available methods with arbitrary
cluster colors.](02-compare-methods_files/figure-html/fig-methods-1.png)

Toy axial view clustered by available methods (K=6; method-specific
spatial settings). Colors indicate cluster IDs (arbitrary).

![Panels showing axial slices for available methods with arbitrary
cluster colors.](02-compare-methods_files/figure-html/fig-methods-2.png)

Toy axial view clustered by available methods (K=6; method-specific
spatial settings). Colors indicate cluster IDs (arbitrary).

![Panels showing axial slices for available methods with arbitrary
cluster colors.](02-compare-methods_files/figure-html/fig-methods-3.png)

Toy axial view clustered by available methods (K=6; method-specific
spatial settings). Colors indicate cluster IDs (arbitrary).

``` r
par(mfrow = c(1, 1))
```

## Compare

``` r
if (length(results_ok) >= 1) {
  comparison <- do.call(compare_cluster4d, results_ok)
  comparison
}
comparison

# Summarize basic facts from outputs
if (length(results_ok) >= 1) {
  sizes <- lapply(results_ok, function(x) table(x$cluster))
  data.frame(
    method = names(sizes),
    n_clusters = sapply(results_ok, function(x) x$n_clusters),
    min_size = sapply(sizes, min),
    max_size = sapply(sizes, max),
    mean_size = sapply(sizes, function(t) round(mean(t), 1))
  )
}
```

## Notes

- SNIC assigns voxels in a single pass via a priority queue (see
  [`snic()`](../reference/snic.md)); clusters are connected by
  construction. Runtime depends on input size and queue operations.
- Supervoxels uses iterative reassignment with spatial/feature kernels
  (see [`supervoxels()`](../reference/supervoxels.md)); more iterations
  can change results and runtime.
- SLIC uses local search windows and can preserve the requested K when
  `preserve_k = TRUE` (see
  [`cluster4d_slic()`](../reference/cluster4d_slic.md)).

## See also

- Validate & compare: articles/validate-compare.html#checks
- Method deep dives: articles/method-deep-dives.html#slic
- Performance & memory: articles/performance-memory.html#scaling
