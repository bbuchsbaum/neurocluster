# G3S Optimization Checklist

This file tracks speed/accuracy optimization work for `cluster4d_g3s()` with a strict keep-or-revert policy.

## Baseline (pre-optimization reference)

- Source run: `Rscript inst/benchmarks/bench_methods.R`
- G3S aggregate (best ARI per dataset/method):
  - median elapsed: `0.076s`
  - mean elapsed: `0.139s`
  - median ARI: `0.880`
  - mean ARI: `0.784`
  - error rows: `0`

## Current Kept State

- Source run: `Rscript inst/benchmarks/bench_methods.R` (latest)
- G3S aggregate (best ARI per dataset/method):
  - median elapsed: `0.097s`
  - mean elapsed: `0.180s`
  - median ARI: `0.880`
  - mean ARI: `0.784`
  - error rows: `0`
- Note: full-suite wall times are noisy across runs on this machine; keep/revert was driven by targeted checks below.

## Keep/Revert Gate

For each patch:

1. Run targeted correctness checks:
   - `devtools::test(filter = "g3s")`
2. Run benchmark comparison:
   - minimum: focused G3S benchmark slice on shared synthetic datasets
   - preferred before merge: full `inst/benchmarks/bench_methods.R`
3. Keep only if:
   - correctness passes, and
   - speed improves clearly (or is neutral), and
   - ARI is maintained (no material regression)
4. If ambiguous, do not keep.

## Priority Queue

1. Reuse a single `kNN` graph between seed selection and propagation.
2. Reduce priority-queue churn in `g3s_propagate_cpp` via best-cost pruning.
3. Boundary-only refinement (instead of full-voxel scan per iteration).
4. Optional OpenMP for refinement and centroid update loops.
5. Reduce copy/layout overhead in feature matrix handoff.
6. Adaptive `alpha` policy for geometry-dependent robustness.
7. Stronger seed fallback (farthest-point instead of raw remainder fill).

## Progress Log

- [x] Item 1: shared `kNN`
  - Status: KEPT
  - Evidence: targeted micro-benchmark comparing old-vs-new seeding/propagation prep path showed deterministic output parity and `~1.50x` geomean speedup (`blobs_small`, `touching_blobs`, `z_layers_soft`).
- [x] Item 2: queue pruning
  - Status: REVERTED
  - Reason: correctness regression (`block_small` returned 4 labels due unlabeled `0` spillover; ARI collapsed).
- [x] Item 3: boundary-only refinement
  - Status: REVERTED
  - Reason: no unambiguous speed win in benchmark gate; improvement signal was ambiguous.
- [x] Item 4: OpenMP refinement
  - Status: KEPT
  - Evidence: full benchmark gate retained ARI (`median=0.880`, `mean=0.784`) and improved runtime vs pre-item4 kept state (`mean 0.184s -> 0.178s` on best-per-dataset G3S rows; median in same range with run-to-run jitter).
- [x] Item 5: data layout/copy reduction
  - Status: REVERTED
  - Reason: no unambiguous speed win in full benchmark gate; one attempt increased mean runtime on larger cases.
- [x] Item 6: adaptive alpha
  - Status: REVERTED
  - Reason: A/B gate (`adaptive_alpha=TRUE` vs `FALSE` on benchmark datasets) showed no clear speed gain and material quality regressions.
  - Evidence (`inst/benchmarks/g3s_adaptive_ab.csv`):
    - median speedup: `0.972x`
    - geometric mean speedup: `0.998x`
    - mean ARI delta: `-0.01754`
    - ARI nonnegative: `11/21` cases
- [x] Item 7: seed fallback improvement
  - Status: REVERTED
  - Reason: farthest-point fallback changed behavior too aggressively; quality and speed tradeoff was mixed with significant regressions on several dataset/param pairs.
  - Evidence (`inst/benchmarks/g3s_seed_fallback_ab.csv`, old vs new monkey-patch A/B):
    - median speedup (new/old): `1.000x`
    - geometric mean speedup (new/old): `0.914x`
    - mean ARI delta (new-old): `-0.03799`
    - ARI nonnegative: `12/21` cases

## Post-Checklist Experiments

- Spatial kNN backend switch (`FNN::get.knn` -> configurable backend with `RANN::nn2` default in G3S)
  - Status: KEPT
  - Files: `R/g3s.R`, `tests/testthat/test_g3s.R`
  - Evidence (`inst/benchmarks/g3s_knn_backend_ab_seeded.csv`, seeded A/B):
    - median speedup (RANN/FNN): `1.000x`
    - geometric mean speedup (RANN/FNN): `1.029x`
    - mean ARI delta (RANN-FNN): `0.00000`
    - ARI nonnegative: `21/21` cases
