# Correlation-SLIC for fMRI: Design Notes (neurocluster)

## Purpose

Implementation notes for adding a new 4D clustering method that is SLIC-like in 3D space and correlation-driven in fMRI feature space, optimized for whole-brain scale.

Primary goal: preserve SLIC locality/compactness while replacing expensive high-dimensional correlation comparisons with a low-dimensional embedding that approximates correlation.

Scope of this note:
- algorithm design,
- speed-critical data flow,
- neuroim2 integration choices,
- R/Rcpp boundary,
- connectivity guarantees,
- defaults/tuning guidance,
- staged implementation plan.

## Design Constraints in This Package

This package already uses `neuroim2` for data structures and extraction paths (`NeuroVec`, `NeuroVol`, `series()`, `index_to_coord()`, `spacing()`, `ClusteredNeuroVol`).

For this method, we should:
- keep user-facing API consistent with `cluster4d_*` wrappers,
- use `neuroim2` objects at the R boundary,
- keep heavy compute in C++ (`src/`) via Rcpp,
- avoid introducing a separate I/O/data model just for this method,
- remain compatible with existing `cluster4d_result` creation flow.

## High-Level Algorithm

1. Extract masked voxel time series `X` from 4D data (`N x T`).
2. Build a small embedding per voxel (`d=32..128`) where dot products approximate Pearson correlation.
3. Run SLIC-style local iterative clustering in 3D with feature term + spatial compactness term.
4. Enforce connectivity so each final label is a single connected component (no islands).
5. Return labels + centers in standard `cluster4d_result` format.

Core distance used during clustering:

`D(v,k) = (1 - <f_v, mu_f_k>) + alpha * ||s_v - mu_s_k||^2 / S^2`

where:
- `f_v` and `mu_f_k` are unit-norm embedding vectors,
- `s_v` and `mu_s_k` are voxel and center spatial coordinates,
- `S ~= (N/K)^(1/3)` is expected spacing,
- `alpha` controls compactness.

## Correlation -> Dot Product Embedding

For demeaned time series, correlation equals cosine similarity:
- `rho(x,y) = <x,y> / (||x|| ||y||)`.

Direct repeated correlation in `T` dimensions is too expensive. Replace with a hash-based random projection (CountSketch-like):
- each timepoint `t` hashes to bin `h[t] in {0..d-1}` and sign `s[t] in {-1,+1}`,
- accumulate into `u_j = sum_{t: h[t]=j} s[t] * (x_t - mean(x))`,
- normalize `u` to unit length, yielding `f`.

This makes similarity comparisons O(`d`) instead of O(`T`) during iterative clustering.

## One-Pass Demeaning Trick

To avoid two passes per voxel:

Accumulate for each voxel:
- `sum_x = sum_t x_t`
- `u_raw[j] = sum_{t: h[t]=j} s[t] * x_t`

Precompute once per run:
- `sumS[j] = sum_{t: h[t]=j} s[t]`

Then:
- `mean = sum_x / T`
- `u[j] = u_raw[j] - mean * sumS[j]`

Finally normalize `u`.

This keeps embedding build linear in `N*T` with low constants.

## Multi-Run Strategy

Preferred default: average embeddings across runs, then cluster once.

For run `r`, compute unit embedding `f_v^(r)`. Aggregate:
- `g_v = normalize(sum_r w_r * f_v^(r))`

Advantages:
- avoids label-matching across runs,
- parallelizable run-wise,
- stable and simple.

Alternative for harder run heterogeneity:
- concatenate run embeddings then normalize (`[f^(1),...,f^(R)]`),
- higher compute/memory.

## SLIC-Like Local Clustering Details

### Initialization

Use grid seeding over mask bounding box with approximate step `S`.

For each seed point:
- snap to nearest in-mask voxel within small search radius (1-2 voxels),
- initialize center feature from that voxel embedding.

If seed count differs from `K`:
- downsample or top-up deterministically.

### Local Candidate Restriction (critical for speed)

Use spatial bins of size `S`:
- assign centers to bins each iteration,
- each voxel only evaluates centers in its bin + 26 neighboring bins,
- candidate count stays small and roughly constant.

### Update

After assignment:
- spatial center = mean voxel coordinates,
- feature center = mean embedding, then renormalize to unit norm.

Use 3-6 iterations with optional early stop when label changes are small.

## Connectivity Enforcement (No Islands)

Final pass:
1. Find connected components within each label via BFS/DFS on masked voxels (6-neighbor default).
2. Keep largest component per label.
3. Reassign smaller components to best neighboring label.
4. Optionally merge components below `min_size`.

Recommended reassignment rule:
- primary: neighbor label with highest boundary contact,
- optional tie-break: better feature-center similarity.

Complexity is O(`N`) and guarantees contiguous final parcels.

## neuroim2 Integration Plan (R-Side)

Follow existing wrapper style (`cluster4d_*`, `slic4d_supervoxels`):

1. Validate inputs with `validate_cluster4d_inputs()`.
2. Use `mask_idx <- which(mask > 0)`.
3. Use `coords <- as.matrix(index_to_coord(mask, mask_idx))`.
4. Use `feat <- t(as.matrix(series(vec, mask_idx)))` for `N x T`.
5. Call C++ core with compact arguments:
   - `feat` (`N x T`, numeric),
   - `coords` (`N x 3`),
   - mask linear indices (0-based),
   - dims/spacings,
   - `K`, `d`, `alpha`, `iters`, `connectivity`, `min_size`, thread count,
   - optional run weights / seed.
6. Rebuild `ClusteredNeuroVol` and produce `cluster4d_result` via `create_cluster4d_result()`.

Important: use `neuroim2` for object semantics at boundaries, but keep all iterative kernels in C++.

## C++ Data Layout and Structures

For performance and simpler loops:

- `labels[N]` as `int32_t`.
- `coords_x/y/z[N]` as `float` (or `int16` if using voxel grid units only).
- `mask_lin_idx[N]` as `int32_t`.
- `id_map[V]` as `int32_t` mapping full linear index -> masked voxel id (`-1` outside).
- embedding accumulation in dimension-major `U[d*N]` for hash updates.
- clustering feature matrix in voxel-major `feat[N*d]` for fast per-voxel dot products.
- center arrays:
  - spatial: `cx[K], cy[K], cz[K]`,
  - features: contiguous `cf[K*d]`.

## Parallelization Strategy

OpenMP targets:
- embedding build over voxel chunks,
- assignment step over voxels,
- update step via thread-local accumulators then reduce.

Guidelines:
- avoid atomics on large `(K*d)` arrays,
- reuse allocated buffers across iterations,
- avoid thread creation/destruction inside tight loops,
- keep memory access contiguous where possible.

## Numerical and Performance Defaults

Starting defaults:
- `d = 64`,
- `iters = 5`,
- `alpha = 0.5`,
- `connectivity = 6`,
- `min_size = max(10, round(0.25 * S^3))`.

Speed-oriented variant:
- stage 1: `d=32`, `iters=3`,
- stage 2 refinement: `d=64`, `iters=2` initialized from stage-1 labels/centers.

## Complexity

Let:
- `N` masked voxels,
- `T` timepoints,
- `d` embedding dim,
- `K` clusters,
- `I` iterations,
- `C` local candidate centers per voxel.

Costs:
- embedding: O(`N*T`) dominant,
- clustering: O(`I*N*C*d`) with bounded `C` via bins,
- connectivity: O(`N`).

This supports thousands of clusters while remaining practical.

## Implementation Milestones

1. Add new method wrapper (e.g., `cluster4d_corrslic()` and `cluster4d(..., method=...)` wiring).
2. Add minimal C++ core:
   - embedding,
   - local SLIC iterations,
   - connectivity enforcement.
3. Add correctness tests:
   - dimensions/labels/counts,
   - deterministic seed behavior,
   - connectivity guarantee,
   - edge cases (`K=1`, tiny masks, small T).
4. Add performance smoke benchmark on synthetic data.
5. Tune defaults and expose advanced parameters.

## Testing Checklist

Functional:
- returns `cluster4d_result` with expected fields,
- cluster ids contiguous and 1-based in R,
- no unlabeled masked voxels,
- exactly one connected component per label after enforcement.

Stability:
- deterministic output with fixed seed,
- robust to low-variance voxels and constant time series,
- handles single-timepoint input gracefully (fallback behavior defined).

Performance:
- benchmark against current `slic`/`flash3d`/`snic` on representative `N,T`.

## Open Design Questions

1. Method naming in public API:
- add as new `method` (e.g., `"corr_slic"`),
- or fold into current `"slic"` with feature mode switch.

2. Spatial units in distance term:
- voxel index units (fast, simple),
- or physical mm units using `spacing()` (more anatomically faithful).

3. Embedding build source:
- from `series()` matrix (`N x T`) first implementation (simpler integration),
- later optimize with chunked/block extraction if memory pressure appears.

4. Connectivity neighborhood default:
- 6-neighbor for stricter topology,
- 26-neighbor for smoother shapes.

## Immediate Recommendation

Implement this method first with:
- `neuroim2` extraction (`series`, `index_to_coord`) for consistency,
- C++ kernels for embedding + local SLIC + connectivity,
- average-embedding multi-run consensus,
- conservative defaults (`d=64`, `iters=5`, `alpha=0.5`, 6-neighbor).

Then profile and optimize hot loops (dot products, update reduction, buffer reuse) before broad API expansion.
