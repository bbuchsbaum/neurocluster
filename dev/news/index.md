# Changelog

## neurocluster 0.1.0.9000

### Performance fixes

- The sparse commute solver now uses `eigencore (>= 1.0.3)` from CRAN as
  a required import, replacing optional RSpectra. Large masks no longer
  need a separate solver installation. Sparse solves use a fixed local
  RNG seed, check original-Laplacian residuals and orthogonality, and
  disable implicit dense fallback. The minimum R version is now 4.1, as
  required by eigencore.
- [`rena_plus()`](https://bbuchsbaum.github.io/neurocluster/dev/reference/rena_plus.md)
  (`cluster4d(method = "rena_plus")`) is no longer quadratic in the
  number of voxels. The ReNA coarsening walked each row of the
  column-compressed adjacency matrix with Armadillo’s row iterator,
  which costs O(N) per row; it now walks the columns of the transpose,
  which visits the same neighbours in the same order. Labels are
  bit-identical to the previous implementation. On the 29,532-voxel
  neuroim2 mask (20 timepoints, `n_clusters = 200`) run time drops from
  118 s to 0.5 s
  ([\#3](https://github.com/bbuchsbaum/neurocluster/issues/3)).
- [`commute_cluster()`](https://bbuchsbaum.github.io/neurocluster/dev/reference/commute_cluster.md)
  (`cluster4d(method = "commute")`) gains
  `eigen_solver = c("auto", "dense", "sparse")`. The dense path formed
  an N-by-N Laplacian and a full
  [`eigen()`](https://rdrr.io/r/base/eigen.html) (O(N^2) memory, O(N^3)
  time; about 7 GB per dense copy on a 30,000-voxel mask) although only
  the `ncomp` smallest eigenpairs are used. Above 3,000 voxels `"auto"`
  now computes only those eigenpairs with shift-invert Lanczos; at or
  below 3,000 voxels results are unchanged. The initial RSpectra
  implementation finished the 29,532-voxel mask in about 5 s
  ([\#3](https://github.com/bbuchsbaum/neurocluster/issues/3)); that
  timing does not characterize the replacement eigencore backend. Dense
  and sparse solvers are compared against numerical oracles in the
  tests.

### Continuous integration

- The `Native sanitizers` workflow’s asan job had failed since August
  2026 without compiling the package: `tools/run-native-sanitizers.sh`
  preloaded the ASan runtime into GNU make, whose own `strndup()` calls
  ASan reported as heap-buffer-overflows under `strict_string_checks=1`.
  The runtime is now preloaded only into the R process that runs the
  smoke cases, and the script fails if ASan is not active there. No
  overflow was found in package code
  ([\#2](https://github.com/bbuchsbaum/neurocluster/issues/2)).

### Documentation and interface fixes

- [`?cluster4d`](https://bbuchsbaum.github.io/neurocluster/dev/reference/cluster4d.md)
  now lists all 13 `method` values (adding `corr_slic`, `rena_plus` and
  `commute`) and no longer labels `g3s` as “recommended”; method order
  is stated not to be a ranking. The “Algorithm Comparison” table now
  covers every method, reports measured speed tiers (SNIC is listed as
  slow on whole-brain masks because its k-NN gradient seeding dominates
  run time), and states which common arguments each method uses. The
  quick
  [`cluster4d()`](https://bbuchsbaum.github.io/neurocluster/dev/reference/cluster4d.md)
  example now runs during checks
  ([\#1](https://github.com/bbuchsbaum/neurocluster/issues/1)).
- [`cluster4d_rena()`](https://bbuchsbaum.github.io/neurocluster/dev/reference/cluster4d_rena.md)
  and
  [`cluster4d_rena_plus()`](https://bbuchsbaum.github.io/neurocluster/dev/reference/cluster4d_rena_plus.md)
  no longer show an inactive `spatial_weight = 0.5` default in their
  usage. `spatial_weight` is no longer a formal argument; supplying it
  (now through `...`) is still an error. Code that passed a fourth
  positional argument intending `spatial_weight` now binds it to
  `max_iterations`
  ([`cluster4d_rena()`](https://bbuchsbaum.github.io/neurocluster/dev/reference/cluster4d_rena.md))
  or `r`
  ([`cluster4d_rena_plus()`](https://bbuchsbaum.github.io/neurocluster/dev/reference/cluster4d_rena_plus.md))
  ([\#1](https://github.com/bbuchsbaum/neurocluster/issues/1)).
