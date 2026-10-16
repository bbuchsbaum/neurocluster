# neurocluster 0.1.0.9000

## Documentation and interface fixes

* `?cluster4d` now lists all 13 `method` values (adding `corr_slic`,
  `rena_plus` and `commute`) and no longer labels `g3s` as "recommended";
  method order is stated not to be a ranking. The "Algorithm Comparison"
  table now covers every method, reports measured speed tiers (SNIC is
  listed as slow on whole-brain masks because its k-NN gradient seeding
  dominates run time), and states which common arguments each method uses.
  The quick `cluster4d()` example now runs during checks (#1).
* `cluster4d_rena()` and `cluster4d_rena_plus()` no longer show an inactive
  `spatial_weight = 0.5` default in their usage. `spatial_weight` is no
  longer a formal argument; supplying it (now through `...`) is still an
  error. Code that passed a fourth positional argument intending
  `spatial_weight` now binds it to `max_iterations` (`cluster4d_rena()`) or
  `r` (`cluster4d_rena_plus()`) (#1).
