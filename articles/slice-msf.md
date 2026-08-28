# Evaluate Slice-MSF (experimental)

> **Status: experimental and not recommended automatically.** Slice-MSF
> remains available for explicit evaluation, but
> [`suggest_cluster4d_params()`](https://bbuchsbaum.github.io/neurocluster/reference/suggest_cluster4d_params.md)
> will not select it. Earlier implementations could collapse distinct
> feature vectors to zero-distance edges through reliability weighting
> and could produce singleton-dominated exact-K partitions. Those paths
> have been rejected or replaced, have regression gates, and pass local
> package recertification at this source revision. Broader real-data
> validation and hosted release evidence are still pending. Validate the
> method against a named estimand on your own data before considering
> downstream use.

## Which interface are you evaluating?

The direct and unified interfaces answer different cluster-count
questions.
[`slice_msf()`](https://bbuchsbaum.github.io/neurocluster/reference/slice_msf.md)
defaults to its natural Felzenszwalb-Huttenlocher (FH) segmentation.
`cluster4d(..., method = "slice_msf")` always requests exactly
`n_clusters` and runs the shared connected exact-K repair when
necessary.

| Control | Direct [`slice_msf()`](https://bbuchsbaum.github.io/neurocluster/reference/slice_msf.md) default | Unified [`cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md) Slice-MSF default |
|----|---:|---:|
| Count policy | natural (`target_k_global = -1`) | exact requested `n_clusters` |
| DCT mode count | `r` = 12 | `r` = 12 |
| FH control | `compactness = 5` | `spatial_weight = 0.2`, mapped to compactness 2 |
| Minimum size | 80 | automatic from expected parcel size, between 2 and 20 |
| Base runs | 3 | 1 |
| Consensus | on | off |
| Cross-slice edges | on | on |
| Reliability weight | `gamma = 0` | `gamma = 0` |
| Axial sketch smoothing | `z_mult = 0` | `z_mult = 0` |

These differences are large enough that results from the two defaults
are not interchangeable. State the entry point and effective controls in
any evaluation.

## What similarity does Slice-MSF use?

For each included voxel, Slice-MSF linearly detrends and z-scores the
time series, projects it onto non-DC DCT-II modes, applies any
run-specific positive mode weights, and normalizes the sketch to unit
length. If `z_mult > 0`, each mode is smoothed along the axial direction
and the resulting voxel sketch is renormalized. An adjacent edge then
has feature dissimilarity

``` math
d(i,j) = 1 - \operatorname{clamp}(u_i^\top u_j, -1, 1).
```

A zero-information sketch remains the zero vector rather than inventing
a direction; its edge dissimilarity is therefore one to every neighbor.

The signed adjacent-pair temporal correlation returned by
[`slice_msf_single()`](https://bbuchsbaum.github.io/neurocluster/reference/slice_msf_single.md)
is diagnostic only; it does not weight or mask an edge. `gamma` is
retained for API compatibility but must be zero. A positive value fails
closed with condition class `slice_msf_unsupported_gamma` rather than
changing the geometry.

The unified `spatial_weight` is also easy to misread. For Slice-MSF it
is not a convex feature/spatial blend. The wrapper computes
`compactness = 10 * spatial_weight`, and the direct implementation then
uses FH scale `2 / (compactness + 1)`. Larger values therefore reduce
the FH component scale; feature distance remains the cosine
dissimilarity above.

## Compare natural and requested exact K

Use one deterministic fixture and hold all controls except the count
policy fixed. The known labels here are an evaluation aid, not a claim
that this synthetic scenario represents an empirical dataset.

``` r

syn <- generate_synthetic_volume(
  scenario = "gaussian_blobs", dims = c(12, 12, 4),
  n_clusters = 4, n_time = 32, noise_sd = 0.05, seed = 42
)
```

``` r

natural <- slice_msf(
  syn$vec, syn$mask, target_k_global = -1,
  r = 8, compactness = 4, min_size = 2,
  num_runs = 1, consensus = FALSE,
  stitch_z = TRUE, nbhd = 8, gamma = 0, z_mult = 0, seed = 17
)

exact <- cluster4d(
  syn$vec, syn$mask, n_clusters = 4, method = "slice_msf",
  spatial_weight = 0.4, min_size = 2, r = 8,
  num_runs = 1, consensus = FALSE,
  stitch_z = TRUE, gamma = 0, z_mult = 0, seed = 17,
  verbose = FALSE
)
```

| mode | natural_K | requested_K | final_K | smallest | largest | repair_direction | repair_ran | ARI_to_fixture_truth |
|:---|---:|---:|---:|---:|---:|:---|:---|---:|
| natural | 6 | NA | 6 | 2 | 220 | none | FALSE | 0.96 |
| requested exact K | 6 | 4 | 4 | 44 | 231 | merge | TRUE | 1.00 |

Count, size, repair, and fixture-specific accuracy diagnostics. {.table}

Natural mode applies no count repair. Exact mode first makes each label
a connected component, then uses adjacent Ward merges when it is over
target or deterministic feature-weighted minimum-spanning-tree
bisections when it is under target. Candidate ties are resolved in
stable label and voxel order. The result is therefore reproducible for
the same data, mask, controls, and package build.

The table reports ARI alongside count and size diagnostics so a
favorable score cannot hide a structurally invalid partition. On real
data there is no known truth label; replace this column with a
registered endpoint such as held-out homogeneity, cross-session
stability, or agreement with an external reference.

## How is feasibility enforced?

For exact K, each final parcel must be connected under the selected
graph and contain at least `min_size` voxels. A disconnected mask
imposes a minimum K equal to its number of connected components. For
component sizes $`n_c`$, the minimum-size constraint imposes maximum K
$`\sum_c \lfloor n_c / \text{min_size} \rfloor`$. A request outside
those bounds fails instead of silently relaxing K or parcel size.

``` r

infeasible <- tryCatch(
  cluster4d(
    syn$vec, syn$mask, n_clusters = 4, method = "slice_msf",
    spatial_weight = 0.4, min_size = 200, r = 8,
    num_runs = 1, consensus = FALSE, seed = 17, verbose = FALSE
  ),
  cluster4d_exact_k_infeasible = identity
)

data.frame(
  reason = infeasible$reason,
  requested_K = infeasible$target_k,
  minimum_K = infeasible$minimum_k,
  maximum_K = infeasible$maximum_k,
  voxels = infeasible$n_voxels,
  minimum_size = infeasible$min_cluster_size,
  natural_K = infeasible$current_k
)
#>                 reason requested_K minimum_K maximum_K voxels minimum_size
#> 1 minimum_cluster_size           4         1         2    576          200
#>   natural_K
#> 1         2
```

The structured condition exposes `reason`, `target_k`, `minimum_k`,
`maximum_k`, `n_voxels`, `min_cluster_size`, and `current_k`. Catch that
class when a pipeline needs to record infeasibility separately from
implementation errors. A global exact-K request also requires
`stitch_z = TRUE`; independent slices do not define one connected
volumetric repair graph.

## Which provenance should be retained?

Both count modes populate `metadata$exact_k_repair`. Retain at least:

- `natural_k`, `requested_k`, `direction`, and `ran`;
- `pre_sizes`, `post_sizes`, `min_cluster_size`, and `connectivity`;
- `effective_gamma`, `effective_z_mult`, `seed`, `num_runs`, and
  `consensus`;
- `operations` when exact-K repair ran.

The operation list can be long, so reports should summarize it rather
than print every merge or split. Preserve the full result object as the
audit artifact. For multiple runs, also retain `metadata$ensemble$runs`,
which records the DCT frequencies and weights for each seeded base
partition. Consensus weights runs uniformly; temporal-smoothness
diagnostics do not alter those weights.

## What should you test before using it?

Treat the method as a candidate estimand implementation, not a speed
shortcut. Before use, define the included mask, connectivity, count
policy, minimum size, feature controls, seed, and scientific endpoint.
Then verify determinism, connectedness, cluster-size constraints,
perturbation behavior, and held-out or external validity on
representative data. A successful synthetic example is a smoke test, not
release evidence.

Continue to [Build consensus and stitch across
slices](https://bbuchsbaum.github.io/neurocluster/articles/consensus-stitch.md)
for the ensemble/topology controls and [Validate and compare
clusterings](https://bbuchsbaum.github.io/neurocluster/articles/validate-compare.md)
for a comparison workflow.
