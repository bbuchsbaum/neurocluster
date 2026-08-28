# Choose the right neurocluster interface

`neurocluster` has one main workflow and several specialist surfaces. If
you are starting an analysis, use
[`cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md):
it gives every method the same input checks and the same result
contract. Move to a direct method function only when you need controls
that the unified interface intentionally does not standardize.

## What is the shortest successful path?

| Task | Start here | What comes back |
|----|----|----|
| Fit one spatial clustering | [`cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md) | A validated-shape `cluster4d_result` |
| Check a fitted result against its inputs | [`validate_cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/validate_cluster4d.md) | Errors, warnings, and a compact summary |
| Compare compatible fits | [`compare_cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/compare_cluster4d.md) | Explicit size, dispersion, coherence, or agreement estimands |
| Choose a defensible setting | [`suggest_cluster4d_params()`](https://bbuchsbaum.github.io/neurocluster/reference/suggest_cluster4d_params.md) | A starting point to validate on your own data |
| Create a reproducible example | [`generate_synthetic_volume()`](https://bbuchsbaum.github.io/neurocluster/reference/generate_synthetic_volume.md) | Data, mask, and known latent labels |

The complete first-use journey is in [Getting
started](https://bbuchsbaum.github.io/neurocluster/articles/getting-started.md).
Parameter effects and validation are covered by [Choose
parameters](https://bbuchsbaum.github.io/neurocluster/articles/choose-parameters.md)
and [Validate and compare
clusterings](https://bbuchsbaum.github.io/neurocluster/articles/validate-compare.md).

## Which methods share the unified interface?

``` r

methods <- eval(formals(cluster4d)$method)
paste(methods, collapse = ", ")
#> [1] "supervoxels, snic, slic, corr_slic, brs_slic, slice_msf, flash3d, g3s, rena, rena_plus, mcl, acsc, commute"
```

The list above is read from the installed function, so it cannot
silently drift from the current API. Method names select
implementations; they are not a performance ranking or endorsement. In
particular, `slice_msf` remains experimental, is excluded from automatic
recommendations, and should only be used for explicit evaluation under
its [documented
contract](https://bbuchsbaum.github.io/neurocluster/articles/slice-msf.md).
See [Compare clustering
methods](https://bbuchsbaum.github.io/neurocluster/articles/compare-methods.md)
for a scoped comparison and [Method deep
dives](https://bbuchsbaum.github.io/neurocluster/articles/method-deep-dives.md)
for algorithmic distinctions.

## When should you use a direct function?

Typed method-specific wrappers expose native controls while retaining
the standard result contract. Examples include
[`cluster4d_snic()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_snic.md),
[`cluster4d_slic()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_slic.md),
[`cluster4d_g3s()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_g3s.md),
[`cluster4d_rena()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_rena.md),
[`cluster4d_mcl()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_mcl.md),
and
[`cluster4d_acsc()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_acsc.md).

The current
[`cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md)
documentation classifies
[`supervoxels()`](https://bbuchsbaum.github.io/neurocluster/reference/supervoxels.md),
[`snic()`](https://bbuchsbaum.github.io/neurocluster/reference/snic.md),
[`slic4d_supervoxels()`](https://bbuchsbaum.github.io/neurocluster/reference/slic4d_supervoxels.md),
[`slice_msf()`](https://bbuchsbaum.github.io/neurocluster/reference/slice_msf.md),
and
[`supervoxels_flash3d()`](https://bbuchsbaum.github.io/neurocluster/reference/supervoxels_flash3d.md)
as legacy, deprecated entry points. They remain exported for
compatibility, but new code should use
[`cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md)
or a typed `cluster4d_*()` wrapper. Arguments and defaults need not
match across lifecycle generations. The mapping guide is [Direct method
APIs and
cluster4d](https://bbuchsbaum.github.io/neurocluster/articles/legacy-apis.md).

## What should you inspect after fitting?

Use the result rather than reconstructing state from method internals:

``` r

fit$labels          # labels in included-mask voxel order
fit$actual_k        # number of labels actually produced
fit$centers         # cluster means in the original feature space
fit$coord_centers   # cluster centroids in physical millimetres
fit$parameters      # recorded fitting parameters
fit$provenance      # geometry, mask, feature, and label provenance
```

[Read clustering
results](https://bbuchsbaum.github.io/neurocluster/articles/objects-results.md)
explains these fields and shows the validation receipt that should
accompany a saved or exported fit.
