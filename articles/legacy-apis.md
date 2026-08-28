# Direct method APIs and cluster4d

The package keeps
[`cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md),
typed `cluster4d_*()` wrappers, and older method-native functions. The
current
[`cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md)
documentation classifies
[`supervoxels()`](https://bbuchsbaum.github.io/neurocluster/reference/supervoxels.md),
[`snic()`](https://bbuchsbaum.github.io/neurocluster/reference/snic.md),
[`slic4d_supervoxels()`](https://bbuchsbaum.github.io/neurocluster/reference/slic4d_supervoxels.md),
[`slice_msf()`](https://bbuchsbaum.github.io/neurocluster/reference/slice_msf.md),
and
[`supervoxels_flash3d()`](https://bbuchsbaum.github.io/neurocluster/reference/supervoxels_flash3d.md)
as legacy/deprecated entry points. They remain exported for
compatibility, but new code should use the unified interface or its
typed wrappers. Equivalence should never be assumed unless mapped
parameters are identical.

[`slice_msf()`](https://bbuchsbaum.github.io/neurocluster/reference/slice_msf.md)
and
[`cluster4d_slice_msf()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_slice_msf.md)
are additionally experimental. The typed wrapper standardizes the result
and requests exact K, but neither entry point is currently recommended
automatically or for production use. See the [Slice-MSF evaluation
contract](https://bbuchsbaum.github.io/neurocluster/articles/slice-msf.md)
before evaluating either one.

## What does the unified interface map?

| Unified method | Legacy entry point | Preferred typed wrapper | Important mapping |
|----|----|----|----|
| `"supervoxels"` | [`supervoxels()`](https://bbuchsbaum.github.io/neurocluster/reference/supervoxels.md) | [`cluster4d_supervoxels()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_supervoxels.md) | `alpha = 1 - spatial_weight` |
| `"snic"` | [`snic()`](https://bbuchsbaum.github.io/neurocluster/reference/snic.md) | [`cluster4d_snic()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_snic.md) | `compactness = 10 * spatial_weight` |
| `"slic"` | [`slic4d_supervoxels()`](https://bbuchsbaum.github.io/neurocluster/reference/slic4d_supervoxels.md) | [`cluster4d_slic()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_slic.md) | `compactness = 20 * spatial_weight` |
| `"slice_msf"` (experimental) | [`slice_msf()`](https://bbuchsbaum.github.io/neurocluster/reference/slice_msf.md) | [`cluster4d_slice_msf()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_slice_msf.md) | natural versus exact K, connectivity, and FH scale are translated; positive `gamma` is rejected |
| `"flash3d"` | [`supervoxels_flash3d()`](https://bbuchsbaum.github.io/neurocluster/reference/supervoxels_flash3d.md) | [`cluster4d_flash3d()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_flash3d.md) | `lambda_s = spatial_weight`; temporal weight is complementary unless supplied |
| `"g3s"` | none (`g3s()` does not exist) | [`cluster4d_g3s()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d_g3s.md) | `alpha = 1 - spatial_weight` |

Those statements describe the current wrappers. Defaults and
method-specific arguments can still differ, so record the actual
[`cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md)
or typed-wrapper call—not an informal translation.

## Can a mapping be checked rather than trusted?

This small compatibility example checks the exact legacy SNIC mapping on
one generated volume. It explains an existing script; it is not a
recommendation to start new code with
[`snic()`](https://bbuchsbaum.github.io/neurocluster/reference/snic.md).

``` r

unified <- cluster4d(
  syn$vec, syn$mask, n_clusters = 4, method = "snic",
  spatial_weight = 0.4
)
direct <- snic(syn$vec, syn$mask, K = 4, compactness = 4)
```

For this input and package version, the labels are identical. That is a
check of this mapping—not a blanket promise that every unified/direct
pair is bit-identical.

## Which interface should you publish?

Use
[`cluster4d()`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md)
when methods are compared under a shared contract. Use a typed
`cluster4d_*()` wrapper when the protocol depends on a native parameter.
Treat the legacy functions as compatibility surfaces. In every case,
save the result’s `parameters` and `provenance`, validate it against the
original `vec` and `mask`, and report the actual rather than merely
requested cluster count. See [Read clustering
results](https://bbuchsbaum.github.io/neurocluster/articles/objects-results.md).
