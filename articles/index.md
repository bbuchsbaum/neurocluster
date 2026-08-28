# Articles

### Tutorials

- [Getting started with
  neurocluster](https://bbuchsbaum.github.io/neurocluster/articles/getting-started.md):

  A first spatially constrained clustering, from data to a checked
  result.

- [Compare clustering
  methods](https://bbuchsbaum.github.io/neurocluster/articles/compare-methods.md):

  Compare valid partitions using explicit, shared estimands.

- [End-to-end: From NIfTI to
  clusters](https://bbuchsbaum.github.io/neurocluster/articles/end-to-end-export.md):

  Read 4D data and a mask, fit a checked partition, and round-trip its
  NIfTI label map.

### How-tos

- [Choose parameters for your
  data](https://bbuchsbaum.github.io/neurocluster/articles/choose-parameters.md):

  Learn how parcel count and spatial weighting change the estimand.

- [Use parallel execution
  responsibly](https://bbuchsbaum.github.io/neurocluster/articles/speed-parallel.md):

  Select per-call parallel controls and verify what the method actually
  used.

- [Visualize and export
  results](https://bbuchsbaum.github.io/neurocluster/articles/visualize-export.md):

  Read parcel maps, inspect temporal summaries, and verify NIfTI export.

- [Build consensus and stitch across
  slices](https://bbuchsbaum.github.io/neurocluster/articles/consensus-stitch.md):

  Understand the two independent Slice-MSF controls for ensembles and
  volumetric topology.

- [Evaluate Slice-MSF
  (experimental)](https://bbuchsbaum.github.io/neurocluster/articles/slice-msf.md):

  Inspect Slice-MSF’s feature geometry, natural and exact-K modes,
  feasibility gates, and provenance before deciding whether to evaluate
  it.

- [Validate and compare
  clusterings](https://bbuchsbaum.github.io/neurocluster/articles/validate-compare.md):

  Separate structural validity, single-partition diagnostics, and
  pairwise agreement.

- [Direct method APIs and
  cluster4d](https://bbuchsbaum.github.io/neurocluster/articles/legacy-apis.md):

  Understand what the unified interface maps and what remains
  method-specific.

### Explanations

- [Understand spatially constrained
  clustering](https://bbuchsbaum.github.io/neurocluster/articles/spatially-constrained-clustering.md):

  See how feature similarity, distance, and connectivity shape a parcel
  map.

- [How the clustering methods
  differ](https://bbuchsbaum.github.io/neurocluster/articles/method-deep-dives.md):

  Compare algorithmic work units and controls without turning them into
  universal rankings.

- [Reason about performance and
  memory](https://bbuchsbaum.github.io/neurocluster/articles/performance-memory.md):

  Define what can be estimated before measuring method-specific runtime
  and memory.

### Benchmarks

- [Benchmark with a reproducible
  receipt](https://bbuchsbaum.github.io/neurocluster/articles/benchmarks.md):

  Run a scoped micro-comparison and report enough context to reproduce
  it.

### Reference Maps

- [Choose the right neurocluster
  interface](https://bbuchsbaum.github.io/neurocluster/articles/api-overview.md):

  Start with cluster4d, then move to validation, comparison, or direct
  method APIs.

- [Read, validate, and save clustering
  results](https://bbuchsbaum.github.io/neurocluster/articles/objects-results.md):

  Interpret the cluster4d_result contract and preserve its provenance.
