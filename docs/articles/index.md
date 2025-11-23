# Articles

### Tutorials

- [Getting started with neurocluster](01-getting-started.md):

  Quick tour of the core workflow.

- [Compare clustering methods](02-compare-methods.md):

  Run multiple methods on the same data and compare outputs.

- [End-to-end: From NIfTI to clusters](03-end-to-end-export.md):

  Load data, cluster, visualize, and export results.

### How-tos

- [Choose parameters for your data](10-choose-parameters.md):

  Use suggestions and tune K, spatial weighting, and connectivity.

- [Speed up and parallelize](11-speed-parallel.md):

  Threads, grain size, and method choices for big data.

- [Visualize and export results](12-visualize-export.md):

  Plot cluster slices and write NIfTI outputs.

- [Consensus and slice stitching](13-consensus-stitch.md):

  Use slice_msf with consensus and stitch clusters across z-slices.

- [Validate and compare clusterings](14-validate-compare.md):

  Quality checks, cluster sizes, and basic overlaps.

- [Legacy APIs and backward compatibility](15-legacy-apis.md):

  Map old interfaces to the unified cluster4d entry point.

### Explanations

- [Spatially constrained
  clustering](20-spatially-constrained-clustering.md):

  Concepts: spatial vs. feature weighting, connectivity, and artifacts.

- [Method deep dives](21-method-deep-dives.md):

  Supervoxels, SNIC, SLIC, slice-MSF, and FLASH-3D—ideas and tradeoffs.

- [Performance and memory tradeoffs](22-performance-memory.md):

  Scaling rules, memory estimates, and method choices.

### Reference Maps

- [API overview](30-api-overview.md):

  Group exported functions by task and link to details.

- [Objects and results](31-objects-results.md):

  cluster4d_result structure, S3 methods, and neuroim2 objects.
