# Adjacency-preserving exact-K repair

Internal machinery for converting a masked-voxel partition to exactly K
connected labels. Over-target partitions use only adjacent Ward merges;
under-target partitions split deterministic non-articulation leaves from
a spanning tree. Disconnected pieces of an input label are normalized
before either operation.

## Usage

``` r
force_exact_k(
  labels,
  feature_mat,
  K_target,
  mask,
  connectivity = 6L,
  graph_info = NULL
)
```

## Arguments

- labels:

  Positive integer labels in masked-voxel order.

- feature_mat:

  Numeric matrix with voxels in rows.

- K_target:

  Desired number of connected labels.

- mask:

  NeuroVol defining voxel geometry and inclusion.

- connectivity:

  Grid connectivity, one of 6, 18, or 26.

- graph_info:

  Optional prebuilt `.exact_k_graph()` result for the same mask and
  connectivity. Supplying it avoids rebuilding the masked-grid graph,
  which callers that already hold one would otherwise pay for twice. The
  receipt is rejected if its dimensions, included voxels, or
  connectivity do not match `mask`.

## Value

Contiguous positive integer labels with exactly `K_target` connected
components, or a `cluster4d_exact_k_infeasible` condition.
