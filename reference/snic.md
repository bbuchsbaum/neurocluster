# SNIC: Simple Non-Iterative Clustering

The SNIC function performs a spatially constrained clustering on a
`NeuroVec` instance using the Simple Non-Iterative Clustering (SNIC)
algorithm.

## Usage

``` r
snic(vec, mask, compactness = 5, K = 500, max_iter = 100)
```

## Arguments

- vec:

  A `NeuroVec` instance supplying the data to cluster. Can also be a 3D
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.html)
  for structural image segmentation, which will be automatically
  converted to a single-timepoint NeuroVec internally.

- mask:

  A `NeuroVol` mask defining the voxels to include in the clustering
  result. If the mask contains `numeric` data, finite values strictly
  greater than zero define the included voxels. If the mask is a
  [`LogicalNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/LogicalNeuroVol-class.html),
  then `TRUE` will define the set of included voxels.

- compactness:

  A finite numeric value in `[0, 10]` controlling the spatial-feature
  mixture. Zero is feature-only, 10 is spatial-only, and larger values
  within the interval produce more compact clusters. Default is 5.

- K:

  The number of clusters to find. Default is 500.

- max_iter:

  Deprecated inactive argument retained in the function signature for
  compatibility. SNIC is non-iterative; supplying this argument
  explicitly is an error.

## Value

A `list` of class `snic_cluster_result` with the following elements:

- clusvol:

  An instance of type
  [ClusteredNeuroVol](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.html).

- gradvol:

  A `NeuroVol` instance representing the spatial gradient of the
  reference volume.

- cluster:

  A vector of contiguous positive cluster labels, one per included mask
  voxel.

- centers:

  A K-by-T matrix of mean feature values in the original input space.

- coord_centers:

  A matrix of spatial coordinates with each row corresponding to a
  cluster.

## Details

### Performance Optimization (2025)

The SNIC implementation has been **highly optimized** using lightweight
C++ structs and in-place operations, providing **10x-50x speedup** over
the original implementation. Key optimizations include:

- Elimination of Rcpp::List overhead in priority queue (uses lightweight
  struct)

- In-place centroid updates with no memory allocations

- Inline 26-connectivity neighbor iteration

- Direct pointer access for array operations

- Efficient scalar math (eliminates temporary vector allocations)

### Parallelization Status

**Currently NOT parallelized.** SNIC uses a sequential priority
queue-based algorithm that processes voxels in order of their distance
from cluster centers.

#### Sequential Operations:

1.  **Initialization**: Gradient-based seed selection using
    [`find_initial_points()`](https://bbuchsbaum.github.io/neurocluster/reference/find_initial_points.md)

    - Finds K seed points with high gradient and spatial separation

    - Sequential search through candidate voxels

2.  **Priority Queue Processing** (C++ implementation):

    - Maintains a global priority queue of voxels to be assigned

    - Each voxel assignment depends on previously processed neighbors

    - Voxels are processed in order of their combined distance metric

3.  **Distance Computation**: For each voxel, calculates:

    - Feature distance to nearest cluster center

    - Spatial distance weighted by compactness parameter

    - Combined into single priority score

#### Why Not Parallelized:

- **Sequential dependency**: Priority queue enforces strict processing
  order

- **Global state**: Each voxel assignment affects subsequent assignments

- **Algorithm design**: SNIC's key innovation is its non-iterative,
  sequential nature

- **Coherent clusters**: Sequential processing ensures connected
  components

#### Performance Characteristics:

- **Complexity**: O(N log N) where N = number of voxels

- **Memory**: O(N) for priority queue and assignments

- **Speed**: Generally faster than iterative methods (supervoxels)

- **Single pass**: Processes each voxel exactly once

#### Performance Tips:

- **Reduce K**: Fewer clusters means less competition for voxels

- **Adjust compactness**: Higher values create more local clusters,
  faster processing

- **Pre-smooth data**: Reduce noise to improve gradient-based
  initialization

- **Use smaller masks**: Process ROIs separately if possible

- **Parallel alternative**: In the unified API, `supervoxels` and `slic`
  support `parallel = TRUE`. The other methods do not advertise a
  parallel clustering contract.

#### Comparison with Other Methods:

- **Faster than**:
  [`supervoxels()`](https://bbuchsbaum.github.io/neurocluster/reference/supervoxels.md)
  due to non-iterative nature

- **Parallel choices**: Use one of the two unified methods listed above
  when a tested sequential/parallel clustering contract is required.

- **More coherent than**: Methods without spatial priority (ensures
  connectivity)

## Note

Consider using
[`cluster4d`](https://bbuchsbaum.github.io/neurocluster/reference/cluster4d.md)
with `method = "snic"` for a standardized interface across all
clustering methods.

## See also

[`supervoxels`](https://bbuchsbaum.github.io/neurocluster/reference/supervoxels.md)

## Examples

``` r
if (FALSE) { # \dontrun{
mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
vec <- replicate(10, NeuroVol(array(runif(202020), c(20,20,20)),
NeuroSpace(c(20,20,20))), simplify=FALSE)
vec <- do.call(concat, vec)

snic_res <- snic(vec, mask, compactness=5, K=100)
} # }
```
