# Recursive Nearest Agglomeration (ReNA) Clustering

Performs spatially-constrained hierarchical clustering using recursive
1-nearest neighbor aggregation. ReNA is a fast, linear-time
agglomerative algorithm that avoids the "percolation" problem of
standard hierarchical methods by using 1-NN graphs to ensure balanced
cluster sizes.

## Usage

``` r
rena(
  bvec,
  mask,
  K = 100,
  connectivity = 26,
  max_iterations = 50,
  verbose = FALSE,
  exact_k = TRUE
)
```

## Arguments

- bvec:

  A `NeuroVec` instance supplying the 4D data to cluster.

- mask:

  A `NeuroVol` mask defining the voxels to include in clustering. If
  numeric, nonzero values define included voxels. If logical, TRUE
  values define included voxels.

- K:

  The number of clusters to find (default 100).

- connectivity:

  Neighborhood connectivity (6 or 26). Default 26. 6 = face neighbors
  only, 26 = face + edge + corner neighbors.

- max_iterations:

  Maximum number of recursion iterations (default 50). Algorithm stops
  when K clusters are reached or max_iterations is hit.

- verbose:

  Logical; whether to print progress messages. Default FALSE.

- exact_k:

  Logical; whether to use edge pruning to ensure exactly K clusters.
  Default TRUE. If FALSE, may produce slightly more or fewer than K
  clusters.

## Value

A `list` of class `rena_cluster_result` (inheriting from
`cluster_result`) with the following elements:

- clusvol:

  An instance of type
  [ClusteredNeuroVol](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.html).

- cluster:

  A vector of cluster indices equal to the number of voxels in the mask.

- centers:

  A matrix of cluster centers with each row representing the feature
  vector for a cluster.

- coord_centers:

  A matrix of spatial coordinates with each row corresponding to a
  cluster.

- n_clusters:

  Actual number of clusters found.

- method:

  Clustering method used ("rena").

- parameters:

  List of all parameters used.

- metadata:

  List containing iteration count and convergence information.

## Details

### Algorithm Description

ReNA recursively aggregates features using a 1-Nearest Neighbor (1-NN)
graph to ensure clusters remain roughly balanced in size and spatially
compact. The algorithm proceeds in iterations:

1.  **Connectivity Constraint**: Clustering respects the spatial
    topology graph G. Distance is only calculated between spatially
    adjacent voxels.

2.  **Similarity Calculation**: For all connected edges, compute squared
    Euclidean distance in feature space.

3.  **1-NN Subgraph**: For every voxel, identify its single nearest
    neighbor. This forms a directed subgraph Q. The 1-NN graph is less
    likely to "percolate" (create giant components) compared to k-NN
    graphs where k \>= 2.

4.  **Connected Components**: Extract weakly connected components of Q.
    These components become the clusters for the current iteration.

5.  **Reduction**: Features within each component are averaged to form a
    single "super-feature" for the next iteration. The graph G is
    contracted so edges between new clusters exist if any constituent
    voxels were connected.

6.  **Stopping Condition**: Repeat recursively until the number of
    clusters drops to the desired K. If exact_k = TRUE, edges are pruned
    by distance to ensure exactly K clusters.

### Performance Characteristics

- **Complexity**: O(N log N) per iteration where N = number of voxels

- **Memory**: O(N) for graph and assignments

- **Speed**: Generally faster than iterative methods (supervoxels)

- **Iterations**: Typically 5-15 iterations to reach K clusters from N
  voxels

### Advantages over Other Methods

- **No percolation**: Unlike single-linkage, doesn't form one giant
  cluster

- **Balanced sizes**: 1-NN graph encourages roughly equal cluster sizes

- **Topology-aware**: Respects spatial structure via connectivity graph

- **Deterministic**: No random initialization

- **Fast convergence**: Linear-time operations per iteration

### Parallelization Status

**Currently NOT parallelized.** ReNA uses a sequential recursive
algorithm that processes the graph in a hierarchical manner. Each
iteration depends on the previous iteration's clustering result.

#### Sequential Operations:

1.  **Distance Computation**: Pairwise distances for connected voxels
    (C++)

2.  **1-NN Finding**: Each voxel finds its nearest neighbor (C++)

3.  **Component Detection**: Union-Find algorithm for connected
    components (C++)

4.  **Feature Aggregation**: Mean pooling within components (C++)

5.  **Graph Contraction**: Building reduced graph for next iteration
    (C++)

#### Why Not Parallelized:

- **Sequential dependency**: Each iteration requires previous
  iteration's result

- **Fast per-iteration**: Linear-time operations make parallelization
  overhead high

- **Small iteration count**: Typically 5-15 iterations total

- **C++ optimization**: Core operations already optimized in C++

#### Performance Tips:

- **Reduce connectivity**: Use connectivity=6 instead of 26 for larger
  graphs

- **Use smaller K**: Fewer target clusters means fewer iterations

- **Pre-smooth data**: Reduces noise, improves cluster coherence

- **Alternative for parallelism**: Consider
  [`slice_msf()`](slice_msf.md) or `flash3d()`

## Note

Consider using [`cluster4d`](cluster4d.md) with `method = "rena"` for a
standardized interface across all clustering methods.

## References

Hoyos-Idrobo, A., Varoquaux, G., Kahn, J., & Thirion, B. (2019).
Recursive Nearest Agglomeration (ReNA): Fast clustering for
approximation of structured signals. *Pattern Recognition*, 94, 17-28.

## See also

[`cluster4d`](cluster4d.md) for unified interface,
[`supervoxels`](supervoxels.md) for iterative heat kernel clustering,
[`snic`](snic.md) for non-iterative priority queue clustering

## Examples

``` r
if (FALSE) { # \dontrun{
# Small synthetic example
library(neuroim2)
mask <- NeuroVol(array(1, c(10,10,10)), NeuroSpace(c(10,10,10)))
vec <- replicate(20,
                 NeuroVol(array(runif(10*10*10), c(10,10,10)),
                          NeuroSpace(c(10,10,10))),
                 simplify=FALSE)
vec <- do.call(concat, vec)

# Run ReNA clustering
rena_res <- rena(vec, mask, K=50, connectivity=26)
print(rena_res$n_clusters)

# With verbose output
rena_res <- rena(vec, mask, K=50, verbose=TRUE)
} # }
```
