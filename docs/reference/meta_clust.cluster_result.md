# Meta Clustering for Cluster Results

The meta_clust function performs meta clustering on a given clustering
result by applying hierarchical clustering on the cluster centers.

## Usage

``` r
# S3 method for class 'cluster_result'
meta_clust(x, cuts = NULL, ...)
```

## Arguments

- x:

  A clustering result, typically an object of class `"cluster_result"`.

- cuts:

  Integer vector specifying the number of cluster cuts to consider.
  Default is `NULL`, which generates cuts at 2, 5, and 10 clusters (or
  fewer depending on the number of input clusters).

- ...:

  Additional arguments:

  algo

  :   A character string indicating the clustering algorithm to use.
      Default is "hclust" (hierarchical clustering).

  dist_method

  :   Character string: "correlation" (default) or "euclidean" for
      distance calculation between cluster centers.

  hclust_method

  :   A character string specifying the agglomeration method to use for
      hierarchical clustering. Default is "ward.D".

## Value

A list containing:

- cvols:

  A list of
  [`ClusteredNeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVol-class.html)
  instances.

- cuts:

  The number of cluster cuts.

- cutmat:

  A matrix representing the cluster assignments for each cut.

- hclus:

  The hierarchical clustering result.

- original_result:

  The original clustering result (optional reference).

## Details

### Parallelization Status

**Currently NOT parallelized.** Meta-clustering performs hierarchical
clustering on cluster centers using sequential R functions.

#### Sequential Operations:

1.  **Distance Matrix Computation**:

    - Computes correlation distance (1 - cor) between all cluster
      centers

    - O(K²) pairwise correlations where K = number of input clusters

    - Sequential double loop in base R

2.  **Hierarchical Clustering**:
    [`hclust()`](https://rdrr.io/r/stats/hclust.html)

    - Standard hierarchical agglomeration (Ward.D by default)

    - Sequential merging of closest clusters

    - O(K² log K) complexity

3.  **Tree Cutting**: [`cutree()`](https://rdrr.io/r/stats/cutree.html)

    - Cuts dendrogram at multiple heights

    - Creates nested cluster assignments

    - Linear in number of clusters

#### Why Not Parallelized:

- **Small scale**: Usually operates on hundreds of clusters, not
  thousands of voxels

- **R built-ins**: Uses standard
  [`hclust()`](https://rdrr.io/r/stats/hclust.html) which is sequential

- **Fast enough**: Typically completes in seconds even for large K

- **Memory efficient**: Only stores K×K distance matrix

#### Potential for Parallelization:

- Distance matrix computation could use parallel correlation

- Some hierarchical clustering variants support parallelization

- Multiple cuts could be computed in parallel

#### Performance Characteristics:

- **Input size**: K clusters from initial clustering (typically
  100-1000)

- **Complexity**: O(K²) for distances, O(K² log K) for clustering

- **Memory**: O(K²) for distance matrix

- **Speed**: Usually \< 1 second for K \< 500

#### Performance Tips:

- **Pre-reduce clusters**: Start with fewer initial clusters if
  meta-clustering is slow

- **Use appropriate linkage**: Ward.D is slower but often gives better
  results

- **Limit cuts**: Fewer cut levels = faster processing

- **Alternative**: Use consensus clustering
  ([`merge_clus()`](merge_clus.md)) for different approach

#### Use Cases:

Meta-clustering is useful for:

- Creating multi-resolution parcellations

- Hierarchical organization of functional regions

- Reducing large numbers of clusters to interpretable groups

- Finding stable cluster boundaries across scales

## See also

[`hclust`](https://rdrr.io/r/stats/hclust.html),
[`cutree`](https://rdrr.io/r/stats/cutree.html)
