# K-nearest-neighbor shrink

Replace each voxel by the mean of its k nearest neighbors in its local
spatial neighborhood.

## Usage

``` r
knn_shrink(bvec, mask, k = 5, connectivity = 27)
```

## Arguments

- bvec:

  A
  [`NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html)
  instance (the data).

- mask:

  A
  [`NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol-class.html)
  mask defining the voxels to include. If numeric, nonzero = included.

- k:

  The number of nearest neighbors to average over.

- connectivity:

  The number of spatial neighbors to include in the search around each
  voxel.

## Value

A `SparseNeuroVec` or similar object with the smoothed data.

## Examples

``` r
if (FALSE) { # \dontrun{
mask <- neuroim2::NeuroVol(array(1, c(20,20,20)), neuroim2::NeuroSpace(c(20,20,20)))
bvec <- replicate(10,
                  neuroim2::NeuroVol(array(runif(20*20*20), c(20,20,20)),
                           neuroim2::NeuroSpace(c(20,20,20))),
                  simplify=FALSE)
bvec <- do.call(neuroim2::concat, bvec)
sbvec <- knn_shrink(bvec, mask, k=3)
} # }
```
