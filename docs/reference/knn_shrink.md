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
mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
#> Error in NeuroVol(array(1, c(20, 20, 20)), NeuroSpace(c(20, 20, 20))): could not find function "NeuroVol"
bvec <- replicate(10,
                  NeuroVol(array(runif(20*20*20), c(20,20,20)),
                           NeuroSpace(c(20,20,20))),
                  simplify=FALSE)
#> Error in NeuroVol(array(runif(20 * 20 * 20), c(20, 20, 20)), NeuroSpace(c(20,     20, 20))): could not find function "NeuroVol"
bvec <- do.call(concat, bvec)
#> Error: object 'bvec' not found

sbvec <- knn_shrink(bvec, mask, k=3)
#> Error: object 'bvec' not found
```
