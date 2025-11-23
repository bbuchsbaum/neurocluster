# Tesselate a Mask Volume into K Clusters using K-means

This function tesselates a given mask volume into K clusters using
k-means clustering applied to spatial coordinates. It returns a
clustered mask volume object.

## Usage

``` r
tesselate(mask, K = 100)
```

## Arguments

- mask:

  A `NeuroVol` object representing the mask volume.

- K:

  An integer value specifying the number of clusters (default: 100).

  If `K` exceeds the number of nonzero voxels, a warning is issued and
  `K` is set to the number of nonzero voxels.

## Value

An instance of `ClusteredNeuroVol` representing the clustered mask
volume.

## Examples

``` r
# Assuming you have a NeuroVol object 'mask' and you want to create 150 clusters
clustered_volume <- tesselate(mask, K = 150)
#> Error: object 'mask' not found
```
