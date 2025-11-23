# Visualize and export results

## Plot

``` r
plot(res)                 # default: axial/sagittal/coronal
plot(res, slice = c(2,2,2), view = "axial")
```

## Export

``` r
writeVol(res$clusvol, "clusters.nii.gz")
```

## Single cluster

``` r
cl_id <- 1
vol <- NeuroVol(array(0, dim(mask)), space(mask))
vol[res$cluster == cl_id] <- 1
writeVol(vol, sprintf("cluster_%d.nii.gz", cl_id))
```
