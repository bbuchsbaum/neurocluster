# End-to-end: From NIfTI to clusters

## I/O and setup

``` r
# Example with synthetic data (safe to run)
mask <- NeuroVol(array(1, c(6,6,6)), NeuroSpace(c(6,6,6)))
vec  <- NeuroVec(array(rnorm(6*6*6*12), c(6,6,6,12)), NeuroSpace(c(6,6,6,12)))

# Real data (uncomment and supply paths)
# mask <- read_vol("brain_mask.nii")
# vec  <- read_vec("fmri_data.nii.gz")
```

## Cluster

``` r
res <- cluster4d(vec, mask, n_clusters = 5, method = "slic", max_iterations = 5)
```

## Inspect and visualize

``` r
print(res)
summary(res)
plot(res)
```

## Export NIfTI

``` r
writeVol(res$clusvol, "clusters.nii.gz")
```

## Next steps

- Tune parameters: articles/choose-parameters.html#suggestions
- Compare methods: articles/compare-methods.html#run-methods
- Validate results: articles/validate-compare.html#checks
