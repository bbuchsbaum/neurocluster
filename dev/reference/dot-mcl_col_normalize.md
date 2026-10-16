# Sparse Markov Clustering (MCL) Backend

Internal sparse MCL implementation for graph clustering of masked
voxels. The implementation keeps iterates sparse via per-column pruning
to remain practical on neuroimaging graphs.

## Usage

``` r
.mcl_col_normalize(mat)
```
