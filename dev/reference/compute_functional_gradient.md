# Compute Functional Gradient from Feature Matrix

Calculates the functional gradient for each voxel as the average
dissimilarity to its spatial neighbors in feature space. Low gradient
values indicate stable, homogeneous regions (good seed locations); high
values indicate boundaries.

## Usage

``` r
compute_functional_gradient(feature_mat, coords, k_neighbors = 26)
```

## Arguments

- feature_mat:

  Numeric matrix (N x M) of features (voxels x dimensions).

- coords:

  Numeric matrix (N x 3) of spatial coordinates.

- k_neighbors:

  Integer; number of nearest spatial neighbors to use. Default: 26.

## Value

Numeric vector of length N containing gradient values for each voxel.

## Details

For each voxel i, the functional gradient is: \$\$grad(i) = \frac{1}{k}
\sum\_{j \in N_k(i)} (1 - f_i \cdot f_j)\$\$

where \\N_k(i)\\ are the k nearest spatial neighbors and \\f_i \cdot
f_j\\ is the dot product (cosine similarity for unit-length vectors).

## Examples

``` r
if (FALSE) { # \dontrun{
features <- matrix(rnorm(500 * 15), 500, 15)
features <- t(apply(features, 1, function(x) x / sqrt(sum(x^2))))
coords <- matrix(rnorm(500 * 3), 500, 3)

grad <- compute_functional_gradient(features, coords, k_neighbors = 26)
hist(grad, breaks = 50, main = "Functional Gradient Distribution")
} # }
```
