# Estimate Optimal Number of SVD Components

Uses the elbow method or cumulative variance to suggest an appropriate
number of components for SVD compression.

## Usage

``` r
estimate_optimal_components(
  feature_mat,
  max_components = 50,
  method = c("variance", "elbow"),
  variance_target = 0.95
)
```

## Arguments

- feature_mat:

  Numeric matrix (N x T)

- max_components:

  Maximum components to test (default: 50)

- method:

  Either "elbow" (find elbow in scree plot) or "variance" (return number
  needed for 95% variance). Default: "variance".

- variance_target:

  Target variance for "variance" method (default: 0.95)

## Value

Integer; suggested number of components

## Examples

``` r
if (FALSE) { # \dontrun{
data <- matrix(rnorm(1000 * 300), 1000, 300)
n_opt <- estimate_optimal_components(data)
print(n_opt)
} # }
```
