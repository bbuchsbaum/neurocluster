# commute_cluster isSymmetric Error Fix

## Problem

The [`commute_cluster()`](reference/commute_cluster.md) function was
failing with the following error:

    Error in `UseMethod("isSymmetric")`:
      no applicable method for 'isSymmetric' applied to an object of class
      "c('dgCMatrix', 'CsparseMatrix', 'dsparseMatrix', 'generalMatrix',
         'dMatrix', 'sparseMatrix', 'Matrix')"

## Root Cause

The issue occurred at
[R/commute_cluster.R:274](R/commute_cluster.R#L274):

``` r
sym_ok <- if (inherits(W, "Matrix")) Matrix::isSymmetric(W, checkDN = FALSE) else isSymmetric(W)
```

While the code correctly checked `inherits(W, "Matrix")` and attempted
to call
[`Matrix::isSymmetric()`](https://rdrr.io/pkg/Matrix/man/isSymmetric-methods.html),
the function was not available in the package namespace because it
wasn’t explicitly imported.

The imports section only had:

``` r
#' @importFrom Matrix t
```

This meant that when
[`Matrix::isSymmetric()`](https://rdrr.io/pkg/Matrix/man/isSymmetric-methods.html)
was called, R couldn’t find the function, causing the error to fall
through to
[`base::isSymmetric()`](https://rdrr.io/r/base/isSymmetric.html), which
doesn’t have a method for sparse matrix classes.

## Solution

Added `isSymmetric` to the Matrix imports in
[R/commute_cluster.R:195](R/commute_cluster.R#L195):

``` r
#' @importFrom Matrix t isSymmetric
```

Then regenerated package documentation to update NAMESPACE:

``` bash
Rscript -e "devtools::document()"
R CMD INSTALL .
```

## Verification

After the fix: - ✅ [`commute_cluster()`](reference/commute_cluster.md)
runs successfully - ✅ Test `test_clustering.R:64` now passes - ✅
`isSymmetric` properly imported in NAMESPACE (line 55) - ✅ Sparse
matrix symmetry checking works correctly

### Test Results

``` r
library(neurocluster)
mask <- NeuroVol(array(1, c(5,5,5)), NeuroSpace(c(5,5,5)))
vec <- NeuroVec(array(rnorm(5*5*5*10), c(5,5,5,10)), NeuroSpace(c(5,5,5,10)))

result <- commute_cluster(vec, mask, K=10, ncomp=5)
# ✓ Completes successfully
# ✓ Returns 10 clusters
# ✓ Returns "commute_time_cluster_result" object
```

## Files Modified

1.  **R/commute_cluster.R** (line 195)
    - Added `isSymmetric` to `@importFrom Matrix` directive
2.  **NAMESPACE** (auto-generated, line 55)
    - Added `importFrom(Matrix,isSymmetric)`

## Lesson Learned

When using functions from the Matrix package with the
`Package::function()` syntax in roxygen2-documented R packages, the
function must still be explicitly imported via `@importFrom` directives.
Simply using the namespace prefix (`Matrix::`) is not sufficient for
package exports without the corresponding import declaration.
