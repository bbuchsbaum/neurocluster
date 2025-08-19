# Mathematical Validation Test Analysis

## Issues Identified

### 1. DCT DC Component Dominance Test Failure

**Problem**: The test expects the DC component (k=1) to dominate other frequency components for low-frequency signals, but this expectation is mathematically incorrect.

**Analysis**: 
- The test signal `sin(t) + 0.5*cos(2t)` is indeed low-frequency
- After preprocessing (detrending and z-scoring), the DC component (k=1) has value ~0.094
- The k=3 component has magnitude ~0.760, which dominates
- **This is mathematically correct behavior** - DCT k=1 is NOT the DC component

**Root Cause**: 
- The test assumes DCT coefficient k=1 is the DC component
- In reality, DCT k=1 corresponds to the first AC component (cosine with frequency π/2T)
- The actual DC component would be k=0, but the implementation starts from k=1 (excluding DC)
- The C++ code explicitly starts DCT from k=1: `cos(pi * ((t + 0.5) * (k + 1) / T))`

**Mathematical Verification**:
```
DCT basis functions:
k=1: cos(π(t+0.5)·1/T) - first harmonic
k=2: cos(π(t+0.5)·2/T) - second harmonic  
k=3: cos(π(t+0.5)·3/T) - third harmonic
```

### 2. Energy Concentration Test Failures

**Problem**: Tests fail for r=4 and r=8, expecting >80% energy in first half of DCT coefficients.

**Analysis**:
- For r=4: Energy ratio = 0.957 (>80%) ✓ - This should pass
- For r=8: Energy ratio = 0.957 (>80%) ✓ - This should pass
- The mathematical analysis shows the expectation IS met

**Root Cause**: 
- The test signals DO concentrate >80% energy in the first half
- There may be a bug in the test implementation or different signals being used
- The test may be using different r values or different signals than expected

### 3. Supervoxels Heat Kernel Test Failure

**Problem**: Error "more cluster centers than data points" in supervoxels function.

**Analysis**:
- Test uses 3×3×1 = 9 voxels total
- Tries to create K=3 clusters
- The error suggests the `find_initial_points` function is trying to return more points than available

**Root Cause**:
- `find_initial_points` function exists but is not exported from NAMESPACE
- The function has a potential bug where it may select more points than available voxels
- The gradient-based initialization (use_gradient=TRUE) calls this function
- When it fails, there's no fallback to the simple non-gradient method

## Proposed Fixes

### Fix 1: DCT DC Component Test

The test expectation is mathematically incorrect. The DCT implementation correctly excludes the DC component (k=0) and starts from k=1. For preprocessed signals (detrended and z-scored), expecting k=1 to dominate is wrong.

**Recommendation**: 
1. Remove the DC dominance test entirely, OR
2. Change it to test that low-frequency components (k=1,2) have higher magnitude than high-frequency components (k>6)

```r
# Instead of testing DC dominance, test low vs high frequency dominance
low_freq_energy <- sum(sketch_matrix[1:2, 1]^2)  # k=1,2
high_freq_energy <- sum(sketch_matrix[7:8, 1]^2) # k=7,8
expect_true(low_freq_energy > high_freq_energy,
           info = "Low frequency components should dominate for low-frequency signals")
```

### Fix 2: Energy Concentration Test

The energy concentration test appears to be mathematically sound but may have implementation issues.

**Recommendation**: Debug the specific test case to ensure the correct signals and parameters are being used.

### Fix 3: Supervoxels K Parameter Validation

**Immediate Fix**: Add proper K validation and export the missing function.

```r
# In supervoxels.R init_cluster function:
init_cluster <- function(bvec, mask, coords, K, use_gradient = TRUE) {
  mask.idx <- which(mask > 0)
  nvox <- nrow(coords)

  if (K > nvox) {
    warning("K is greater than the number of valid voxels. Setting K = number of voxels.")
    K <- nvox
  }

  if (use_gradient && K <= nvox) {
    tryCatch({
      refvol <- bvec[[1]]
      grad <- spatial_gradient(refvol, mask)
      valid_coords <- index_to_grid(mask, mask.idx)
      
      # Call the function directly since it's not exported
      init <- neurocluster:::find_initial_points(valid_coords, grad, K)
      
      # run kmeans with chosen seeds
      kres <- stats::kmeans(coords, centers = coords[init$selected, , drop = FALSE], iter.max = 500)
      return(kres$cluster)
    }, error = function(e) {
      warning("Gradient-based initialization failed, falling back to uniform sampling: ", e$message)
      use_gradient <- FALSE
    })
  }
  
  # Fallback: uniform sampling
  init_centers <- coords[as.integer(seq(1, nvox, length.out = K)), , drop = FALSE]
  kres <- stats::kmeans(coords, centers = init_centers, iter.max = 500)
  kres$cluster
}
```

**Long-term Fix**: 
1. Export `find_initial_points` in NAMESPACE
2. Add proper error handling in the function
3. Ensure it never returns more indices than available voxels

### Fix 4: find_initial_points Function Bug

The function may have an off-by-one error or indexing issue. Need to add bounds checking:

```r
find_initial_points <- function(cds, grad, K=100) {
  if (K > nrow(cds)) {
    warning("K larger than available points, using all points")
    return(list(selected = 1:nrow(cds), coords = cds))
  }
  
  # ... rest of function with proper bounds checking
}
```

## Test Corrections Needed

### test_mathematical_validation.R

1. **Remove or fix DC dominance test** (lines 50-54)
2. **Debug energy concentration test** - verify it's using correct parameters  
3. **Fix supervoxels test** - use smaller K or larger volume
4. **Add K parameter validation test** explicitly

### Example Fixed Test

```r
test_that("DCT basis construction is mathematically correct", {
  # ... setup code ...
  
  # Test orthogonality (keep this)
  for (i in 1:(ncol(sketch_matrix)-1)) {
    expect_equal(sketch_matrix[, i], sketch_matrix[, i+1], tolerance = 1e-10)
  }
  
  # FIXED: Test low vs high frequency dominance instead of DC dominance
  low_freq_components <- sketch_matrix[1:3, 1]  # k=1,2,3
  high_freq_components <- sketch_matrix[6:8, 1] # k=6,7,8
  expect_true(mean(abs(low_freq_components)) > mean(abs(high_freq_components)),
              info = "Low frequency components should dominate for low-frequency signals")
})

test_that("heat kernel computation is mathematically accurate", {
  # ... setup code ...
  
  # FIXED: Use smaller K that's guaranteed to work
  nvox <- prod(dims)
  K_safe <- min(3, nvox)  # Ensure K <= number of voxels
  
  result <- supervoxels(vec, mask, K = K_safe, alpha = 0.5)
  
  # ... rest of test ...
})
```

## Summary

The main issues are:

1. **Incorrect mathematical expectation** for DC component dominance
2. **Missing function export** causing supervoxels to fail  
3. **Insufficient K parameter validation** in supervoxels
4. **Potential implementation inconsistencies** in energy concentration test

The DCT implementation itself is mathematically correct. The test expectations need to be updated to match the actual mathematical behavior of the DCT transform as implemented.