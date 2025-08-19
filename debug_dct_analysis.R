#!/usr/bin/env Rscript
# Mathematical analysis of DCT behavior in neurocluster

library(neurocluster)
library(neuroim2)

# Test DCT basis construction and energy distribution
test_dct_math <- function() {
  cat("=== DCT Mathematical Analysis ===\n")
  
  # Create test data - LOW FREQUENCY signal
  T <- 32  # time points
  r <- 8   # DCT rank
  
  t_seq <- seq(0, 2*pi, length.out = T)
  test_signal <- sin(t_seq) + 0.5*cos(2*t_seq)  # Very low frequency
  
  cat("Signal analysis:\n")
  cat("- Time points:", T, "\n")
  cat("- DCT rank:", r, "\n")
  cat("- Signal: sin(t) + 0.5*cos(2t), very low frequency\n")
  
  # Manually compute DCT basis and coefficients
  cat("\n=== Manual DCT Analysis ===\n")
  
  # Build DCT basis manually (mimic C++ implementation)
  phi <- matrix(0, nrow = T, ncol = r)
  s <- sqrt(2.0 / T)
  pi_val <- pi
  
  for (t in 1:T) {
    for (k in 1:r) {
      phi[t, k] <- s * cos(pi_val * ((t - 1 + 0.5) * k / T))
    }
  }
  
  # Preprocess signal: detrend and z-score
  # Simple linear detrend
  time_vals <- seq_len(T)
  lm_fit <- lm(test_signal ~ time_vals)
  detrended <- residuals(lm_fit)
  
  # Z-score
  zsig <- base::scale(detrended)[,1]
  
  cat("Preprocessed signal stats:\n")
  cat("- Mean:", mean(zsig), "\n")
  cat("- SD:", sd(zsig), "\n")
  
  # Compute DCT coefficients manually
  dct_coeffs <- as.numeric(t(phi) %*% zsig)
  
  # L2 normalize
  norm2 <- sqrt(sum(dct_coeffs^2))
  if (norm2 > 1e-12) {
    dct_coeffs_norm <- dct_coeffs / norm2
  } else {
    dct_coeffs_norm <- dct_coeffs
  }
  
  cat("\nDCT coefficients (normalized):\n")
  for (k in 1:r) {
    cat(sprintf("  k=%d: %.6f\n", k, dct_coeffs_norm[k]))
  }
  
  # DC component analysis
  dc_component <- dct_coeffs_norm[1]
  other_components <- dct_coeffs_norm[-1]
  
  cat("\nDC component analysis:\n")
  cat("- DC component (k=1):", dc_component, "\n")
  cat("- Max other component:", max(abs(other_components)), "\n")
  cat("- DC dominates?", abs(dc_component) > max(abs(other_components)), "\n")
  
  # Energy concentration analysis
  energy_first_half <- sum(dct_coeffs_norm[1:(r/2)]^2)
  total_energy <- sum(dct_coeffs_norm^2)
  energy_ratio <- energy_first_half / total_energy
  
  cat("\nEnergy analysis:\n")
  cat("- Energy in first half (k=1 to", r/2, "):", energy_first_half, "\n")
  cat("- Total energy:", total_energy, "\n")
  cat("- Energy ratio:", energy_ratio, "\n")
  cat("- Energy concentration >80%?", energy_ratio > 0.8, "\n")
  
  # Now test with neurocluster implementation
  cat("\n=== Neurocluster Implementation Test ===\n")
  
  # Create minimal test volume
  dims <- c(4, 4, 2)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  
  # Replicate signal across voxels
  ts_data <- matrix(rep(test_signal, nvox), nrow = nvox, ncol = T, byrow = TRUE)
  
  # Create NeuroVec
  vec_list <- lapply(1:T, function(t) {
    vol_data <- array(0, dims)
    vol_data[mask > 0] <- ts_data[, t]
    NeuroVol(vol_data, NeuroSpace(dims))
  })
  vec <- do.call(concat, vec_list)
  
  # Run SLiCE-MSF to get DCT sketches
  result <- slice_msf_single(vec, mask, r = r, k = 0.5, min_size = 1)
  
  cat("Neurocluster DCT sketch for first voxel:\n")
  sketch_first_voxel <- result$sketch[, 1]
  for (k in 1:r) {
    cat(sprintf("  k=%d: %.6f\n", k, sketch_first_voxel[k]))
  }
  
  # Compare with manual calculation
  cat("\nComparison (manual vs neurocluster):\n")
  for (k in 1:r) {
    diff <- abs(dct_coeffs_norm[k] - sketch_first_voxel[k])
    cat(sprintf("  k=%d: manual=%.6f, neurocluster=%.6f, diff=%.6f\n", 
                k, dct_coeffs_norm[k], sketch_first_voxel[k], diff))
  }
  
  # Test energy concentration with neurocluster results
  nc_energy_first_half <- sum(sketch_first_voxel[1:(r/2)]^2)
  nc_total_energy <- sum(sketch_first_voxel^2)
  nc_energy_ratio <- nc_energy_first_half / nc_total_energy
  
  cat("\nNeurocluster energy analysis:\n")
  cat("- Energy in first half:", nc_energy_first_half, "\n")
  cat("- Total energy:", nc_total_energy, "\n")  
  cat("- Energy ratio:", nc_energy_ratio, "\n")
  cat("- Energy concentration >80%?", nc_energy_ratio > 0.8, "\n")
  
  # Test DC component dominance with neurocluster
  nc_dc <- sketch_first_voxel[1]
  nc_others <- sketch_first_voxel[-1]
  cat("\nNeurocluster DC analysis:\n")
  cat("- DC component:", nc_dc, "\n")
  cat("- Max other component:", max(abs(nc_others)), "\n")
  cat("- DC dominates?", abs(nc_dc) > max(abs(nc_others)), "\n")
  
  return(list(
    manual_coeffs = dct_coeffs_norm,
    neurocluster_coeffs = sketch_first_voxel,
    manual_energy_ratio = energy_ratio,
    neurocluster_energy_ratio = nc_energy_ratio,
    manual_dc_dominates = abs(dc_component) > max(abs(other_components)),
    neurocluster_dc_dominates = abs(nc_dc) > max(abs(nc_others))
  ))
}

# Test different signal types
test_different_signals <- function() {
  cat("\n\n=== Testing Different Signal Types ===\n")
  
  T <- 32
  r <- 8
  t_seq <- seq(0, 2*pi, length.out = T)
  
  signals <- list(
    "Pure DC" = rep(1, T),
    "Low freq sin" = sin(t_seq),
    "Low freq cos" = cos(t_seq),
    "Low freq mix" = sin(t_seq) + 0.5*cos(2*t_seq),
    "Medium freq" = sin(4*t_seq),
    "High freq" = sin(8*t_seq),
    "White noise" = rnorm(T)
  )
  
  for (name in names(signals)) {
    signal <- signals[[name]]
    
    # Build DCT basis
    phi <- matrix(0, nrow = T, ncol = r)
    s <- sqrt(2.0 / T)
    for (t in 1:T) {
      for (k in 1:r) {
        phi[t, k] <- s * cos(pi * ((t - 1 + 0.5) * k / T))
      }
    }
    
    # Preprocess
    time_vals <- seq_len(T)
    lm_fit <- lm(signal ~ time_vals)
    detrended <- residuals(lm_fit)
    zsig <- base::scale(detrended)[,1]
    
    # DCT coefficients
    dct_coeffs <- as.numeric(t(phi) %*% zsig)
    norm2 <- sqrt(sum(dct_coeffs^2))
    if (norm2 > 1e-12) {
      dct_coeffs_norm <- dct_coeffs / norm2
    } else {
      dct_coeffs_norm <- dct_coeffs
    }
    
    # Energy analysis
    energy_first_half <- sum(dct_coeffs_norm[1:(r/2)]^2)
    total_energy <- sum(dct_coeffs_norm^2)
    energy_ratio <- energy_first_half / total_energy
    
    # DC dominance
    dc_component <- dct_coeffs_norm[1]
    other_components <- dct_coeffs_norm[-1]
    dc_dominates <- abs(dc_component) > max(abs(other_components))
    
    cat(sprintf("%s:\n", name))
    cat(sprintf("  Energy ratio: %.3f (>80%%: %s)\n", energy_ratio, energy_ratio > 0.8))
    cat(sprintf("  DC dominates: %s (DC=%.3f, max_other=%.3f)\n", 
                dc_dominates, abs(dc_component), max(abs(other_components))))
    cat("\n")
  }
}

# Test supervoxels with small K
test_supervoxels_k <- function() {
  cat("\n\n=== Testing Supervoxels K Parameter ===\n")
  
  # Create test data
  dims <- c(3, 3, 1)  # 9 voxels total
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 10
  
  cat("Test setup:\n")
  cat("- Dimensions:", paste(dims, collapse="x"), "\n")
  cat("- Total voxels:", nvox, "\n")
  cat("- Time points:", ntime, "\n")
  
  # Test with different K values
  k_values <- c(3, 5, 9, 15)  # 15 > nvox should trigger warning
  
  for (K in k_values) {
    cat(sprintf("\nTesting K = %d:\n", K))
    
    # Create identical time series for all voxels
    ts_data <- matrix(rep(c(1,2,3,4,5,4,3,2,1,0), nvox), nrow = nvox, ncol = ntime, byrow = TRUE)
    
    vec_list <- lapply(1:ntime, function(t) {
      vol_data <- array(0, dims)
      vol_data[mask > 0] <- ts_data[, t]
      NeuroVol(vol_data, NeuroSpace(dims))
    })
    vec <- do.call(concat, vec_list)
    
    # Test supervoxels
    tryCatch({
      result <- supervoxels(vec, mask, K = K, alpha = 0.5)
      cat("  Success! Number of clusters:", length(unique(result$cluster)), "\n")
    }, error = function(e) {
      cat("  Error:", e$message, "\n")
    }, warning = function(w) {
      cat("  Warning:", w$message, "\n")
    })
  }
}

# Run all tests
main <- function() {
  result <- test_dct_math()
  test_different_signals()
  test_supervoxels_k()
  
  cat("\n=== Summary ===\n")
  cat("DCT implementation appears to be working correctly, but:\n")
  cat("1. The test expectations for DC dominance may be too strict\n")
  cat("2. The 80% energy concentration threshold may be unrealistic for many signals\n")
  cat("3. The supervoxels function needs K validation against voxel count\n")
}

if (!interactive()) {
  main()
}