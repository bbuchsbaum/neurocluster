#!/usr/bin/env Rscript
# Debug energy concentration test

library(neurocluster)
library(neuroim2)

debug_energy_concentration <- function() {
  cat("=== Debug Energy Concentration Test ===\n")
  
  dims <- c(4, 4, 1)
  mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
  nvox <- prod(dims)
  ntime <- 32
  
  # Create test signal with known frequency content
  t_seq <- seq(0, 2*pi, length.out = ntime)
  
  # Very low-frequency signal - closer to DC with gentle variation
  low_freq_signal <- 1 + 0.3*sin(0.5*t_seq) + 0.1*cos(0.8*t_seq)
  
  cat("Testing very low-frequency signal: 1 + 0.3*sin(0.5*t) + 0.1*cos(0.8*t)\n")
  cat("Time points:", ntime, "\n")
  cat("Voxels:", nvox, "\n")
  
  # Test with different DCT ranks
  r_values <- c(4, 8, 16)
  
  for (r in r_values) {
    cat(sprintf("\n--- Testing r = %d ---\n", r))
    
    # Test low-frequency signal
    ts_data_low <- matrix(rep(low_freq_signal, nvox), nrow = nvox, ncol = ntime, byrow = TRUE)
    
    vec_list_low <- lapply(1:ntime, function(t) {
      vol_data <- array(0, dims)
      vol_data[mask > 0] <- ts_data_low[, t]
      NeuroVol(vol_data, NeuroSpace(dims))
    })
    vec_low <- do.call(concat, vec_list_low)
    
    result_low <- slice_msf_single(vec_low, mask, r = r, k = 0.5, min_size = 1)
    
    # For low-frequency signals, DCT should capture most energy in first few coefficients
    sketches_low <- result_low$sketch[, 1]  # First voxel's sketch
    
    cat("DCT coefficients:\n")
    for (k in 1:r) {
      cat(sprintf("  k=%d: %.6f\n", k, sketches_low[k]))
    }
    
    energy_first_half <- sum(sketches_low[1:(r/2)]^2)
    total_energy <- sum(sketches_low^2)
    energy_ratio <- energy_first_half / total_energy
    
    cat(sprintf("Energy analysis:\n"))
    cat(sprintf("  First half coefficients (k=1 to %d): %s\n", r/2, paste(1:(r/2), collapse=",")))
    cat(sprintf("  Energy in first half: %.6f\n", energy_first_half))
    cat(sprintf("  Total energy: %.6f\n", total_energy))
    cat(sprintf("  Energy ratio: %.6f\n", energy_ratio))
    cat(sprintf("  Passes 80%% threshold: %s\n", energy_ratio > 0.8))
    
    # Detailed breakdown
    cat("Individual coefficient energies:\n")
    for (k in 1:r) {
      coeff_energy <- sketches_low[k]^2
      cat(sprintf("  k=%d: %.6f (%.1f%%)\n", k, coeff_energy, 100*coeff_energy/total_energy))
    }
  }
}

# Run the debug
debug_energy_concentration()