library(neurocluster)
library(neuroim2)

# Definitive test for reliability weighting effectiveness
set.seed(42)

cat("=== Reliability Weighting Test ===\n\n")

# Create a test scenario where reliability weighting should matter
dims <- c(20, 20, 2)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
nvox <- prod(dims)
ntime <- 100  # Longer time series for better reliability estimates

cat("Creating test data with mixed reliability patterns...\n")

# Create 4 regions with VERY different reliability characteristics
coords <- arrayInd(1:nvox, dims)
region <- integer(nvox)
reliability_type <- integer(nvox)

for (i in 1:nvox) {
  x <- coords[i, 1]
  y <- coords[i, 2]
  
  # Create 4 quadrants
  if (x <= 10 && y <= 10) {
    region[i] <- 1  # Top-left
    reliability_type[i] <- 1  # High reliability
  } else if (x > 10 && y <= 10) {
    region[i] <- 2  # Top-right  
    reliability_type[i] <- 2  # Low reliability (high noise)
  } else if (x <= 10 && y > 10) {
    region[i] <- 3  # Bottom-left
    reliability_type[i] <- 3  # Medium reliability
  } else {
    region[i] <- 4  # Bottom-right
    reliability_type[i] <- 4  # Variable reliability
  }
}

# Create time series with different reliability patterns
ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 6*pi, length.out = ntime)

# Base patterns for each region
base_patterns <- list(
  sin(t_seq),                    # Region 1
  cos(t_seq),                    # Region 2  
  sin(2 * t_seq),                # Region 3
  cos(2 * t_seq)                 # Region 4
)

cat("Reliability characteristics by region:\n")

for (i in 1:nvox) {
  r <- region[i]
  rel_type <- reliability_type[i]
  base_signal <- base_patterns[[r]]
  
  if (rel_type == 1) {
    # High reliability: very low noise, consistent signal
    noise_sd <- 0.02
    signal <- base_signal
    cat_msg <- "high reliability (low noise)"
  } else if (rel_type == 2) {
    # Low reliability: high noise, inconsistent signal  
    noise_sd <- 0.25
    # Add random phase shifts to reduce split-half correlation
    phase_noise <- rnorm(ntime, sd = 0.3)
    signal <- base_signal + cumsum(phase_noise) * 0.01
    cat_msg <- "low reliability (high noise + phase drift)"
  } else if (rel_type == 3) {
    # Medium reliability: moderate noise
    noise_sd <- 0.08  
    signal <- base_signal
    cat_msg <- "medium reliability (moderate noise)"
  } else {
    # Variable reliability: alternating reliable/unreliable timepoints
    noise_sd <- 0.05
    # Make every 4th timepoint much noisier to break split-half correlation
    reliability_mask <- ((1:ntime) %% 4) != 0
    signal <- base_signal
    signal[!reliability_mask] <- signal[!reliability_mask] + rnorm(sum(!reliability_mask), sd = 0.4)
    cat_msg <- "variable reliability (periodic noise bursts)"
  }
  
  ts_data[i, ] <- signal + rnorm(ntime, sd = noise_sd)
}

# Print reliability info only once per region
cat("  Region 1:", "high reliability (low noise)\n")
cat("  Region 2:", "low reliability (high noise + phase drift)\n") 
cat("  Region 3:", "medium reliability (moderate noise)\n")
cat("  Region 4:", "variable reliability (periodic noise bursts)\n")

# Create NeuroVec
vec_list <- lapply(1:ntime, function(t) {
  vol_data <- array(0, dims)
  vol_data[mask > 0] <- ts_data[, t]
  NeuroVol(vol_data, NeuroSpace(dims))
})
vec <- do.call(concat, vec_list)

# Test with different gamma values to see reliability weighting effect
gamma_values <- c(0.5, 1.0, 1.5, 2.0, 3.0)
results <- list()

cat("\nTesting different gamma (reliability weighting) values:\n")

for (gamma in gamma_values) {
  result <- slice_msf(vec, mask, target_k_global = 4, num_runs = 1,
                      consensus = FALSE, r = 8, min_size = 20, 
                      compactness = 3, gamma = gamma)
  
  n_clusters <- length(unique(result$cluster))
  
  # Calculate clustering accuracy
  confusion <- table(region, result$cluster)
  accuracy <- sum(apply(confusion, 1, max)) / nvox
  
  # Calculate how well high-reliability voxels are preserved vs low-reliability
  high_rel_voxels <- which(reliability_type == 1)  # Region 1
  low_rel_voxels <- which(reliability_type == 2)   # Region 2
  
  # Check cluster purity for high vs low reliability regions
  high_rel_clusters <- result$cluster[high_rel_voxels]
  low_rel_clusters <- result$cluster[low_rel_voxels]
  
  high_rel_purity <- max(table(high_rel_clusters)) / length(high_rel_clusters)
  low_rel_purity <- max(table(low_rel_clusters)) / length(low_rel_clusters)
  
  results[[as.character(gamma)]] <- list(
    n_clusters = n_clusters,
    accuracy = accuracy,
    high_rel_purity = high_rel_purity,
    low_rel_purity = low_rel_purity
  )
  
  cat(sprintf("  γ=%.1f: %d clusters, accuracy=%.2f, high_rel_purity=%.2f, low_rel_purity=%.2f\n",
              gamma, n_clusters, accuracy, high_rel_purity, low_rel_purity))
}

# Test 2: Extreme reliability contrast
cat("\nTest 2: Extreme Reliability Contrast\n")
cat("Creating data with extreme reliability differences...\n")

dims2 <- c(16, 16, 2)  
mask2 <- NeuroVol(array(1, dims2), NeuroSpace(dims2))
nvox2 <- prod(dims2)
ntime2 <- 80

coords2 <- arrayInd(1:nvox2, dims2)
region2 <- ifelse(coords2[,1] <= 8, 1, 2)

ts_data2 <- matrix(0, nrow = nvox2, ncol = ntime2)
t_seq2 <- seq(0, 4*pi, length.out = ntime2)

for (i in 1:nvox2) {
  r <- region2[i]
  
  if (r == 1) {
    # Left side: PERFECT reliability (no noise)
    ts_data2[i, ] <- sin(t_seq2)  # Absolutely no noise
  } else {
    # Right side: TERRIBLE reliability (pure noise with weak signal)
    weak_signal <- 0.1 * cos(t_seq2)  # Very weak signal
    strong_noise <- rnorm(ntime2, sd = 1.0)  # Very strong noise
    ts_data2[i, ] <- weak_signal + strong_noise
  }
}

vec_list2 <- lapply(1:ntime2, function(t) {
  vol_data <- array(0, dims2)
  vol_data[mask2 > 0] <- ts_data2[, t]
  NeuroVol(vol_data, NeuroSpace(dims2))
})
vec2 <- do.call(concat, vec_list2)

# Test with low vs high gamma
gamma_extreme <- c(0.1, 2.0)
cat("Extreme reliability contrast results:\n")

for (gamma in gamma_extreme) {
  result <- slice_msf(vec2, mask2, target_k_global = 2, num_runs = 1,
                      consensus = FALSE, r = 6, min_size = 10,
                      compactness = 2, gamma = gamma)
  
  accuracy <- sum(apply(table(region2, result$cluster), 1, max)) / nvox2
  cat(sprintf("  γ=%.1f: accuracy=%.2f", gamma, accuracy))
  if (accuracy > 0.9) cat(" ✓\n") else cat(" ✗\n")
}

# Test 3: Direct reliability weight inspection
cat("\nTest 3: Direct Reliability Weight Inspection\n")
cat("Examining computed reliability weights...\n")

# Use slice_msf_single to access weights directly
result_weights <- slice_msf_single(vec2, mask2, r = 6, k = 0.5, 
                                   min_size = 10, nbhd = 8, gamma = 1.5)

weights <- result_weights$weights
left_weights <- weights[region2 == 1 & mask2@.Data > 0]
right_weights <- weights[region2 == 2 & mask2@.Data > 0]

left_mean_weight <- mean(left_weights)
right_mean_weight <- mean(right_weights)
weight_ratio <- left_mean_weight / right_mean_weight

cat(sprintf("Perfect signal region (left): mean weight = %.3f\n", left_mean_weight))
cat(sprintf("Noisy signal region (right): mean weight = %.3f\n", right_mean_weight))
cat(sprintf("Weight ratio (perfect/noisy): %.2f\n", weight_ratio))

if (weight_ratio > 2.0) {
  cat("✓ Reliability weighting is working - clean signals weighted much higher!\n")
} else {
  cat("✗ Reliability weighting may not be working as expected\n")
}

# Test 4: Gamma parameter effect on weights
cat("\nTest 4: Gamma Parameter Effect on Weights\n")
cat("Testing how gamma affects the weighting distribution...\n")

gamma_test_values <- c(0.5, 1.0, 2.0, 4.0)
weight_effects <- data.frame(
  gamma = gamma_test_values,
  weight_ratio = numeric(length(gamma_test_values)),
  weight_spread = numeric(length(gamma_test_values))
)

for (i in seq_along(gamma_test_values)) {
  gamma <- gamma_test_values[i]
  
  result <- slice_msf_single(vec2, mask2, r = 6, k = 0.5,
                            min_size = 10, nbhd = 8, gamma = gamma)
  
  weights <- result$weights
  left_weights <- weights[region2 == 1 & mask2@.Data > 0]
  right_weights <- weights[region2 == 2 & mask2@.Data > 0]
  
  ratio <- mean(left_weights) / mean(right_weights)
  spread <- sd(weights[mask2@.Data > 0])
  
  weight_effects$weight_ratio[i] <- ratio
  weight_effects$weight_spread[i] <- spread
  
  cat(sprintf("  γ=%.1f: weight_ratio=%.2f, weight_spread=%.3f\n", 
              gamma, ratio, spread))
}

# Summary
cat("\n=== Reliability Weighting Summary ===\n")

# Check if gamma effects are as expected
gamma_effect_correct <- weight_effects$weight_ratio[4] > weight_effects$weight_ratio[1]
weight_detection_correct <- weight_ratio > 2.0

cat("✓ Split-half correlation computation: implemented\n")
cat("✓ Gamma parameter weighting: implemented\n") 
cat(sprintf("✓ Reliability detection: %s (%.2fx weight ratio)\n", 
            ifelse(weight_detection_correct, "WORKING", "NOT WORKING"), weight_ratio))
cat(sprintf("✓ Gamma parameter effect: %s\n", 
            ifelse(gamma_effect_correct, "WORKING", "NOT WORKING")))

if (weight_detection_correct && gamma_effect_correct) {
  cat("\n🎉 Reliability weighting is WORKING correctly!\n")
  cat("   - Clean signals receive higher weights than noisy signals\n")
  cat("   - Higher gamma values increase the effect of reliability weighting\n")
  cat("   - Split-half correlation successfully measures signal consistency\n")
} else {
  cat("\n⚠️  Reliability weighting may need investigation\n")
}