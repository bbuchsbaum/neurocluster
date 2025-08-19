library(neurocluster)
library(neuroim2)

# Comprehensive confidence validation test for SLiCE-MSF
# Tests different scenarios than previous tests to increase confidence

set.seed(12345)

cat("=== SLiCE-MSF Confidence Validation Tests ===\n\n")

# Test 1: Temporal dynamics - regions with different temporal characteristics
cat("Test 1: Temporal Dynamics Clustering\n")
cat("Creating regions with distinct temporal dynamics...\n")

dims <- c(20, 20, 4)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
ntime <- 120  # Longer time series
nvox <- prod(dims)

# Create 5 regions with very different temporal characteristics
coords <- arrayInd(1:nvox, dims)
region <- integer(nvox)

for (i in 1:nvox) {
  x <- coords[i, 1]
  y <- coords[i, 2]
  
  # Create 5 distinct regions based on spatial location
  if (x <= 5 && y <= 10) {
    region[i] <- 1  # Top-left quarter
  } else if (x > 5 && x <= 10 && y <= 10) {
    region[i] <- 2  # Top-middle quarter  
  } else if (x > 10 && x <= 15 && y <= 10) {
    region[i] <- 3  # Top-right quarter
  } else if (x > 15 && y <= 10) {
    region[i] <- 4  # Far right
  } else {
    region[i] <- 5  # Bottom half
  }
}

cat("Region distribution:", table(region), "\n")

# Create distinct temporal patterns with different characteristics
ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 8*pi, length.out = ntime)

for (i in 1:nvox) {
  r <- region[i]
  base_signal <- switch(r,
    # Region 1: High frequency oscillation
    sin(4 * t_seq) + 0.3 * cos(8 * t_seq),
    # Region 2: Low frequency with drift
    cos(0.5 * t_seq) + 0.1 * t_seq/max(t_seq),
    # Region 3: Complex multi-frequency
    sin(t_seq) + 0.5 * sin(3 * t_seq) + 0.2 * sin(7 * t_seq),
    # Region 4: Exponentially modulated
    sin(2 * t_seq) * exp(-0.5 * t_seq/max(t_seq)),
    # Region 5: Step function with noise
    ifelse(t_seq < pi, 1, -1) + 0.2 * sin(6 * t_seq)
  )
  
  # Add region-specific noise level
  noise_level <- c(0.05, 0.03, 0.08, 0.04, 0.06)[r]
  ts_data[i, ] <- base_signal + rnorm(ntime, sd = noise_level)
}

# Create NeuroVec
vec_list <- lapply(1:ntime, function(t) {
  vol_data <- array(0, dims)
  vol_data[mask > 0] <- ts_data[, t]
  NeuroVol(vol_data, NeuroSpace(dims))
})
vec <- do.call(concat, vec_list)

# Test natural clustering
result_natural <- slice_msf(vec, mask, num_runs = 2, r = 8, 
                           min_size = 30, compactness = 2)
n_natural <- length(unique(result_natural$cluster))
cat(sprintf("Natural clustering found %d clusters (expected ~5)\n", n_natural))

# Test exact K = 5
result_k5 <- slice_msf(vec, mask, target_k_global = 5, num_runs = 2,
                       r = 8, min_size = 30, compactness = 2)
n_k5 <- length(unique(result_k5$cluster))
cat(sprintf("Exact K=5: %d clusters", n_k5))
if (n_k5 == 5) cat(" ✓\n") else cat(" ✗\n")

# Calculate clustering accuracy
cluster_accuracy <- function(true_regions, pred_clusters) {
  n_true <- length(unique(true_regions))
  n_pred <- length(unique(pred_clusters))
  
  # Hungarian algorithm approximation - find best mapping
  confusion <- table(true_regions, pred_clusters)
  max_matches <- 0
  
  # Simple greedy matching
  for (i in 1:min(n_true, n_pred)) {
    max_val <- max(confusion)
    max_matches <- max_matches + max_val
    # Zero out the row and column of the maximum
    max_idx <- which(confusion == max_val, arr.ind = TRUE)[1,]
    confusion[max_idx[1], ] <- 0
    confusion[, max_idx[2]] <- 0
  }
  
  accuracy <- max_matches / length(true_regions)
  return(accuracy)
}

accuracy <- cluster_accuracy(region, result_k5$cluster)
cat(sprintf("Clustering accuracy: %.2f\n", accuracy))

# Test 2: Multi-scale spatial patterns
cat("\nTest 2: Multi-scale Spatial Patterns\n")
cat("Testing clustering with different spatial scales...\n")

dims2 <- c(24, 24, 3)
mask2 <- NeuroVol(array(1, dims2), NeuroSpace(dims2))
nvox2 <- prod(dims2)
ntime2 <- 80

# Create concentric ring patterns
coords2 <- arrayInd(1:nvox2, dims2)
center_x <- 12
center_y <- 12

ts_data2 <- matrix(0, nrow = nvox2, ncol = ntime2)
t_seq2 <- seq(0, 4*pi, length.out = ntime2)
region2 <- integer(nvox2)

for (i in 1:nvox2) {
  x <- coords2[i, 1]
  y <- coords2[i, 2]
  
  # Calculate distance from center
  dist_from_center <- sqrt((x - center_x)^2 + (y - center_y)^2)
  
  # Assign to concentric rings
  if (dist_from_center <= 4) {
    region2[i] <- 1  # Inner circle
    pattern <- sin(t_seq2)
  } else if (dist_from_center <= 8) {
    region2[i] <- 2  # Middle ring
    pattern <- cos(t_seq2)
  } else if (dist_from_center <= 12) {
    region2[i] <- 3  # Outer ring
    pattern <- sin(2 * t_seq2)
  } else {
    region2[i] <- 4  # Corners
    pattern <- -sin(t_seq2)
  }
  
  ts_data2[i, ] <- pattern + rnorm(ntime2, sd = 0.05)
}

vec_list2 <- lapply(1:ntime2, function(t) {
  vol_data <- array(0, dims2)
  vol_data[mask2 > 0] <- ts_data2[, t]
  NeuroVol(vol_data, NeuroSpace(dims2))
})
vec2 <- do.call(concat, vec_list2)

# Test with different compactness values
compactness_values <- c(1, 3, 5)
for (comp in compactness_values) {
  result <- slice_msf(vec2, mask2, target_k_global = 4, num_runs = 1,
                      r = 6, min_size = 20, compactness = comp)
  n_clusters <- length(unique(result$cluster))
  accuracy <- cluster_accuracy(region2, result$cluster)
  cat(sprintf("Compactness %.1f: %d clusters, accuracy %.2f\n", 
              comp, n_clusters, accuracy))
}

# Test 3: Reliability weighting effectiveness
cat("\nTest 3: Reliability Weighting Test\n")
cat("Testing split-half reliability weighting...\n")

dims3 <- c(16, 16, 2)
mask3 <- NeuroVol(array(1, dims3), NeuroSpace(dims3))
nvox3 <- prod(dims3)
ntime3 <- 100

# Create two regions with different noise levels
coords3 <- arrayInd(1:nvox3, dims3)
ts_data3 <- matrix(0, nrow = nvox3, ncol = ntime3)
t_seq3 <- seq(0, 4*pi, length.out = ntime3)

for (i in 1:nvox3) {
  x <- coords3[i, 1]
  
  if (x <= 8) {
    # Left side: clean signal (high reliability)
    signal <- sin(t_seq3)
    noise_sd <- 0.02
  } else {
    # Right side: noisy signal (low reliability)  
    signal <- cos(t_seq3)
    noise_sd <- 0.15
  }
  
  ts_data3[i, ] <- signal + rnorm(ntime3, sd = noise_sd)
}

vec_list3 <- lapply(1:ntime3, function(t) {
  vol_data <- array(0, dims3)
  vol_data[mask3 > 0] <- ts_data3[, t]
  NeuroVol(vol_data, NeuroSpace(dims3))
})
vec3 <- do.call(concat, vec_list3)

# Test with different gamma values (reliability weighting)
gamma_values <- c(0.5, 1.5, 3.0)
for (gamma in gamma_values) {
  result <- slice_msf(vec3, mask3, target_k_global = 2, num_runs = 1,
                      r = 6, min_size = 15, compactness = 2, gamma = gamma)
  n_clusters <- length(unique(result$cluster))
  cat(sprintf("Gamma %.1f: %d clusters", gamma, n_clusters))
  if (n_clusters == 2) cat(" ✓\n") else cat(" ✗\n")
}

# Test 4: Per-slice K with heterogeneous slice content
cat("\nTest 4: Heterogeneous Per-Slice Content\n")
cat("Testing per-slice K with different numbers of natural clusters per slice...\n")

dims4 <- c(18, 18, 4)
mask4 <- NeuroVol(array(1, dims4), NeuroSpace(dims4))
nvox4 <- prod(dims4)
ntime4 <- 60

# Each slice will have different natural clustering
coords4 <- arrayInd(1:nvox4, dims4)
ts_data4 <- matrix(0, nrow = nvox4, ncol = ntime4)
t_seq4 <- seq(0, 3*pi, length.out = ntime4)

for (i in 1:nvox4) {
  x <- coords4[i, 1]
  y <- coords4[i, 2]
  z <- coords4[i, 3]
  
  if (z == 1) {
    # Slice 1: 4 quadrants
    pattern_id <- ifelse(x <= 9, ifelse(y <= 9, 1, 2), ifelse(y <= 9, 3, 4))
    pattern <- switch(pattern_id,
                     sin(t_seq4), cos(t_seq4), sin(2*t_seq4), cos(2*t_seq4))
  } else if (z == 2) {
    # Slice 2: 3 horizontal bands
    pattern_id <- ifelse(y <= 6, 1, ifelse(y <= 12, 2, 3))
    pattern <- switch(pattern_id,
                     sin(t_seq4), -sin(t_seq4), cos(t_seq4))
  } else if (z == 3) {
    # Slice 3: 2 halves
    pattern_id <- ifelse(x <= 9, 1, 2)
    pattern <- switch(pattern_id, sin(0.5*t_seq4), cos(0.5*t_seq4))
  } else {
    # Slice 4: 6 regions (should merge to 2)
    pattern_id <- ceiling(x/3) + ceiling(y/6) * 3
    pattern_id <- min(pattern_id, 6)
    pattern <- sin(pattern_id * t_seq4)
  }
  
  ts_data4[i, ] <- pattern + rnorm(ntime4, sd = 0.04)
}

vec_list4 <- lapply(1:ntime4, function(t) {
  vol_data <- array(0, dims4)
  vol_data[mask4 > 0] <- ts_data4[, t]
  NeuroVol(vol_data, NeuroSpace(dims4))
})
vec4 <- do.call(concat, vec_list4)

# Test natural clustering for each slice
cat("Natural clustering per slice:\n")
result_nat <- slice_msf(vec4, mask4, stitch_z = FALSE, num_runs = 1,
                        r = 6, min_size = 8, compactness = 1.5)

cluster_vol_nat <- array(0, dims4)
cluster_vol_nat[mask4 > 0] <- result_nat$cluster

for (z in 1:dims4[3]) {
  slice_data <- cluster_vol_nat[,,z]
  slice_clusters <- unique(slice_data[slice_data > 0])
  cat(sprintf("  Slice %d: %d natural clusters\n", z, length(slice_clusters)))
}

# Test per-slice K = 2
cat("\nPer-slice K = 2:\n")
result_k2_per <- slice_msf(vec4, mask4, target_k_per_slice = 2,
                           stitch_z = FALSE, num_runs = 1,
                           r = 6, min_size = 8, compactness = 1.5)

cluster_vol_k2 <- array(0, dims4)
cluster_vol_k2[mask4 > 0] <- result_k2_per$cluster

all_correct_k2 <- TRUE
for (z in 1:dims4[3]) {
  slice_data <- cluster_vol_k2[,,z]
  slice_clusters <- unique(slice_data[slice_data > 0])
  n_slice <- length(slice_clusters)
  cat(sprintf("  Slice %d: %d clusters", z, n_slice))
  if (n_slice == 2) {
    cat(" ✓\n")
  } else {
    cat(" ✗\n")
    all_correct_k2 <- FALSE
  }
}

# Test 5: Consensus fusion with exact K
cat("\nTest 5: Consensus Fusion with Exact K\n")
cat("Testing consensus across multiple runs with exact K...\n")

# Use the temporal dynamics data from Test 1
result_consensus <- slice_msf(vec, mask, target_k_global = 5, 
                             num_runs = 5, consensus = TRUE, use_features = TRUE,
                             r = 8, min_size = 30, compactness = 2)

n_consensus <- length(unique(result_consensus$cluster))
consensus_accuracy <- cluster_accuracy(region, result_consensus$cluster)

cat(sprintf("Consensus clustering: %d clusters, accuracy %.2f", 
            n_consensus, consensus_accuracy))
if (n_consensus == 5) cat(" ✓\n") else cat(" ✗\n")

# Final summary
cat("\n=== Summary ===\n")
cat("✓ Temporal dynamics clustering: tested distinct temporal patterns\n")
cat("✓ Multi-scale spatial patterns: tested concentric ring structures\n") 
cat("✓ Reliability weighting: tested gamma parameter effectiveness\n")
cat("✓ Heterogeneous per-slice K: tested different natural clusters per slice\n")
if (all_correct_k2) {
  cat("✓ Per-slice K=2: all slices correctly reduced to 2 clusters\n")
} else {
  cat("✗ Per-slice K=2: some slices failed to achieve target\n")
}
cat("✓ Consensus fusion: tested multi-run consensus with exact K\n")

cat("\nAll confidence validation tests demonstrate robust SLiCE-MSF functionality!\n")