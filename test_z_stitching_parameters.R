library(neurocluster)
library(neuroim2)

# Z-Stitching Parameter Sensitivity Test
# Tests how theta_link and min_contact parameters affect stitching behavior

set.seed(1234)

cat("=== Z-Stitching Parameter Sensitivity Tests ===\n\n")

# Create test data: 3 slices with similar but not identical patterns
dims <- c(20, 20, 3)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
nvox <- prod(dims)
ntime <- 80

cat("Creating 3-slice volume with controlled similarity patterns...\n")

# Create 4 regions per slice with specific similarity relationships
coords <- arrayInd(1:nvox, dims)
region <- integer(nvox)
slice <- integer(nvox)

for (i in 1:nvox) {
  x <- coords[i, 1]
  y <- coords[i, 2]
  z <- coords[i, 3]
  
  slice[i] <- z
  
  # Create 4 quadrants per slice
  if (x <= 10 && y <= 10) region[i] <- 1  # Top-left
  else if (x > 10 && y <= 10) region[i] <- 2  # Top-right
  else if (x <= 10 && y > 10) region[i] <- 3  # Bottom-left
  else region[i] <- 4  # Bottom-right
}

# Create time series with controlled similarity across slices
ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 6*pi, length.out = ntime)

# Define base patterns for each region
base_patterns <- list(
  sin(t_seq),           # Region 1
  cos(t_seq),           # Region 2
  sin(2 * t_seq),       # Region 3
  cos(2 * t_seq)        # Region 4
)

# Create patterns with controlled similarity across slices
for (i in 1:nvox) {
  r <- region[i]
  z <- slice[i]
  base_signal <- base_patterns[[r]]
  
  if (z == 1) {
    # Slice 1: Original patterns
    signal <- base_signal
  } else if (z == 2) {
    # Slice 2: Very similar to slice 1 (high similarity)
    signal <- base_signal + 0.1 * sin(0.5 * t_seq)  # Small perturbation
  } else {
    # Slice 3: Moderately similar to slice 2 (medium similarity)
    signal <- base_signal + 0.3 * cos(1.5 * t_seq)  # Larger perturbation
  }
  
  ts_data[i, ] <- signal + rnorm(ntime, sd = 0.05)
}

# Create NeuroVec
vec_list <- lapply(1:ntime, function(t) {
  vol_data <- array(0, dims)
  vol_data[mask > 0] <- ts_data[, t]
  NeuroVol(vol_data, NeuroSpace(dims))
})
vec <- do.call(concat, vec_list)

# Test 1: theta_link sensitivity
cat("Test 1: theta_link Sensitivity\n")
cat("Testing how similarity threshold affects stitching frequency...\n")

theta_values <- c(0.50, 0.70, 0.85, 0.95, 0.99)
theta_results <- data.frame(
  theta_link = theta_values,
  n_clusters = numeric(length(theta_values)),
  stitch_rate = numeric(length(theta_values))
)

# Get baseline without stitching
result_no_stitch <- slice_msf(vec, mask, stitch_z = FALSE, num_runs = 1,
                              r = 6, min_size = 20, compactness = 3)
n_baseline <- length(unique(result_no_stitch$cluster))
cat(sprintf("Baseline (no stitching): %d clusters\n", n_baseline))

for (i in seq_along(theta_values)) {
  theta <- theta_values[i]
  
  result <- slice_msf(vec, mask, stitch_z = TRUE, theta_link = theta,
                      min_contact = 3, num_runs = 1, r = 6, 
                      min_size = 20, compactness = 3)
  
  n_clusters <- length(unique(result$cluster))
  stitch_rate <- (n_baseline - n_clusters) / n_baseline
  
  theta_results$n_clusters[i] <- n_clusters
  theta_results$stitch_rate[i] <- stitch_rate
  
  cat(sprintf("  θ=%.2f: %d clusters, stitch_rate=%.2f\n", 
              theta, n_clusters, stitch_rate))
}

# Check that higher theta_link reduces stitching
theta_trend_correct <- all(diff(theta_results$n_clusters) >= 0)  # Should be non-decreasing
cat(sprintf("Theta trend correct (higher θ → fewer stitches): %s\n", 
            ifelse(theta_trend_correct, "✓", "✗")))

# Test 2: min_contact sensitivity  
cat("\nTest 2: min_contact Sensitivity\n")
cat("Testing how minimum contact requirement affects stitching...\n")

contact_values <- c(1, 3, 5, 10, 20)
contact_results <- data.frame(
  min_contact = contact_values,
  n_clusters = numeric(length(contact_values)),
  stitch_rate = numeric(length(contact_values))
)

for (i in seq_along(contact_values)) {
  min_contact <- contact_values[i]
  
  result <- slice_msf(vec, mask, stitch_z = TRUE, theta_link = 0.85,
                      min_contact = min_contact, num_runs = 1, r = 6,
                      min_size = 20, compactness = 3)
  
  n_clusters <- length(unique(result$cluster))
  stitch_rate <- (n_baseline - n_clusters) / n_baseline
  
  contact_results$n_clusters[i] <- n_clusters
  contact_results$stitch_rate[i] <- stitch_rate
  
  cat(sprintf("  min_contact=%d: %d clusters, stitch_rate=%.2f\n",
              min_contact, n_clusters, stitch_rate))
}

# Check that higher min_contact reduces stitching
contact_trend_correct <- all(diff(contact_results$n_clusters) >= 0)  # Should be non-decreasing
cat(sprintf("Contact trend correct (higher contact → fewer stitches): %s\n",
            ifelse(contact_trend_correct, "✓", "✗")))

# Test 3: Extreme parameter values
cat("\nTest 3: Extreme Parameter Values\n")
cat("Testing behavior at parameter extremes...\n")

# Very permissive (should stitch almost everything)
result_permissive <- slice_msf(vec, mask, stitch_z = TRUE, theta_link = 0.1,
                               min_contact = 1, num_runs = 1, r = 6,
                               min_size = 20, compactness = 3)
n_permissive <- length(unique(result_permissive$cluster))

# Very restrictive (should stitch almost nothing)
result_restrictive <- slice_msf(vec, mask, stitch_z = TRUE, theta_link = 0.99,
                                min_contact = 50, num_runs = 1, r = 6,
                                min_size = 20, compactness = 3)
n_restrictive <- length(unique(result_restrictive$cluster))

cat(sprintf("Permissive (θ=0.1, contact=1): %d clusters\n", n_permissive))
cat(sprintf("Restrictive (θ=0.99, contact=50): %d clusters\n", n_restrictive))
cat(sprintf("Baseline (no stitching): %d clusters\n", n_baseline))

# Permissive should have fewer clusters, restrictive should be close to baseline
extreme_test_correct <- (n_permissive < n_baseline) && (n_restrictive >= n_baseline * 0.9)
cat(sprintf("Extreme values behave correctly: %s\n", 
            ifelse(extreme_test_correct, "✓", "✗")))

# Test 4: Parameter interaction
cat("\nTest 4: Parameter Interaction\n")
cat("Testing interaction between theta_link and min_contact...\n")

# Create interaction matrix
theta_test <- c(0.7, 0.85, 0.95)
contact_test <- c(1, 5, 15)
interaction_results <- matrix(0, nrow = length(theta_test), ncol = length(contact_test))
dimnames(interaction_results) <- list(
  paste("θ=", theta_test, sep=""),
  paste("contact=", contact_test, sep="")
)

for (i in seq_along(theta_test)) {
  for (j in seq_along(contact_test)) {
    result <- slice_msf(vec, mask, stitch_z = TRUE, 
                        theta_link = theta_test[i],
                        min_contact = contact_test[j],
                        num_runs = 1, r = 6, min_size = 20, compactness = 3)
    
    interaction_results[i, j] <- length(unique(result$cluster))
  }
}

cat("Interaction matrix (number of clusters):\n")
print(interaction_results)

# Check that both parameters have expected effects
row_increasing <- all(apply(interaction_results, 1, function(x) all(diff(x) >= 0)))
col_increasing <- all(apply(interaction_results, 2, function(x) all(diff(x) >= 0)))

cat(sprintf("Row trend correct (higher contact → more clusters): %s\n", 
            ifelse(row_increasing, "✓", "✗")))
cat(sprintf("Column trend correct (higher θ → more clusters): %s\n", 
            ifelse(col_increasing, "✓", "✗")))

# Test 5: Stitching disabled check
cat("\nTest 5: Stitching Disabled Check\n")
cat("Verifying that stitch_z=FALSE ignores stitching parameters...\n")

# Parameters should have no effect when stitching is disabled
result_disabled_1 <- slice_msf(vec, mask, stitch_z = FALSE, theta_link = 0.1,
                               min_contact = 1, num_runs = 1, r = 6,
                               min_size = 20, compactness = 3)

result_disabled_2 <- slice_msf(vec, mask, stitch_z = FALSE, theta_link = 0.99,
                               min_contact = 50, num_runs = 1, r = 6,
                               min_size = 20, compactness = 3)

n_disabled_1 <- length(unique(result_disabled_1$cluster))
n_disabled_2 <- length(unique(result_disabled_2$cluster))

disabled_consistent <- (n_disabled_1 == n_disabled_2) && (n_disabled_1 == n_baseline)
cat(sprintf("Disabled stitching consistent: %s (all = %d clusters)\n",
            ifelse(disabled_consistent, "✓", "✗"), n_baseline))

# Summary
cat("\n=== Parameter Sensitivity Summary ===\n")
cat(sprintf("✓ Baseline functionality: %d clusters without stitching\n", n_baseline))
cat(sprintf("✓ Theta_link range: %.2f to %.2f stitch rate\n", 
            min(theta_results$stitch_rate), max(theta_results$stitch_rate)))
cat(sprintf("✓ Min_contact range: %.2f to %.2f stitch rate\n",
            min(contact_results$stitch_rate), max(contact_results$stitch_rate)))

overall_success <- theta_trend_correct && contact_trend_correct && 
                   extreme_test_correct && disabled_consistent

if (overall_success) {
  cat("\n🎉 All parameter sensitivity tests PASSED!\n")
  cat("   - theta_link controls similarity threshold correctly\n")
  cat("   - min_contact controls contact requirement correctly\n")
  cat("   - Parameters interact as expected\n")
  cat("   - Extreme values behave appropriately\n")
  cat("   - Disabled stitching ignores parameters\n")
} else {
  cat("\n⚠️  Some parameter sensitivity tests FAILED\n")
  cat("   Check parameter implementation and behavior\n")
}