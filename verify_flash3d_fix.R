#!/usr/bin/env Rscript
# Verification script for FLASH-3D z-axis striping fix
# This script tests that the spatial-balanced seed selection eliminates striping artifacts

library(neurocluster)
library(neuroim2)

cat("FLASH-3D Fix Verification Script\n")
cat("=================================\n\n")

# Load test data
cat("Loading test data...\n")
if (file.exists("testdata/rscan01.nii.gz") && file.exists("testdata/mask.nii")) {
  vec <- read_vec("testdata/rscan01.nii.gz")
  mask <- read_vol("testdata/mask.nii")

  cat("Data dimensions:", dim(vec), "\n")
  cat("Mask voxels:", sum(mask > 0), "\n\n")

  # Test parameters that previously showed striping
  test_configs <- list(
    list(K = 100, bits = 64, dctM = 16, lambda_s = 0.6, lambda_t = 1.0, label = "Config 1: K=100, bits=64, dctM=16"),
    list(K = 150, bits = 128, dctM = 32, lambda_s = 0.1, lambda_t = 1.0, label = "Config 2: K=150, bits=128, dctM=32"),
    list(K = 200, bits = 64, dctM = 24, lambda_s = 10.0, lambda_t = 1.0, label = "Config 3: K=200, bits=64, dctM=24, high spatial weight")
  )

  for (config in test_configs) {
    cat("\n", rep("=", 60), "\n", sep = "")
    cat("Testing:", config$label, "\n")
    cat(rep("=", 60), "\n", sep = "")

    result <- cluster4d(vec, mask, K = config$K, method = "flash3d",
                       bits = config$bits, dctM = config$dctM,
                       lambda_s = config$lambda_s, lambda_t = config$lambda_t,
                       verbose = TRUE)

    # Analyze z-coordinate distribution
    coords <- result$coord_centers
    z_coords <- coords[, 3]

    cat("\nCluster Spatial Distribution Analysis:\n")
    cat("--------------------------------------\n")
    cat("Total clusters:", nrow(coords), "\n")
    cat("Z-coordinate range: [", min(z_coords), ",", max(z_coords), "]\n", sep = "")
    cat("Z-coordinate mean:", round(mean(z_coords), 2), "\n")
    cat("Z-coordinate SD:", round(sd(z_coords), 2), "\n")

    # Check for z-ordering in cluster IDs
    z_ordered <- z_coords[order(seq_along(z_coords))]
    spearman_cor <- cor(seq_along(z_coords), z_coords, method = "spearman")

    cat("\nZ-Ordering Test:\n")
    cat("----------------\n")
    cat("Spearman correlation (cluster ID vs z-position):", round(spearman_cor, 3), "\n")
    if (abs(spearman_cor) > 0.7) {
      cat("⚠️  WARNING: Strong correlation detected (|r| > 0.7)\n")
      cat("   This suggests striping may still be present.\n")
    } else if (abs(spearman_cor) > 0.3) {
      cat("⚠️  CAUTION: Moderate correlation detected (|r| > 0.3)\n")
      cat("   Some ordering bias may remain.\n")
    } else {
      cat("✓ PASS: Low correlation (|r| < 0.3)\n")
      cat("   Cluster IDs appear spatially unordered.\n")
    }

    # Z-distribution histogram
    cat("\nZ-Distribution by Quintile:\n")
    z_range <- max(z_coords) - min(z_coords)
    for (q in 1:5) {
      lower <- min(z_coords) + (q - 1) * z_range / 5
      upper <- min(z_coords) + q * z_range / 5
      n_in_quintile <- sum(z_coords >= lower & z_coords < upper)
      cat(sprintf("  Q%d [%.1f-%.1f]: %d clusters (%.1f%%)\n",
                  q, lower, upper, n_in_quintile,
                  100 * n_in_quintile / length(z_coords)))
    }

    # Spatial uniformity test (chi-square)
    expected <- length(z_coords) / 5
    observed <- sapply(1:5, function(q) {
      lower <- min(z_coords) + (q - 1) * z_range / 5
      upper <- min(z_coords) + q * z_range / 5
      sum(z_coords >= lower & z_coords < upper)
    })
    chi_sq <- sum((observed - expected)^2 / expected)
    p_value <- pchisq(chi_sq, df = 4, lower.tail = FALSE)

    cat("\nUniformity Test (Chi-square):\n")
    cat("-----------------------------\n")
    cat("Chi-square statistic:", round(chi_sq, 2), "\n")
    cat("P-value:", round(p_value, 4), "\n")
    if (p_value > 0.05) {
      cat("✓ PASS: Z-distribution is spatially uniform (p > 0.05)\n")
    } else {
      cat("⚠️  FAIL: Z-distribution is non-uniform (p < 0.05)\n")
      cat("   Striping artifacts may be present.\n")
    }
  }

  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Verification Complete\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("\nSummary:\n")
  cat("--------\n")
  cat("If all tests show:\n")
  cat("  - Low Spearman correlation (|r| < 0.3)\n")
  cat("  - Uniform z-distribution (p > 0.05)\n")
  cat("Then the fix has successfully eliminated striping artifacts.\n\n")

} else {
  cat("ERROR: Test data not found.\n")
  cat("Please ensure testdata/rscan01.nii.gz and testdata/mask.nii exist.\n")
  cat("\nAlternatively, test with your own data:\n")
  cat("  result <- cluster4d(your_vec, your_mask, K = 100, method = 'flash3d')\n")
  cat("  coords <- result$coord_centers\n")
  cat("  cor(1:nrow(coords), coords[, 3], method = 'spearman')  # Should be near 0\n")
}
