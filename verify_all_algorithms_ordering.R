#!/usr/bin/env Rscript
# Systematic verification of z-ordering across all clustering algorithms
#
# This script tests whether cluster IDs are systematically ordered along the z-axis
# for all clustering methods in neurocluster package.

library(neurocluster)
library(neuroim2)

cat("Cluster ID Z-Ordering Verification\n")
cat("===================================\n\n")

# Check if test data exists
if (!file.exists("testdata/rscan01.nii.gz") || !file.exists("testdata/mask.nii")) {
  cat("Test data not found. Using synthetic data instead.\n\n")

  # Create synthetic 3D data with clear z-structure
  mask <- NeuroVol(array(1, c(20, 20, 20)), NeuroSpace(c(20, 20, 20)))
  vec <- replicate(50,
                   NeuroVol(array(rnorm(20*20*20), c(20, 20, 20)),
                            NeuroSpace(c(20, 20, 20))),
                   simplify = FALSE)
  vec <- do.call(concat, vec)

} else {
  # Load test data
  cat("Loading test data...\n")
  vec <- read_vec("testdata/rscan01.nii.gz")
  mask <- read_vol("testdata/mask.nii")
}

cat("Data dimensions:", dim(vec), "\n")
cat("Mask voxels:", sum(mask > 0), "\n\n")

# Function to check z-ordering for a clustering result
check_z_ordering <- function(result, method_name) {
  cat("\n", rep("=", 60), "\n", sep = "")
  cat(method_name, "\n")
  cat(rep("=", 60), "\n", sep = "")

  if (is.null(result$coord_centers) || nrow(result$coord_centers) < 3) {
    cat("  [SKIPPED] Too few clusters for correlation analysis\n")
    return(NULL)
  }

  # Check for spatial ordering
  ordering <- check_cluster_ordering(result)

  cat("Correlation between cluster ID and coordinates:\n")
  cat("  X-axis:", sprintf("%.3f", ordering$cor_x), "\n")
  cat("  Y-axis:", sprintf("%.3f", ordering$cor_y), "\n")
  cat("  Z-axis:", sprintf("%.3f", ordering$cor_z), " ← ",
      if (abs(ordering$cor_z) > 0.7) "[STRONG ORDERING]"
      else if (abs(ordering$cor_z) > 0.3) "[MODERATE ORDERING]"
      else "[WEAK/NO ORDERING]", "\n")

  cat("\n", ordering$interpretation, "\n")

  # Quintile analysis
  z_coords <- result$coord_centers[, 3]
  n_clusters <- nrow(result$coord_centers)
  cluster_ids <- seq_len(n_clusters)

  cat("\nZ-Distribution by Cluster ID Quintile:\n")
  for (q in 1:5) {
    ids_in_quintile <- which(cluster_ids <= q * n_clusters / 5 &
                               cluster_ids > (q - 1) * n_clusters / 5)
    z_mean <- mean(z_coords[ids_in_quintile])
    z_sd <- sd(z_coords[ids_in_quintile])
    cat(sprintf("  Quintile %d (IDs %3d-%-3d): mean z = %6.2f ± %.2f\n",
                q, min(ids_in_quintile), max(ids_in_quintile), z_mean, z_sd))
  }

  if (ordering$is_ordered) {
    cat("\n⚠ ORDERING DETECTED - Randomization recommended before visualization\n")
  } else {
    cat("\n✓ No significant ordering detected\n")
  }

  return(ordering)
}

# Test all algorithms
results <- list()

cat("\n\nTesting clustering algorithms with K=50...\n")
cat(rep("=", 60), "\n", sep = "")

# 1. Supervoxels
cat("\n[1/6] Testing supervoxels...\n")
tryCatch({
  result_supervoxels <- supervoxels(vec, mask, K = 50, sigma1 = 1, sigma2 = 2.5,
                                    iterations = 10, verbose = FALSE)
  results$supervoxels <- check_z_ordering(result_supervoxels, "SUPERVOXELS")
}, error = function(e) {
  cat("  [ERROR]", e$message, "\n")
})

# 2. SNIC
cat("\n[2/6] Testing snic...\n")
tryCatch({
  result_snic <- snic(vec, mask, K = 50, compactness = 5)
  results$snic <- check_z_ordering(result_snic, "SNIC")
}, error = function(e) {
  cat("  [ERROR]", e$message, "\n")
})

# 3. slice_msf
cat("\n[3/6] Testing slice_msf...\n")
tryCatch({
  result_slice_msf <- slice_msf(vec, mask, target_k_global = 50,
                                num_runs = 1, consensus = FALSE)
  results$slice_msf <- check_z_ordering(result_slice_msf, "SLICE_MSF")
}, error = function(e) {
  cat("  [ERROR]", e$message, "\n")
})

# 4. FLASH-3D
cat("\n[4/6] Testing flash3d...\n")
tryCatch({
  result_flash3d <- supervoxels_flash3d(vec, mask, K = 50, bits = 64,
                                        dctM = 16, verbose = FALSE)
  results$flash3d <- check_z_ordering(result_flash3d, "FLASH3D (FIXED)")
}, error = function(e) {
  cat("  [ERROR]", e$message, "\n")
})

# 5. SLIC4D
cat("\n[5/6] Testing slic4d...\n")
tryCatch({
  result_slic4d <- slic4d(vec, mask, K = 50, compactness = 5)
  results$slic4d <- check_z_ordering(result_slic4d, "SLIC4D")
}, error = function(e) {
  cat("  [ERROR]", e$message, "\n")
})

# 6. ACSC
cat("\n[6/6] Testing acsc...\n")
tryCatch({
  result_acsc <- acsc(vec, mask, K = 50)
  results$acsc <- check_z_ordering(result_acsc, "ACSC")
}, error = function(e) {
  cat("  [ERROR]", e$message, "\n")
})

# Summary
cat("\n\n", rep("=", 60), "\n", sep = "")
cat("SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")

cat("\nAlgorithms with Z-Ordering (|r_z| > 0.3):\n")
for (method in names(results)) {
  if (!is.null(results[[method]]) && results[[method]]$is_ordered) {
    cat(sprintf("  ✗ %-15s: r_z = %6.3f (%s ordering)\n",
                toupper(method),
                results[[method]]$cor_z,
                if (abs(results[[method]]$cor_z) > 0.7) "STRONG" else "MODERATE"))
  }
}

cat("\nAlgorithms without Z-Ordering (|r_z| ≤ 0.3):\n")
for (method in names(results)) {
  if (!is.null(results[[method]]) && !results[[method]]$is_ordered) {
    cat(sprintf("  ✓ %-15s: r_z = %6.3f\n",
                toupper(method),
                results[[method]]$cor_z))
  }
}

cat("\n\nRecommendations:\n")
cat("  1. For algorithms with z-ordering: Use randomize_cluster_ids() before visualization\n")
cat("  2. Always use discrete color scales OR randomize IDs when using continuous scales\n")
cat("  3. Report ordering status in publication methods sections\n")
cat("\nSee CLUSTER_ID_ORDERING.md for detailed guidance.\n\n")
