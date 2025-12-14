#!/usr/bin/env Rscript
# Example: Demonstrating Cluster ID Ordering Issue and Solution
#
# This script shows how cluster IDs become spatially ordered and how to fix it

library(neurocluster)
library(neuroim2)

cat("Cluster ID Ordering Example\n")
cat("===========================\n\n")

# Check if test data exists
if (!file.exists("testdata/rscan01.nii.gz") || !file.exists("testdata/mask.nii")) {
  cat("Test data not found. This example demonstrates the concept:\n\n")

  cat("# Load your data\n")
  cat("vec <- read_vec('your_data.nii.gz')\n")
  cat("mask <- read_vol('your_mask.nii')\n\n")

  cat("# Cluster with slice_msf (produces z-ordered IDs)\n")
  cat("result_ordered <- cluster4d(vec, mask, K = 100, method = 'slice_msf')\n\n")

  cat("# Check for spatial ordering\n")
  cat("ordering <- check_cluster_ordering(result_ordered)\n")
  cat("print(ordering)\n\n")

  cat("# If strongly ordered (|r| > 0.3), randomize\n")
  cat("if (ordering$is_ordered) {\n")
  cat("  result_random <- randomize_cluster_ids(result_ordered, seed = 42)\n")
  cat("  cat('✓ Cluster IDs randomized\\n')\n")
  cat("}\n\n")

  cat("# Compare correlations\n")
  cat("cat('Before randomization:',\n")
  cat("    cor(1:nrow(result_ordered$coord_centers),\n")
  cat("        result_ordered$coord_centers[,3],\n")
  cat("        method='spearman'), '\\n')\n\n")

  cat("cat('After randomization:',\n")
  cat("    cor(1:nrow(result_random$coord_centers),\n")
  cat("        result_random$coord_centers[,3],\n")
  cat("        method='spearman'), '\\n')\n")

} else {

  # Load test data
  cat("Loading test data...\n")
  vec <- read_vec("testdata/rscan01.nii.gz")
  mask <- read_vol("testdata/mask.nii")
  cat("Data loaded: ", dim(vec), "\n")
  cat("Mask voxels:", sum(mask > 0), "\n\n")

  # Test with slice_msf
  cat("Running slice_msf clustering (K=50)...\n")
  result_msf <- cluster4d(vec, mask, K = 50, method = "slice_msf",
                         num_runs = 1, verbose = FALSE)
  cat("Clustering complete.\n\n")

  # Check ordering
  cat("Checking for spatial ordering in cluster IDs:\n")
  cat(rep("-", 60), "\n", sep = "")
  ordering <- check_cluster_ordering(result_msf)

  cat("Correlation between cluster ID and coordinates:\n")
  cat("  X-axis: ", sprintf("%.3f", ordering$cor_x), "\n")
  cat("  Y-axis: ", sprintf("%.3f", ordering$cor_y), "\n")
  cat("  Z-axis: ", sprintf("%.3f", ordering$cor_z), " ← **",
      if (abs(ordering$cor_z) > 0.7) "STRONG" else if (abs(ordering$cor_z) > 0.3) "MODERATE" else "WEAK",
      " ordering**\n")
  cat("\n", ordering$interpretation, "\n\n")

  # Randomize
  if (ordering$is_ordered) {
    cat("Applying randomization...\n")
    result_random <- randomize_cluster_ids(result_msf, seed = 42)

    # Check again
    ordering_after <- check_cluster_ordering(result_random)

    cat("\nAfter randomization:\n")
    cat(rep("-", 60), "\n", sep = "")
    cat("Correlation between cluster ID and coordinates:\n")
    cat("  X-axis: ", sprintf("%.3f", ordering_after$cor_x), "\n")
    cat("  Y-axis: ", sprintf("%.3f", ordering_after$cor_y), "\n")
    cat("  Z-axis: ", sprintf("%.3f", ordering_after$cor_z), " ← **",
        if (abs(ordering_after$cor_z) > 0.3) "STILL ORDERED" else "FIXED",
        "**\n")
    cat("\n", ordering_after$interpretation, "\n\n")

    cat("✓ Success: Cluster IDs no longer spatially ordered\n\n")
  } else {
    cat("No strong ordering detected - randomization not needed.\n\n")
  }

  # Visualization comparison
  cat("\nVisualization Impact:\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("\nWith ORDERED cluster IDs:\n")
  cat("  - Continuous color scale shows artificial gradient\n")
  cat("  - Inferior regions (low z) → cool colors (blue/purple)\n")
  cat("  - Superior regions (high z) → warm colors (yellow/red)\n")
  cat("  - This gradient is MISLEADING - it's a computational artifact!\n\n")

  cat("With RANDOMIZED cluster IDs:\n")
  cat("  - Each cluster appears as a distinct color patch\n")
  cat("  - No artificial gradients\n")
  cat("  - Spatial structure reflects actual clustering, not ID ordering\n\n")

  cat("Recommendation:\n")
  cat("  ALWAYS use randomize_cluster_ids() before visualization\n")
  cat("  OR use discrete color palettes (e.g., RColorBrewer::Set3)\n\n")

  # Z-distribution analysis
  cat("\nZ-Distribution of Cluster Centroids:\n")
  cat(rep("=", 60), "\n", sep = "")

  z_coords <- result_msf$coord_centers[, 3]
  z_range <- range(z_coords)

  cat("\nBefore randomization (cluster IDs 1-", nrow(result_msf$coord_centers), "):\n", sep = "")
  for (q in 1:5) {
    ids_in_quintile <- which((1:nrow(result_msf$coord_centers)) <= q * nrow(result_msf$coord_centers) / 5 &
                             (1:nrow(result_msf$coord_centers)) > (q - 1) * nrow(result_msf$coord_centers) / 5)
    z_mean <- mean(z_coords[ids_in_quintile])
    cat(sprintf("  Cluster IDs %2d-%-2d: mean z = %.1f\n",
                min(ids_in_quintile), max(ids_in_quintile), z_mean))
  }
  cat("\n  ↑ Notice how mean z increases with cluster ID!\n")
  cat("  This creates artificial inferior → superior gradient\n\n")

  cat("After randomization, cluster IDs are decoupled from z-position.\n")
}

cat("\n", rep("=", 60), "\n", sep = "")
cat("Summary\n")
cat(rep("=", 60), "\n", sep = "")
cat("\n1. Many clustering algorithms produce spatially-ordered cluster IDs\n")
cat("2. This creates misleading gradients with continuous color scales\n")
cat("3. Use check_cluster_ordering() to detect the issue\n")
cat("4. Use randomize_cluster_ids() to fix it\n")
cat("5. Always randomize before visualization or publication\n\n")

cat("For more information, see:\n")
cat("  - CLUSTER_ID_ORDERING.md (comprehensive guide)\n")
cat("  - FLASH3D_FIX_SUMMARY.md (algorithm-level fix example)\n\n")
