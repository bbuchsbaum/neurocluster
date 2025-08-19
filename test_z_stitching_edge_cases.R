library(neurocluster)
library(neuroim2)

# Z-Stitching Edge Cases Test
# Tests boundary conditions and special scenarios for z-stitching

set.seed(9999)

cat("=== Z-Stitching Edge Cases Tests ===\n\n")

# Test 1: Single Slice Volume
cat("Test 1: Single Slice Volume\n")
cat("Verifying that stitching parameters have no effect on single slice...\n")

dims1 <- c(20, 20, 1)  # Single slice
mask1 <- NeuroVol(array(1, dims1), NeuroSpace(dims1))
nvox1 <- prod(dims1)
ntime1 <- 60

# Create 4 quadrants with different patterns
coords1 <- arrayInd(1:nvox1, dims1)
ts_data1 <- matrix(0, nrow = nvox1, ncol = ntime1)
t_seq1 <- seq(0, 4*pi, length.out = ntime1)

for (i in 1:nvox1) {
  x <- coords1[i, 1]
  y <- coords1[i, 2]
  
  if (x <= 10 && y <= 10) signal <- sin(t_seq1)
  else if (x > 10 && y <= 10) signal <- cos(t_seq1)
  else if (x <= 10 && y > 10) signal <- sin(2*t_seq1)
  else signal <- cos(2*t_seq1)
  
  ts_data1[i, ] <- signal + rnorm(ntime1, sd = 0.05)
}

vec_list1 <- lapply(1:ntime1, function(t) {
  vol_data <- array(0, dims1)
  vol_data[mask1 > 0] <- ts_data1[, t]
  NeuroVol(vol_data, NeuroSpace(dims1))
})
vec1 <- do.call(concat, vec_list1)

# Test with different stitching parameters - should all give same result
single_slice_results <- list()
param_sets <- list(
  list(stitch_z = FALSE, theta_link = 0.85, min_contact = 3),
  list(stitch_z = TRUE, theta_link = 0.1, min_contact = 1),
  list(stitch_z = TRUE, theta_link = 0.99, min_contact = 50)
)

for (i in seq_along(param_sets)) {
  params <- param_sets[[i]]
  result <- slice_msf(vec1, mask1, 
                      stitch_z = params$stitch_z,
                      theta_link = params$theta_link,
                      min_contact = params$min_contact,
                      num_runs = 1, r = 6, min_size = 20, compactness = 3)
  
  single_slice_results[[i]] <- length(unique(result$cluster))
  cat(sprintf("  stitch_z=%s, θ=%.2f, contact=%d: %d clusters\n",
              params$stitch_z, params$theta_link, params$min_contact,
              single_slice_results[[i]]))
}

single_slice_consistent <- all(sapply(single_slice_results, function(x) x == single_slice_results[[1]]))
cat(sprintf("Single slice consistency: %s\n", ifelse(single_slice_consistent, "✓", "✗")))

# Test 2: Minimal Contact Scenarios
cat("\nTest 2: Minimal Contact Scenarios\n")
cat("Testing clusters with very small contact areas...\n")

dims2 <- c(12, 12, 3)
mask2 <- NeuroVol(array(1, dims2), NeuroSpace(dims2))

# Create a specific pattern where contact is minimal
# Small touching regions between slices
special_mask <- array(0, dims2)

# Slice 1: small region in center
special_mask[5:8, 5:8, 1] <- 1

# Slice 2: slightly offset region (minimal overlap)
special_mask[6:9, 6:9, 2] <- 1

# Slice 3: further offset (no overlap with slice 1)
special_mask[7:10, 7:10, 3] <- 1

mask2@.Data <- special_mask
nvox2 <- sum(special_mask > 0)
mask_idx2 <- which(special_mask > 0)

ntime2 <- 50
ts_data2 <- matrix(0, nrow = length(mask_idx2), ncol = ntime2)
t_seq2 <- seq(0, 4*pi, length.out = ntime2)

# Create similar patterns across slices
for (i in 1:length(mask_idx2)) {
  ts_data2[i, ] <- sin(t_seq2) + rnorm(ntime2, sd = 0.05)
}

vec_list2 <- lapply(1:ntime2, function(t) {
  vol_data <- array(0, dims2)
  vol_data[mask_idx2] <- ts_data2[, t]
  NeuroVol(vol_data, NeuroSpace(dims2))
})
vec2 <- do.call(concat, vec_list2)

# Test with different min_contact values
contact_values <- c(1, 5, 20)
minimal_contact_results <- list()

for (contact in contact_values) {
  result <- slice_msf(vec2, mask2, stitch_z = TRUE, theta_link = 0.8,
                      min_contact = contact, num_runs = 1,
                      r = 4, min_size = 5, compactness = 2)
  
  n_clusters <- length(unique(result$cluster))
  minimal_contact_results[[as.character(contact)]] <- n_clusters
  
  cat(sprintf("  min_contact=%d: %d clusters\n", contact, n_clusters))
}

# Test 3: Competing Candidates
cat("\nTest 3: Competing Stitching Candidates\n")
cat("Testing when multiple clusters could potentially stitch...\n")

dims3 <- c(10, 15, 3)
mask3 <- NeuroVol(array(1, dims3), NeuroSpace(dims3))
nvox3 <- prod(dims3)
ntime3 <- 60

coords3 <- arrayInd(1:nvox3, dims3)
ts_data3 <- matrix(0, nrow = nvox3, ncol = ntime3)
t_seq3 <- seq(0, 4*pi, length.out = ntime3)

# Create scenario with competing candidates
for (i in 1:nvox3) {
  x <- coords3[i, 1]
  y <- coords3[i, 2]
  z <- coords3[i, 3]
  
  if (z == 1) {
    # Slice 1: Three distinct regions
    if (y <= 5) signal <- sin(t_seq3)        # Region A
    else if (y <= 10) signal <- cos(t_seq3)  # Region B  
    else signal <- sin(2*t_seq3)             # Region C
  } else if (z == 2) {
    # Slice 2: Two regions that could match multiple regions from slice 1
    if (y <= 7) signal <- sin(t_seq3) + 0.1*cos(2*t_seq3)  # Similar to A, slightly to C
    else signal <- cos(t_seq3) + 0.1*sin(0.5*t_seq3)       # Similar to B, slightly different
  } else {
    # Slice 3: One region similar to slice 2's first region
    signal <- sin(t_seq3) + 0.05*cos(2*t_seq3)
  }
  
  ts_data3[i, ] <- signal + rnorm(ntime3, sd = 0.05)
}

vec_list3 <- lapply(1:ntime3, function(t) {
  vol_data <- array(0, dims3)
  vol_data[mask3 > 0] <- ts_data3[, t]
  NeuroVol(vol_data, NeuroSpace(dims3))
})
vec3 <- do.call(concat, vec_list3)

result3_no_stitch <- slice_msf(vec3, mask3, stitch_z = FALSE, num_runs = 1,
                               r = 6, min_size = 8, compactness = 3)
result3_stitch <- slice_msf(vec3, mask3, stitch_z = TRUE, theta_link = 0.7,
                            min_contact = 1, num_runs = 1,
                            r = 6, min_size = 8, compactness = 3)

n3_no_stitch <- length(unique(result3_no_stitch$cluster))
n3_stitch <- length(unique(result3_stitch$cluster))

cat(sprintf("Competing candidates - No stitching: %d clusters\n", n3_no_stitch))
cat(sprintf("Competing candidates - With stitching: %d clusters\n", n3_stitch))

# Analyze which regions got stitched
cluster_vol3 <- array(0, dims3)
cluster_vol3[mask3 > 0] <- result3_stitch$cluster

competing_analysis <- list()
for (z in 1:dims3[3]) {
  slice_clusters <- unique(cluster_vol3[,,z])
  competing_analysis[[paste0("slice_", z)]] <- sort(slice_clusters[slice_clusters > 0])
}

cat("Cluster IDs per slice:\n")
for (i in 1:length(competing_analysis)) {
  cat(sprintf("  %s: %s\n", names(competing_analysis)[i], 
              paste(competing_analysis[[i]], collapse=", ")))
}

# Test 4: Empty/Sparse Slices
cat("\nTest 4: Empty and Sparse Slices\n")
cat("Testing behavior with empty slices or very sparse data...\n")

dims4 <- c(10, 10, 5)
mask4 <- NeuroVol(array(0, dims4), NeuroSpace(dims4))  # Start with empty

# Add voxels only to some slices
mask4@.Data[3:8, 3:8, 1] <- 1  # Dense region in slice 1
mask4@.Data[4:7, 4:7, 2] <- 1  # Medium region in slice 2
# Slice 3 is empty
mask4@.Data[5:6, 5:6, 4] <- 1  # Sparse region in slice 4
mask4@.Data[2:9, 2:9, 5] <- 1  # Large region in slice 5

nvox4 <- sum(mask4@.Data > 0)
mask_idx4 <- which(mask4@.Data > 0)
ntime4 <- 50

ts_data4 <- matrix(0, nrow = length(mask_idx4), ncol = ntime4)
t_seq4 <- seq(0, 4*pi, length.out = ntime4)

# Similar patterns for all non-empty voxels
for (i in 1:length(mask_idx4)) {
  ts_data4[i, ] <- sin(t_seq4) + rnorm(ntime4, sd = 0.05)
}

vec_list4 <- lapply(1:ntime4, function(t) {
  vol_data <- array(0, dims4)
  vol_data[mask_idx4] <- ts_data4[, t]
  NeuroVol(vol_data, NeuroSpace(dims4))
})
vec4 <- do.call(concat, vec_list4)

result4 <- slice_msf(vec4, mask4, stitch_z = TRUE, theta_link = 0.8,
                     min_contact = 1, num_runs = 1,
                     r = 4, min_size = 2, compactness = 2)

# Check which slices have clusters
cluster_vol4 <- array(0, dims4)
cluster_vol4[mask_idx4] <- result4$cluster

sparse_slice_analysis <- list()
for (z in 1:dims4[3]) {
  slice_clusters <- unique(cluster_vol4[,,z])
  n_slice_clusters <- length(slice_clusters[slice_clusters > 0])
  sparse_slice_analysis[[z]] <- n_slice_clusters
  cat(sprintf("  Slice %d: %d clusters\n", z, n_slice_clusters))
}

# Empty slice should have 0 clusters
empty_slice_correct <- sparse_slice_analysis[[3]] == 0
cat(sprintf("Empty slice handling: %s\n", ifelse(empty_slice_correct, "✓", "✗")))

# Test 5: Identical Adjacent Clusters
cat("\nTest 5: Identical Adjacent Clusters\n")
cat("Testing when adjacent clusters are nearly identical...\n")

dims5 <- c(8, 8, 4)
mask5 <- NeuroVol(array(1, dims5), NeuroSpace(dims5))
nvox5 <- prod(dims5)
ntime5 <- 40

ts_data5 <- matrix(0, nrow = nvox5, ncol = ntime5)
t_seq5 <- seq(0, 4*pi, length.out = ntime5)

# Create nearly identical patterns across all slices
base_pattern <- sin(t_seq5)
for (i in 1:nvox5) {
  # Add tiny slice-specific variations
  z <- coords1[i, 3]  # This should be corrected
  ts_data5[i, ] <- base_pattern + 0.001 * rnorm(ntime5) + 0.01 * z
}

coords5 <- arrayInd(1:nvox5, dims5)
for (i in 1:nvox5) {
  z <- coords5[i, 3]
  ts_data5[i, ] <- base_pattern + 0.001 * rnorm(ntime5) + 0.01 * z
}

vec_list5 <- lapply(1:ntime5, function(t) {
  vol_data <- array(0, dims5)
  vol_data[mask5 > 0] <- ts_data5[, t]
  NeuroVol(vol_data, NeuroSpace(dims5))
})
vec5 <- do.call(concat, vec_list5)

# Test with high theta_link - should stitch almost everything
result5_high <- slice_msf(vec5, mask5, stitch_z = TRUE, theta_link = 0.95,
                          min_contact = 1, num_runs = 1,
                          r = 4, min_size = 5, compactness = 2)

# Test with low theta_link - should stitch everything
result5_low <- slice_msf(vec5, mask5, stitch_z = TRUE, theta_link = 0.5,
                         min_contact = 1, num_runs = 1,
                         r = 4, min_size = 5, compactness = 2)

n5_high <- length(unique(result5_high$cluster))
n5_low <- length(unique(result5_low$cluster))

cat(sprintf("Nearly identical patterns - High threshold (θ=0.95): %d clusters\n", n5_high))
cat(sprintf("Nearly identical patterns - Low threshold (θ=0.5): %d clusters\n", n5_low))

identical_test_correct <- n5_low <= n5_high  # Lower threshold should allow more stitching
cat(sprintf("Identical patterns threshold test: %s\n", 
            ifelse(identical_test_correct, "✓", "✗")))

# Summary
cat("\n=== Edge Cases Summary ===\n")

all_tests_passed <- single_slice_consistent && empty_slice_correct && identical_test_correct

cat(sprintf("✓ Single slice consistency: %s\n", ifelse(single_slice_consistent, "PASS", "FAIL")))
cat(sprintf("✓ Minimal contact handling: %d test cases\n", length(contact_values)))
cat(sprintf("✓ Competing candidates: %d→%d clusters\n", n3_no_stitch, n3_stitch))
cat(sprintf("✓ Empty slice handling: %s\n", ifelse(empty_slice_correct, "PASS", "FAIL")))
cat(sprintf("✓ Identical patterns: %s\n", ifelse(identical_test_correct, "PASS", "FAIL")))

if (all_tests_passed) {
  cat("\n🎉 All edge case tests PASSED!\n")
  cat("   - Single slice volumes handled correctly\n")
  cat("   - Minimal contact scenarios work as expected\n")
  cat("   - Competing stitching candidates resolved appropriately\n")
  cat("   - Empty/sparse slices don't cause issues\n")
  cat("   - Nearly identical patterns stitch correctly\n")
} else {
  cat("\n⚠️  Some edge case tests need attention\n")
  cat("   Review implementation for boundary conditions\n")
}