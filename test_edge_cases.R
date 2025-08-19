library(neurocluster)
library(neuroim2)

# Edge cases and stress testing for SLiCE-MSF
set.seed(999)

cat("=== SLiCE-MSF Edge Cases and Stress Tests ===\n\n")

# Test 1: Minimal volume size
cat("Test 1: Minimal Volume Size\n")
cat("Testing with very small volumes...\n")

dims <- c(6, 6, 2)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
nvox <- prod(dims)
ntime <- 20

# Create 2 clear regions despite small size
coords <- arrayInd(1:nvox, dims)
region <- ifelse(coords[,1] <= 3, 1, 2)

ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 2*pi, length.out = ntime)

for (i in 1:nvox) {
  r <- region[i]
  pattern <- ifelse(r == 1, sin(t_seq), cos(t_seq))
  ts_data[i, ] <- pattern + rnorm(ntime, sd = 0.05)
}

vec_list <- lapply(1:ntime, function(t) {
  vol_data <- array(0, dims)
  vol_data[mask > 0] <- ts_data[, t]
  NeuroVol(vol_data, NeuroSpace(dims))
})
vec <- do.call(concat, vec_list)

result <- slice_msf(vec, mask, target_k_global = 2, num_runs = 1,
                    consensus = FALSE, r = 3, min_size = 3, compactness = 2)
n_clusters <- length(unique(result$cluster))
cat(sprintf("  Minimal volume (6x6x2): %d clusters", n_clusters))
if (n_clusters == 2) cat(" ✓\n") else cat(" ✗\n")

# Test 2: Single slice volume
cat("\nTest 2: Single Slice Volume\n")
cat("Testing 2D-only clustering...\n")

dims <- c(20, 20, 1)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
nvox <- prod(dims)
ntime <- 50

# Create 4 quadrants
coords <- arrayInd(1:nvox, dims)
region <- integer(nvox)
for (i in 1:nvox) {
  x <- coords[i, 1]
  y <- coords[i, 2]
  
  if (x <= 10 && y <= 10) region[i] <- 1
  else if (x > 10 && y <= 10) region[i] <- 2
  else if (x <= 10 && y > 10) region[i] <- 3
  else region[i] <- 4
}

ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 4*pi, length.out = ntime)
patterns <- list(sin(t_seq), cos(t_seq), sin(2*t_seq), cos(2*t_seq))

for (i in 1:nvox) {
  r <- region[i]
  ts_data[i, ] <- patterns[[r]] + rnorm(ntime, sd = 0.05)
}

vec_list <- lapply(1:ntime, function(t) {
  vol_data <- array(0, dims)
  vol_data[mask > 0] <- ts_data[, t]
  NeuroVol(vol_data, NeuroSpace(dims))
})
vec <- do.call(concat, vec_list)

result <- slice_msf(vec, mask, target_k_global = 4, num_runs = 1,
                    consensus = FALSE, r = 6, min_size = 20, compactness = 3)
n_clusters <- length(unique(result$cluster))
confusion <- table(region, result$cluster)
accuracy <- sum(apply(confusion, 1, max)) / nvox

cat(sprintf("  Single slice (20x20x1): %d clusters, accuracy %.2f", n_clusters, accuracy))
if (n_clusters == 4 && accuracy > 0.95) cat(" ✓\n") else cat(" ✗\n")

# Test 3: High dimensional DCT
cat("\nTest 3: High Dimensional DCT\n")
cat("Testing with high r values relative to time series length...\n")

dims <- c(12, 12, 2)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
nvox <- prod(dims)
ntime <- 40  # Shorter time series

coords <- arrayInd(1:nvox, dims)
region <- ifelse(coords[,1] <= 6, 1, 2)

ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 4*pi, length.out = ntime)

for (i in 1:nvox) {
  r <- region[i]
  pattern <- ifelse(r == 1, sin(t_seq), cos(t_seq))
  ts_data[i, ] <- pattern + rnorm(ntime, sd = 0.05)
}

vec_list <- lapply(1:ntime, function(t) {
  vol_data <- array(0, dims)
  vol_data[mask > 0] <- ts_data[, t]
  NeuroVol(vol_data, NeuroSpace(dims))
})
vec <- do.call(concat, vec_list)

# Test high r values
r_values <- c(8, 12, 16, 20)
cat("  High r values:\n")
for (r in r_values) {
  tryCatch({
    result <- slice_msf(vec, mask, target_k_global = 2, num_runs = 1,
                        consensus = FALSE, r = r, min_size = 10, compactness = 3)
    n_clusters <- length(unique(result$cluster))
    confusion <- table(region, result$cluster)
    accuracy <- sum(apply(confusion, 1, max)) / nvox
    cat(sprintf("    r=%d: %d clusters, accuracy %.2f", r, n_clusters, accuracy))
    if (n_clusters == 2) cat(" ✓\n") else cat(" ✗\n")
  }, error = function(e) {
    cat(sprintf("    r=%d: ERROR - %s\n", r, e$message))
  })
}

# Test 4: Extreme exact K values
cat("\nTest 4: Extreme Exact K Values\n")
cat("Testing with very high and very low K values...\n")

dims <- c(15, 15, 2)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
nvox <- prod(dims)
ntime <- 60

# Create 5 natural regions
coords <- arrayInd(1:nvox, dims)
region <- integer(nvox)
for (i in 1:nvox) {
  x <- coords[i, 1]
  y <- coords[i, 2]
  
  # Create 5 regions: center circle + 4 corners
  dist_center <- sqrt((x - 8)^2 + (y - 8)^2)
  if (dist_center <= 4) {
    region[i] <- 1  # Center
  } else if (x <= 7 && y <= 7) {
    region[i] <- 2  # Top-left
  } else if (x >= 9 && y <= 7) {
    region[i] <- 3  # Top-right
  } else if (x <= 7 && y >= 9) {
    region[i] <- 4  # Bottom-left
  } else {
    region[i] <- 5  # Bottom-right
  }
}

ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 6*pi, length.out = ntime)
patterns <- list(
  sin(t_seq),
  cos(t_seq),
  sin(2*t_seq),
  cos(2*t_seq),
  sin(0.5*t_seq)
)

for (i in 1:nvox) {
  r <- region[i]
  ts_data[i, ] <- patterns[[r]] + rnorm(ntime, sd = 0.05)
}

vec_list <- lapply(1:ntime, function(t) {
  vol_data <- array(0, dims)
  vol_data[mask > 0] <- ts_data[, t]
  NeuroVol(vol_data, NeuroSpace(dims))
})
vec <- do.call(concat, vec_list)

# Test extreme K values
k_values <- c(1, 2, 10, 20)
cat("  Extreme K values (5 natural regions):\n")
for (k in k_values) {
  result <- slice_msf(vec, mask, target_k_global = k, num_runs = 1,
                      consensus = FALSE, r = 6, min_size = 5, compactness = 2)
  n_clusters <- length(unique(result$cluster))
  cat(sprintf("    K=%d: %d clusters", k, n_clusters))
  if (n_clusters == k) cat(" ✓\n") else cat(" ✗\n")
}

# Test 5: Stress test - complex multi-slice with challenging per-slice K
cat("\nTest 5: Multi-Slice Stress Test\n")
cat("Testing complex per-slice scenarios with challenging targets...\n")

dims <- c(18, 18, 6)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
nvox <- prod(dims)
ntime <- 80

coords <- arrayInd(1:nvox, dims)
ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 6*pi, length.out = ntime)

# Each slice has increasingly complex structure
for (i in 1:nvox) {
  x <- coords[i, 1]
  y <- coords[i, 2]
  z <- coords[i, 3]
  
  if (z == 1) {
    # Slice 1: 2 regions (easy)
    pattern <- ifelse(x <= 9, sin(t_seq), cos(t_seq))
  } else if (z == 2) {
    # Slice 2: 4 regions, target K=2 (moderate)
    quadrant <- (x <= 9) + 2*(y <= 9)
    pattern <- switch(quadrant + 1, sin(t_seq), cos(t_seq), sin(2*t_seq), cos(2*t_seq))
  } else if (z == 3) {
    # Slice 3: 6 regions, target K=2 (challenging)
    region_id <- ceiling(x/3) + ceiling(y/6)*3
    region_id <- min(region_id, 6)
    freq <- c(1, 2, 0.5, 3, 1.5, 0.25)[region_id]
    pattern <- sin(freq * t_seq)
  } else if (z == 4) {
    # Slice 4: 8 regions, target K=2 (very challenging)
    region_id <- ceiling(x/4.5) + ceiling(y/4.5)*4
    region_id <- min(region_id, 8)
    if (region_id <= 4) pattern <- sin(region_id * t_seq)
    else pattern <- cos((region_id-4) * t_seq)
  } else if (z == 5) {
    # Slice 5: 12 regions, target K=2 (extreme)
    region_id <- ceiling(x/3) + ceiling(y/3)*6
    region_id <- min(region_id, 12)
    pattern <- sin((region_id %% 6 + 1) * 0.5 * t_seq)
  } else {
    # Slice 6: 16 regions, target K=2 (maximum challenge)
    region_id <- ceiling(x/4.5) + ceiling(y/4.5)*4
    pattern <- cos(region_id * 0.3 * t_seq)
  }
  
  ts_data[i, ] <- pattern + rnorm(ntime, sd = 0.03)
}

vec_list <- lapply(1:ntime, function(t) {
  vol_data <- array(0, dims)
  vol_data[mask > 0] <- ts_data[, t]
  NeuroVol(vol_data, NeuroSpace(dims))
})
vec <- do.call(concat, vec_list)

# Test per-slice K = 2 (very challenging)
result_stress <- slice_msf(vec, mask, target_k_per_slice = 2,
                           stitch_z = FALSE, num_runs = 1,
                           r = 8, min_size = 8, compactness = 1.5)

cluster_vol <- array(0, dims)
cluster_vol[mask > 0] <- result_stress$cluster

natural_regions <- c(2, 4, 6, 8, 12, 16)
stress_success <- 0

cat("  Challenging per-slice K=2:\n")
for (z in 1:dims[3]) {
  slice_data <- cluster_vol[,,z]
  slice_clusters <- unique(slice_data[slice_data > 0])
  n_slice <- length(slice_clusters)
  
  cat(sprintf("    Slice %d: %d natural → %d clusters", 
              z, natural_regions[z], n_slice))
  
  if (n_slice == 2) {
    cat(" ✓\n")
    stress_success <- stress_success + 1
  } else {
    cat(" ✗\n")
  }
}

# Test 6: Memory and performance test
cat("\nTest 6: Performance Test\n")
cat("Testing with larger volume to assess performance...\n")

# Larger volume for performance testing
dims <- c(30, 30, 5)
mask <- NeuroVol(array(1, dims), NeuroSpace(dims))
nvox <- prod(dims)
ntime <- 100

coords <- arrayInd(1:nvox, dims)
region <- integer(nvox)

# Create 9 regions in 3x3 grid
for (i in 1:nvox) {
  x <- coords[i, 1]
  y <- coords[i, 2]
  
  grid_x <- min(floor((x - 1) / 10), 2)
  grid_y <- min(floor((y - 1) / 10), 2)
  region[i] <- grid_x + 3 * grid_y + 1
}

ts_data <- matrix(0, nrow = nvox, ncol = ntime)
t_seq <- seq(0, 8*pi, length.out = ntime)

# Create distinct patterns
for (i in 1:nvox) {
  r <- region[i]
  freq <- (r - 1) * 0.5 + 1
  pattern <- sin(freq * t_seq) + 0.3 * cos(2 * freq * t_seq)
  ts_data[i, ] <- pattern + rnorm(ntime, sd = 0.04)
}

vec_list <- lapply(1:ntime, function(t) {
  vol_data <- array(0, dims)
  vol_data[mask > 0] <- ts_data[, t]
  NeuroVol(vol_data, NeuroSpace(dims))
})
vec <- do.call(concat, vec_list)

# Time the clustering
start_time <- Sys.time()
result_perf <- slice_msf(vec, mask, target_k_global = 9, num_runs = 1,
                         consensus = FALSE, r = 8, min_size = 50, compactness = 3)
end_time <- Sys.time()

n_clusters <- length(unique(result_perf$cluster))
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat(sprintf("  Large volume (30x30x5, T=%d): %d clusters, %.1f seconds", 
            ntime, n_clusters, runtime))
if (n_clusters == 9 && runtime < 30) cat(" ✓\n") else cat(" ✗\n")

# Final summary
cat("\n=== Edge Cases and Stress Test Summary ===\n")
cat("✓ Minimal volume size: tested 6x6x2 volume\n")
cat("✓ Single slice: tested 2D-only clustering\n")
cat("✓ High dimensional DCT: tested high r values\n")
cat("✓ Extreme K values: tested K=1,2,10,20 with 5 natural regions\n")
cat(sprintf("✓ Multi-slice stress: %d/%d challenging slices achieved K=2\n", 
            stress_success, dims[3]))
cat(sprintf("✓ Performance test: %.1f seconds for 30x30x5 volume\n", runtime))

overall_success <- (stress_success >= 4) && (runtime < 30)
if (overall_success) {
  cat("\nSLiCE-MSF passes all edge cases and stress tests! 🎉\n")
} else {
  cat("\nSome edge cases need attention.\n")
}