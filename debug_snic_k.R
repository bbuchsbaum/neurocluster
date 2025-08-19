library(neurocluster)
library(neuroim2)

# Function from test file
create_test_data <- function(dims = c(10, 10, 10), nvols = 5) {
  # Create a simple mask
  mask_array <- array(1, dims)
  # Remove border to ensure reasonable connectivity
  mask_array[c(1, dims[1]), , ] <- 0
  mask_array[, c(1, dims[2]), ] <- 0
  mask_array[, , c(1, dims[3])] <- 0
  
  mask <- NeuroVol(mask_array, NeuroSpace(dims))
  
  # Create 4D array for NeuroVec (x, y, z, time)
  data_4d <- array(0, c(dims, nvols))
  
  # Add spatial patterns - create distinct regions
  nx <- dims[1]
  ny <- dims[2] 
  nz <- dims[3]
  
  # Region 1: top-left-front
  data_4d[1:(nx/2), 1:(ny/2), 1:(nz/2), ] <- rnorm(length(1:(nx/2)) * length(1:(ny/2)) * length(1:(nz/2)) * nvols, mean = 10, sd = 1)
  
  # Region 2: bottom-right-back  
  data_4d[(nx/2+1):nx, (ny/2+1):ny, (nz/2+1):nz, ] <- rnorm(length((nx/2+1):nx) * length((ny/2+1):ny) * length((nz/2+1):nz) * nvols, mean = -10, sd = 1)
  
  # Add some noise
  data_4d <- data_4d + rnorm(length(data_4d), sd = 0.5)
  
  # Apply mask
  for (t in 1:nvols) {
    data_4d[,,, t] <- data_4d[,,, t] * mask_array
  }
  
  # Create NeuroVec
  space_4d <- NeuroSpace(c(dims, nvols))
  bvec <- NeuroVec(data_4d, space_4d)
  
  list(bvec = bvec, mask = mask)
}

# Test the failing case
test_data <- create_test_data(dims = c(8, 8, 8), nvols = 2)
mask_indices <- which(test_data$mask > 0)

print(paste("Number of voxels in mask:", length(mask_indices)))

# Test with K=5
result_small <- snic(test_data$bvec, test_data$mask, K = 5, compactness = 5)
n_clusters <- length(unique(result_small$cluster))

print(paste("Requested K=5, got", n_clusters, "clusters"))
print(paste("Unique cluster IDs:", paste(sort(unique(result_small$cluster)), collapse=", ")))

# Check cluster sizes
cluster_sizes <- table(result_small$cluster)
print("Cluster sizes:")
print(cluster_sizes)

# Test if some clusters might be very small (isolated voxels)
small_clusters <- sum(cluster_sizes < 5)
print(paste("Number of clusters with < 5 voxels:", small_clusters))