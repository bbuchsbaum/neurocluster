
library(testthat)
library(neurocluster)
library(neuroim2)
library(purrr)
library(assertthat)

# Create simulated data instead of reading files
# This makes tests more portable and avoids path issues
dims <- c(32, 32, 20)
mask_array <- array(FALSE, dims)
center <- dims / 2
for (i in 1:dims[1]) {
  for (j in 1:dims[2]) {
    for (k in 1:dims[3]) {
      if (sum(((c(i,j,k) - center) / (dims/3))^2) <= 1) {
        mask_array[i,j,k] <- TRUE
      }
    }
  }
}

mask <- NeuroVol(mask_array, NeuroSpace(dims, c(3,3,3)))

# Simulate fMRI data using the new simulate_fmri function from neuroim2
bvec <- neuroim2::simulate_fmri(mask, n_time = 100, seed = 42)

test_that("can run supervoxels on a NeuroVec", {
  cres <- supervoxels(bvec, mask, K=100,
                     sigma1=1, sigma2=4,
                     iterations=10)  # Reduced iterations for faster testing
  expect_s3_class(cres, "cluster_result")
  expect_s4_class(cres$clusvol, "ClusteredNeuroVol")
})

test_that("can run tesselate on a mask", {
  kvol <- tesselate(mask, K=50)
  expect_s4_class(kvol, "ClusteredNeuroVol")
})

test_that("can knn_shrink a NeuroVec", {
  sbvec <- knn_shrink(bvec, mask, k=4)
  expect_s4_class(sbvec, "SparseNeuroVec")
})

test_that("can compute spatial gradient", {
  # Get first volume for gradient computation
  first_vol <- vols(bvec)[[1]]
  grad <- spatial_gradient(first_vol, mask)
  expect_s4_class(grad, "NeuroVol")
})

test_that("can run meta_clust on cluster result", {
  # First create a cluster result
  cres <- supervoxels(bvec, mask, K=50, iterations=5)
  
  # Then run meta clustering
  meta_result <- meta_clust(cres, cuts=2)
  expect_type(meta_result, "list")
  expect_true("cvols" %in% names(meta_result))
  expect_true("hclus" %in% names(meta_result))
})

test_that("can run commute_cluster", {
  # Skip if dependencies not available
  skip_if_not_installed("igraph")
  
  cres <- commute_cluster(bvec, mask, K=20, ncomp=10, 
                         sigma1=.73, sigma2=6, alpha=0.5)
  expect_s3_class(cres, "cluster_result")
  expect_s4_class(cres$clusvol, "ClusteredNeuroVol")
})

