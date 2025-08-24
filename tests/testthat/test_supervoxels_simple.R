library(testthat)
library(neurocluster)
library(neuroim2)

test_that("basic supervoxels works without parallel", {
  set.seed(123)
  
  # Create small test data
  mask <- NeuroVol(array(1, c(10,10,10)), NeuroSpace(c(10,10,10), c(1,1,1)))
  
  # Create random data
  data_arr <- array(rnorm(10*10*10*5), c(10,10,10,5))
  vec_list <- lapply(1:5, function(i) {
    NeuroVol(data_arr[,,,i], NeuroSpace(c(10,10,10), c(1,1,1)))
  })
  bvec <- do.call(concat, vec_list)
  
  # Run with parallel disabled
  res <- supervoxels(bvec, mask, K=10, iterations=2, parallel=FALSE)
  
  # Check basic properties
  expect_equal(length(res$cluster), 1000)  # 10x10x10
  expect_true(all(res$cluster %in% 1:10))
  expect_equal(nrow(res$centers), 10)
  expect_equal(nrow(res$coord_centers), 10)
})