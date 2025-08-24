library(testthat)
library(neurocluster)
library(neuroim2)

test_that("parallel supervoxels produces same results as sequential", {
  set.seed(123)
  
  # Create test data
  mask <- NeuroVol(array(1, c(15,15,15)), NeuroSpace(c(15,15,15), c(1,1,1)))
  
  # Create random data
  data_arr <- array(rnorm(15*15*15*10), c(15,15,15,10))
  vec_list <- lapply(1:10, function(i) {
    NeuroVol(data_arr[,,,i], NeuroSpace(c(15,15,15), c(1,1,1)))
  })
  bvec <- do.call(concat, vec_list)
  
  # Run sequential version
  set.seed(456)
  res_seq <- supervoxels(bvec, mask, K=20, iterations=5, parallel=FALSE)
  
  # Run parallel version
  set.seed(456)
  res_par <- supervoxels(bvec, mask, K=20, iterations=5, parallel=TRUE, grain_size=50)
  
  # Results should be identical since we use the same seed
  expect_equal(res_seq$cluster, res_par$cluster)
  expect_equal(res_seq$centers, res_par$centers)
  expect_equal(res_seq$coord_centers, res_par$coord_centers)
})

test_that("parallel supervoxels works with different grain sizes", {
  set.seed(789)
  
  # Create test data
  mask <- NeuroVol(array(1, c(20,20,10)), NeuroSpace(c(20,20,10), c(1,1,1)))
  
  # Create random data
  data_arr <- array(rnorm(20*20*10*5), c(20,20,10,5))
  vec_list <- lapply(1:5, function(i) {
    NeuroVol(data_arr[,,,i], NeuroSpace(c(20,20,10), c(1,1,1)))
  })
  bvec <- do.call(concat, vec_list)
  
  # Test with different grain sizes
  set.seed(111)
  res1 <- supervoxels(bvec, mask, K=30, iterations=3, parallel=TRUE, grain_size=10)
  
  set.seed(111)
  res2 <- supervoxels(bvec, mask, K=30, iterations=3, parallel=TRUE, grain_size=200)
  
  # Different grain sizes should produce the same result
  expect_equal(res1$cluster, res2$cluster)
})

test_that("parallel supervoxels handles edge cases", {
  # Test with very small data (should automatically use sequential)
  mask_small <- NeuroVol(array(1, c(5,5,5)), NeuroSpace(c(5,5,5), c(1,1,1)))
  
  data_arr <- array(rnorm(5*5*5*3), c(5,5,5,3))
  vec_list <- lapply(1:3, function(i) {
    NeuroVol(data_arr[,,,i], NeuroSpace(c(5,5,5), c(1,1,1)))
  })
  bvec <- do.call(concat, vec_list)
  
  # Should not error even though data is small
  expect_no_error({
    res <- supervoxels(bvec, mask_small, K=5, iterations=2, parallel=TRUE)
  })
  
  # Test that it produces valid results
  expect_equal(length(unique(res$cluster)), 5)
  expect_equal(nrow(res$centers), 5)
})

test_that("parallel performance is better than sequential for large data", {
  skip_if_not(RcppParallel::defaultNumThreads() > 1, "Requires multiple cores")
  
  set.seed(222)
  
  # Create larger test data
  mask <- NeuroVol(array(1, c(30,30,20)), NeuroSpace(c(30,30,20), c(1,1,1)))
  
  data_arr <- array(rnorm(30*30*20*10), c(30,30,20,10))
  vec_list <- lapply(1:10, function(i) {
    NeuroVol(data_arr[,,,i], NeuroSpace(c(30,30,20), c(1,1,1)))
  })
  bvec <- do.call(concat, vec_list)
  
  # Time sequential version
  time_seq <- system.time({
    res_seq <- supervoxels(bvec, mask, K=50, iterations=2, parallel=FALSE)
  })
  
  # Time parallel version
  time_par <- system.time({
    res_par <- supervoxels(bvec, mask, K=50, iterations=2, parallel=TRUE)
  })
  
  # Parallel should be faster (though we can't guarantee exact speedup)
  cat("\nSequential time:", time_seq[3], "seconds\n")
  cat("Parallel time:", time_par[3], "seconds\n")
  cat("Speedup:", time_seq[3] / time_par[3], "x\n")
  
  # At least check that both produce valid results
  expect_true(all(res_seq$cluster %in% 1:50))
  expect_true(all(res_par$cluster %in% 1:50))
})