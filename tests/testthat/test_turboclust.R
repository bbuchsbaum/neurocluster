
library(neuroim2)
library(purrr)
library(assertthat)
library(neurocluster)
library(locfit)

bvec <- read_vec("testdata/rscan01.nii.gz")
mask <- read_vol("testdata/mask.nii")

test_that("can run turbo_clust on a NeuroVec", {


  cres <- turbo_cluster(bvec, mask, K=100, sigma1=1, sigma2=4, filter=list(lp=.03, hp=.33), filter_method="bspline")
})

test_that("can run turbo_clust on a NeuroVec with nreps=5", {

  cres <- turbo_cluster(bvec, mask, K=100, sigma1=1, sigma2=1.5, filter=list(lp=.05, hp=.33), filter_method="bspline",
                        sample_frac=.3, nreps=12)
})

test_that("can bandpass filter a NeuroVec using bsplines", {
  sbvec <- filter_vec(bvec, mask, hp=.8, lp=.04, method="bspline")
})

test_that("can knn_shrink a NeuroVec", {
  sbvec <- knn_shrink(bvec, mask, k=4)
})

test_that("can filter, then shrink, then compute gradient of a NeuroVec", {
  sbvec <- filter_vec(bvec, mask, hp=.8, lp=.03, method="bspline")
  sbvec <- knn_shrink(sbvec, mask, k=4)
  grad <- correlation_gradient(sbvec, mask)
})

test_that("meta_clust a turbo_cluster result", {
  cres <- turbo_cluster(bvec, mask, K=100, sigma1=1, sigma2=4, filter=list(lp=.03, hp=.33), filter_method="bspline")
})

test_that("commute_time cluster", {
  cres <- commute_cluster(bvec,mask,K=64, ncomp=100, sigma1=.73, sigma2=6,alpha=0,
                                         filter=list(lp=.03, hp=.33), filter_method="bspline")
})

