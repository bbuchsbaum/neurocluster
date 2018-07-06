
library(neuroim2)
devtools::load_all()
library(locfit)

test_that("can run turbo_clust on a NeuroVec", {
  bvec <- read_vec("testdata/rscan01.nii.gz")
  mask <- read_vol("testdata/mask.nii")

  cres <- turbo_cluster(bvec, mask, K=100, sigma1=1, sigma2=4, filter=list(lp=.05, hp=.33), filter_method="bspline")
})

test_that("can run turbo_clust on a NeuroVec with nreps=5", {

  cres <- turbo_cluster(bvec, mask, K=100, sigma1=1, sigma2=1.5, filter=list(lp=.05, hp=.33), filter_method="bspline",
                        sample_frac=.3, nreps=12)
})
