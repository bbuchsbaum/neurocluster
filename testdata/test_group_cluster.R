library(neuroim2)
library(dplyr)
library(purrr)

path <- "/Users/bbuchsbaum/Dropbox/code/group_cluster/"

fnames <- paste0(c("sub-2005", "sub-2002", "sub-2003", "sub-2004"), "_all_iresp.nii.gz")
mask <- read_vol(paste0(path, "gray_mask_small.nii"))

initclus <- tesselate(mask, 300)
initmeans <- lapply(fnames, function(fn) {
  p <- paste0(path, "/", fn)
  vec <- read_vec(p)
  mask %>% split_clusters(initclus) %>% map(~ rowMeans(series(vec, coords(.)))) %>% bind_rows()
}
