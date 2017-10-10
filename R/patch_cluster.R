
#' @importFrom rflann Neighbour
#' @importFrom hashmap hashmap
patch_cluster <- function(bvec, mask, K=500, radius=8, knn=5, filter=list(lp=-1, hp=-1)) {

  mask.idx <- which(mask > 0)
  grid <- indexToCoord(mask, mask.idx)
  vgrid <- indexToGrid(mask, mask.idx)

  valmat <- series(bvec, mask.idx)

  if (any(filter) > 0) {
    message("turbo_cluster: filtering time series")
    valmat <- filter_mat(valmat, filter$lp, filter$hp)
  }

  slight <- neuroim::Searchlight(mask, radius=radius)
  clist <- foreach(x=slight) %do% {
    x
  }


  nn <- 27
  hmap <- hashmap("-1","-1")
  nabes <- Neighbour(grid, grid, nn)
  lapply(1:length(clist), function(i) {
    print(i)
    s1 <- series(bvec, clist[[i]])
    sapply(nabes$indices[i,2:nn], function(j) {
      key <- paste0(sort(c(i,j)), collapse="-")
      if (hmap$has_key(key)) {
        #message("has ", key)
        hmap[[key]]
      } else {
        rv <- MatrixCorrelation::RV2(s1, series(bvec, clist[[j]]))
        hmap[[key]] <- rv
        rv
      }
    })
  })






}
