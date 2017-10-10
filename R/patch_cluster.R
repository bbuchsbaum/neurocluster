
#' @importFrom rflann Neighbour
#' @importFrom hashmap hashmap
patch_cluster <- function(bvec, mask, K=500, patch_radius=8, connectivity=27,
                          knn=5, filter=list(lp=-1, hp=-1)) {

  mask.idx <- which(mask > 0)
  grid <- indexToCoord(mask, mask.idx)
  vgrid <- indexToGrid(mask, mask.idx)

  valmat <- series(bvec, mask.idx)

  if (any(filter) > 0) {
    message("patch_cluster: filtering time series")
    valmat <- filter_mat(valmat, filter$lp, filter$hp)
  }

  slight <- neuroim::Searchlight(mask, radius=patch_radius, eager=TRUE)
  clist <- lapply(slight, identity)

  nn <- 27
  hmap <- hashmap("-1","-1")
  nabes <- rflann::Neighbour(grid, grid, nn)
  lapply(1:length(clist), function(i) {
    print(i)
    s1 <- series(bvec, clist[[i]])


    jind <- nabes$indices[i,2:nn]
    p1 <- pmin(i, jind)
    p2 <- pmax(i, jind)
    knames <- paste(p1, p2, sep="-")

    rvs <- sapply(1:length(jind), function(j) {
      if (hmap$has_key(knames[j])) {
        hmap[[knames[[j]]]]
      } else {
        MatrixCorrelation::RV2(s1, series(bvec, clist[[jind[j]]]))
      }
    })

    hmap[[knames]] <- rvs


  })






}
