
#' @importFrom rflann Neighbour
#' @importFrom hashmap hashmap
patch_cluster <- function(bvec, mask, K=500, patch_radius=8, connectivity=27,
                          knn=5, filter=list(lp=0, hp=0)) {

  mask.idx <- which(mask > 0)
  grid <- index_to_coord(mask, mask.idx)
  vgrid <- index_to_grid(mask, mask.idx)

  valmat <- series(bvec, mask.idx)

  if (any(filter > 0)) {
    message("patch_cluster: pre-filtering time series")
    valmat <- filter_mat(valmat, filter$lp, filter$hp)
  }

  slight <- neuroim2::searchlight(mask, radius=patch_radius, eager=TRUE)
  clist <- lapply(slight, identity)

  nn <- 50
  hmap <- hashmap("-1","-1")
  nabes <- rflann::Neighbour(grid, grid, nn)

  Rmat <-  Matrix(0, nrow = length(clist), ncol = length(clist), sparse = TRUE)

  lapply(1:length(clist), function(i) {
    print(i)
    s1 <- series(bvec, clist[[i]])

    jind <- nabes$indices[i,2:nn]
    p1 <- pmin(i, jind)
    p2 <- pmax(i, jind)
    knames <- paste(p1, p2, sep="-")

    for (j in jind) {
      if (Rmat[i,j] == 0) {
        Rmat[i,j] <- MatrixCorrelation::RV2(s1, series(bvec, clist[[j]]))
      }
    }

  })

}
