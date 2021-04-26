

# @import neighborweights
# markov_cluster_fit <- function(feature_mat, connectivity=5, inflation=2, expansion=2) {
#   A <- edge_weights(feature_mat, neighbor_mode = "knn",
#                     weight_mode="heat", type="normal", k=connectivity, sigma=.73)
#   D <- colSums(A)
#   An <- A %*% Diagonal(x=1/D)
#
#   An <- cbind(
#     c(.5,.5,0,0,0,0,0),
#     c(.333,.333,.333,0,0,0,0),
#     c(0,.333,.333,.333,0,0,0),
#     c(0,0,.333,.333,.333,0,0),
#     c(0,0,0,.333,.333,.333,0),
#     c(0,0,0,0,.333,.333,.333),
#     c(0,0,0,0,0,.5,.5)
#   )

  An <- as(An, "sparseMatrix")
  iter = 1
  while (iter < 10) {
    print(nnzero(An))
    print(iter)
    An <- An %*% An
    An <- An * An

    D <- colSums(An)
    conv = sum(D) - nrow(An)
    print(conv)
    An <- An %*% Diagonal(x=1/D)

    aind <- which(An > 0, arr.ind=TRUE)
    t <- .01*mean(An[aind])
    message("threshold: ", t)
    An[which(An > 0 & An < t)] <- 0
    print(An)
    iter <- iter+1
  }
}
