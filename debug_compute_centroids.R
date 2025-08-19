library(neurocluster)
library(neuroim2)
library(purrr)

# Create simple test data
feature_mat <- matrix(rnorm(100), nrow=10, ncol=10)
grid <- matrix(1:30, nrow=10, ncol=3)  # 10 spatial points with x,y,z coords  
assignment <- c(1,1,2,2,2,3,3,3,3,3)   # 3 clusters

# Test compute_centroids
centroids <- neurocluster::compute_centroids(feature_mat, grid, assignment, medoid=FALSE)

print("Structure of centroids:")
print(str(centroids))
print("Names of centroids:")
print(names(centroids))
print("Length of centroids:")
print(length(centroids))

if ("center" %in% names(centroids)) {
  print("centroids$center dim:")
  print(dim(centroids$center))
}

if ("centroid" %in% names(centroids)) {
  print("centroids$centroid dim:")
  print(dim(centroids$centroid))
}

# Test what lapply(centroids, `[[`, "center") produces
print("Testing lapply extraction:")
tryCatch({
  center_extract <- lapply(centroids, `[[`, "center")
  print("Center extract length:")
  print(length(center_extract))
  print("Center extract first element:")
  print(center_extract[[1]])
}, error = function(e) {
  print("Error in lapply extraction:")
  print(e)
})