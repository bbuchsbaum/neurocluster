pkgname <- "neurocluster"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('neurocluster')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("commute_cluster")
### * commute_cluster

flush(stderr()); flush(stdout())

### Name: commute_cluster
### Title: Commute Time Clustering
### Aliases: commute_cluster

### ** Examples

mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
vec <- replicate(10, NeuroVol(array(runif(202020), c(20,20,20)),
NeuroSpace(c(20,20,20))), simplify=FALSE)
vec <- do.call(concat, vec)

commute_res <- commute_cluster(vec, mask, K=100)




cleanEx()
nameEx("compute_centroids")
### * compute_centroids

flush(stderr()); flush(stdout())

### Name: compute_centroids
### Title: Compute Centroids of Clusters
### Aliases: compute_centroids

### ** Examples

## Not run: 
##D   # Assuming `feature_mat`, `grid`, and `assignment` are available
##D   centroids <- compute_centroids(feature_mat, grid, assignment)
##D   # To compute medoids instead of means
##D   medoids <- compute_centroids(feature_mat, grid, assignment, medoid=TRUE)
## End(Not run)




cleanEx()
nameEx("hello")
### * hello

flush(stderr()); flush(stdout())

### Name: hello
### Title: Hello, World!
### Aliases: hello

### ** Examples

hello()



cleanEx()
nameEx("knn_shrink")
### * knn_shrink

flush(stderr()); flush(stdout())

### Name: knn_shrink
### Title: K-nearest-neighbor shrink
### Aliases: knn_shrink

### ** Examples

mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
bvec <- replicate(10,
                  NeuroVol(array(runif(20*20*20), c(20,20,20)),
                           NeuroSpace(c(20,20,20))),
                  simplify=FALSE)
bvec <- do.call(concat, bvec)

sbvec <- knn_shrink(bvec, mask, k=3)




cleanEx()
nameEx("merge_clus")
### * merge_clus

flush(stderr()); flush(stdout())

### Name: merge_clus
### Title: Merge Clustering Results Using a Consensus Clustering Algorithm
### Aliases: merge_clus

### ** Examples

# Assuming clustering1, clustering2, and clustering3 are objects of class "cluster_result"
merged_clustering <- merge_clus(clustering1, clustering2, clustering3, method="SE")




cleanEx()
nameEx("slic4d_supervoxels")
### * slic4d_supervoxels

flush(stderr()); flush(stdout())

### Name: slic4d_supervoxels
### Title: Fast 4D SLIC supervoxels (mask-aware, gradient relocation,
###   preserve-K)
### Aliases: slic4d_supervoxels

### ** Examples

## Not run: 
##D   library(neuroim2)
##D   # Basic usage
##D   res <- slic4d_supervoxels(bvec, mask, K = 1000, compactness = 15)
##D   
##D   # With mask-aware seeding and gradient relocation
##D   res <- slic4d_supervoxels(bvec, mask, K = 1000, 
##D                            seed_method = "mask_poisson",
##D                            seed_relocate = "correlation",
##D                            preserve_k = TRUE)
## End(Not run)



cleanEx()
nameEx("slice_msf")
### * slice_msf

flush(stderr()); flush(stdout())

### Name: slice_msf
### Title: SLiCE-MSF: Slice-wise, Low-rank, Minimum-Spanning Forest
###   Clustering
### Aliases: slice_msf

### ** Examples

## Not run: 
##D # Load example data
##D mask <- NeuroVol(array(1, c(64,64,32)), NeuroSpace(c(64,64,32)))
##D vec <- replicate(100, NeuroVol(array(rnorm(64*64*32), c(64,64,32)),
##D                  NeuroSpace(c(64,64,32))), simplify=FALSE)
##D vec <- do.call(concat, vec)
##D 
##D # Example 1: Basic usage (may show z-plane artifacts)
##D result_basic <- slice_msf(vec, mask, 
##D                           num_runs = 1,
##D                           compactness = 5)
##D 
##D # Example 2: Reduced z-artifacts with aggressive stitching
##D # Recommended for minimizing slice boundaries
##D result_smooth <- slice_msf(vec, mask, 
##D                            theta_link = 0.75,      # More aggressive stitching
##D                            min_contact = 2,        # Moderate contact requirement
##D                            compactness = 3,        # Lower compactness
##D                            num_runs = 3,           # Multiple runs
##D                            consensus = TRUE,       # Enable consensus
##D                            use_features = TRUE,    # Feature-based consensus
##D                            lambda = 0.6)           # Balance features/labels
##D 
##D # Example 3: Conservative approach for anatomical boundaries
##D # Preserves natural boundaries but may show more z-artifacts
##D result_conservative <- slice_msf(vec, mask,
##D                                  theta_link = 0.90,    # Conservative stitching
##D                                  min_contact = 5,      # Strict contact requirement
##D                                  compactness = 7,      # Compact clusters
##D                                  min_size = 120)       # Larger minimum size
##D 
##D # Example 4: Exact K targeting with 100 clusters
##D result_exact <- slice_msf(vec, mask, 
##D                          target_k_global = 100,   # Exactly 100 clusters
##D                          use_features = TRUE,     # Required for exact K
##D                          num_runs = 3,
##D                          consensus = TRUE)
##D 
##D # Example 5: Per-slice consistency (useful for group studies)
##D result_per_slice <- slice_msf(vec, mask,
##D                               target_k_per_slice = 50,  # 50 clusters per slice
##D                               stitch_z = FALSE,         # No z-stitching
##D                               compactness = 5)
##D 
##D # Example 6: High-resolution data optimization
##D # For data with voxel size < 2mm
##D result_highres <- slice_msf(vec, mask,
##D                             min_size = 250,         # Larger clusters for high-res
##D                             r = 15,                 # More DCT components
##D                             compactness = 4,
##D                             theta_link = 0.78,
##D                             num_runs = 5,
##D                             consensus = TRUE)
##D 
##D # Example 7: Comparison with true 3D algorithm
##D # If z-artifacts are unacceptable, use supervoxels instead
##D result_3d <- supervoxels(vec, mask, n_supvox = 500, alpha = 0.5)
##D # supervoxels provides smooth 3D parcels without slice artifacts
## End(Not run)




cleanEx()
nameEx("snic")
### * snic

flush(stderr()); flush(stdout())

### Name: snic
### Title: SNIC: Simple Non-Iterative Clustering
### Aliases: snic

### ** Examples

mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
vec <- replicate(10, NeuroVol(array(runif(202020), c(20,20,20)),
NeuroSpace(c(20,20,20))), simplify=FALSE)
vec <- do.call(concat, vec)

snic_res <- snic(vec, mask, compactness=5, K=100)




cleanEx()
nameEx("spatial_gradient")
### * spatial_gradient

flush(stderr()); flush(stdout())

### Name: spatial_gradient
### Title: Spatial Gradient Calculation
### Aliases: spatial_gradient

### ** Examples

mask <- NeuroVol(array(1, c(20,20,20)), NeuroSpace(c(20,20,20)))
input_vol <- NeuroVol(array(runif(202020), c(20,20,20)),
NeuroSpace(c(20,20,20)))

gradient_vol <- spatial_gradient(input_vol, mask)




cleanEx()
nameEx("supervoxel_cluster_time")
### * supervoxel_cluster_time

flush(stderr()); flush(stdout())

### Name: supervoxel_cluster_time
### Title: Supervoxel Clustering in Time
### Aliases: supervoxel_cluster_time

### ** Examples

feature_mat <- matrix(rnorm(100 * 10), 100, 10)
library(future)
plan(multicore)
cres <- supervoxel_cluster_time(t(feature_mat), K=20)



cleanEx()
nameEx("supervoxels_flash3d")
### * supervoxels_flash3d

flush(stderr()); flush(stdout())

### Name: supervoxels_flash3d
### Title: FLASH-3D: Fast 3D Superclustering for fMRI
### Aliases: supervoxels_flash3d

### ** Examples

## Not run: 
##D   # Basic usage
##D   result <- supervoxels_flash3d(vec, mask, K = 100)
##D   
##D   # With custom parameters
##D   result <- supervoxels_flash3d(vec, mask, K = 100, 
##D                                 lambda_s = 0.8, lambda_t = 1.2,
##D                                 bits = 128, dctM = 16)
##D   
##D   # With barrier for anatomy-aware clustering
##D   barrier_vol <- create_anatomical_barrier(mask)
##D   result <- supervoxels_flash3d(vec, mask, K = 100,
##D                                 lambda_g = 0.5, barrier = barrier_vol)
## End(Not run)




cleanEx()
nameEx("tesselate")
### * tesselate

flush(stderr()); flush(stdout())

### Name: tesselate
### Title: Tesselate a Mask Volume into K Clusters using K-means
### Aliases: tesselate

### ** Examples

# Assuming you have a NeuroVol object 'mask' and you want to create 150 clusters
clustered_volume <- tesselate(mask, K = 150)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
