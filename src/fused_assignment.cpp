// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <unordered_set>
#include <set>
#include <map>
#include <cmath>
#include <atomic>

using namespace Rcpp;
using namespace RcppParallel;

// Optimized heat kernel computation
inline double heat_kernel_opt(const double* x1, const double* x2, int len, double sigma) {
  double dist_sq = 0.0;
  for (int i = 0; i < len; i++) {
    double diff = x1[i] - x2[i];
    dist_sq += diff * diff;
  }
  return std::exp(-dist_sq / (2.0 * sigma * sigma));
}

inline double normalized_heat_kernel_opt(const double* x1, const double* x2, int len, double sigma) {
  double dist_sq = 0.0;
  for (int i = 0; i < len; i++) {
    double diff = x1[i] - x2[i];
    dist_sq += diff * diff;
  }
  double norm_dist = dist_sq / (2.0 * len);
  return std::exp(-norm_dist / (2.0 * sigma * sigma));
}

// Fused operation: Find candidates and compute best assignment in one pass
// This eliminates the huge candlist allocation
// [[Rcpp::export]]
IntegerVector fused_assignment(IntegerMatrix nn_index,
                               NumericMatrix nn_dist,
                               IntegerVector curclus,
                               NumericMatrix coords,
                               NumericMatrix data_centroids,
                               NumericMatrix coord_centroids,
                               NumericMatrix data,
                               double dthresh,
                               double sigma1,
                               double sigma2,
                               double alpha) {
  int n = nn_index.nrow();
  IntegerVector out(n);
  int nswitches = 0;
  
  for(int i = 0; i < n; ++i) {
    // Get neighbors for this voxel
    NumericVector D = nn_dist.row(i);
    IntegerVector ind = nn_index.row(i);
    
    // Find unique cluster candidates in neighborhood
    std::unordered_set<int> cand_set;
    cand_set.insert(curclus[i]); // always include current cluster
    
    for (int j = 0; j < D.size(); ++j) {
      if (D[j] < dthresh) {
        int neighbor_label = curclus[ind[j]];
        cand_set.insert(neighbor_label);
      }
    }
    
    // If only one candidate, keep current assignment
    if (cand_set.size() == 1) {
      out[i] = curclus[i];
      continue;
    }
    
    // Find best cluster among candidates
    double best_score = -1.0;
    int best_cluster = curclus[i];
    
    NumericVector voxel_data = data(_, i);
    NumericVector voxel_coord = coords(_, i);
    
    for (int cid : cand_set) {
      int cluster_idx = cid - 1; // Convert to 0-based index
      
      // Compute data similarity
      NumericVector cluster_data = data_centroids(_, cluster_idx);
      double c1 = normalized_heat_kernel_opt(&voxel_data[0], &cluster_data[0], 
                                             voxel_data.size(), sigma1);
      
      // Compute spatial similarity  
      NumericVector cluster_coord = coord_centroids(_, cluster_idx);
      double c2 = heat_kernel_opt(&voxel_coord[0], &cluster_coord[0],
                                  voxel_coord.size(), sigma2);
      
      double score = alpha * c1 + (1.0 - alpha) * c2;
      
      if (score > best_score) {
        best_score = score;
        best_cluster = cid;
      }
    }
    
    out[i] = best_cluster;
    if (out[i] != curclus[i]) {
      nswitches++;
    }
  }
  
  out.attr("nswitches") = nswitches;
  return out;
}

// Parallel version of fused assignment
struct FusedAssignmentWorker : public Worker {
  // Inputs
  const RMatrix<int> nn_index;
  const RMatrix<double> nn_dist;
  const RVector<int> curclus;
  const RMatrix<double> coords;
  const RMatrix<double> data_centroids;
  const RMatrix<double> coord_centroids;
  const RMatrix<double> data;
  const double dthresh;
  const double sigma1;
  const double sigma2;
  const double alpha;
  
  // Output
  RVector<int> out;
  std::atomic<int>& nswitches;
  
  // Constructor
  FusedAssignmentWorker(const IntegerMatrix& nn_index_,
                       const NumericMatrix& nn_dist_,
                       const IntegerVector& curclus_,
                       const NumericMatrix& coords_,
                       const NumericMatrix& data_centroids_,
                       const NumericMatrix& coord_centroids_,
                       const NumericMatrix& data_,
                       double dthresh_,
                       double sigma1_,
                       double sigma2_,
                       double alpha_,
                       IntegerVector& out_,
                       std::atomic<int>& nswitches_)
    : nn_index(nn_index_),
      nn_dist(nn_dist_),
      curclus(curclus_),
      coords(coords_),
      data_centroids(data_centroids_),
      coord_centroids(coord_centroids_),
      data(data_),
      dthresh(dthresh_),
      sigma1(sigma1_),
      sigma2(sigma2_),
      alpha(alpha_),
      out(out_),
      nswitches(nswitches_) {}
  
  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      // Find unique cluster candidates in neighborhood
      std::unordered_set<int> cand_set;
      cand_set.insert(curclus[i]); // always include current cluster
      
      for (int j = 0; j < nn_index.ncol(); ++j) {
        if (nn_dist(i, j) < dthresh) {
          int neighbor_idx = nn_index(i, j);
          int neighbor_label = curclus[neighbor_idx];
          cand_set.insert(neighbor_label);
        }
      }
      
      // If only one candidate, keep current assignment
      if (cand_set.size() == 1) {
        out[i] = curclus[i];
        continue;
      }
      
      // Find best cluster among candidates
      double best_score = -1.0;
      int best_cluster = curclus[i];
      
      for (int cid : cand_set) {
        int cluster_idx = cid - 1; // Convert to 0-based index
        
        // Compute data similarity
        double c1 = normalized_heat_kernel_opt(
          data.column(i).begin(),
          data_centroids.column(cluster_idx).begin(),
          data.nrow(), 
          sigma1
        );
        
        // Compute spatial similarity
        double c2 = heat_kernel_opt(
          coords.column(i).begin(),
          coord_centroids.column(cluster_idx).begin(),
          coords.nrow(),
          sigma2
        );
        
        double score = alpha * c1 + (1.0 - alpha) * c2;
        
        if (score > best_score) {
          best_score = score;
          best_cluster = cid;
        }
      }
      
      out[i] = best_cluster;
      if (out[i] != curclus[i]) {
        nswitches.fetch_add(1);
      }
    }
  }
};

// [[Rcpp::export]]
IntegerVector fused_assignment_parallel(IntegerMatrix nn_index,
                                       NumericMatrix nn_dist,
                                       IntegerVector curclus,
                                       NumericMatrix coords,
                                       NumericMatrix data_centroids,
                                       NumericMatrix coord_centroids,
                                       NumericMatrix data,
                                       double dthresh,
                                       double sigma1,
                                       double sigma2,
                                       double alpha,
                                       int grain_size = 100) {
  int n = nn_index.nrow();
  IntegerVector out(n);
  std::atomic<int> nswitches(0);
  
  // Create and run parallel worker
  FusedAssignmentWorker worker(nn_index, nn_dist, curclus, coords,
                               data_centroids, coord_centroids, data,
                               dthresh, sigma1, sigma2, alpha,
                               out, nswitches);
  
  parallelFor(0, n, worker, grain_size);
  
  out.attr("nswitches") = nswitches.load();
  return out;
}

// Parallel centroid computation
struct CentroidWorker : public Worker {
  // Inputs
  const RVector<int> cluster_ids;
  const RMatrix<double> data;
  const RMatrix<double> coords;
  const int n_clusters;
  
  // Outputs (preallocated)
  RMatrix<double> data_centroids;
  RMatrix<double> coord_centroids;
  RVector<int> cluster_counts;
  
  // Constructor
  CentroidWorker(const IntegerVector& cluster_ids_,
                 const NumericMatrix& data_,
                 const NumericMatrix& coords_,
                 int n_clusters_,
                 NumericMatrix& data_centroids_,
                 NumericMatrix& coord_centroids_,
                 IntegerVector& cluster_counts_)
    : cluster_ids(cluster_ids_),
      data(data_),
      coords(coords_),
      n_clusters(n_clusters_),
      data_centroids(data_centroids_),
      coord_centroids(coord_centroids_),
      cluster_counts(cluster_counts_) {}
  
  // Initialize to zero
  void operator()(std::size_t begin, std::size_t end) {
    // Each thread handles a subset of clusters
    for (std::size_t k = begin; k < end; k++) {
      // Zero out this cluster's accumulators
      for (int j = 0; j < data.nrow(); j++) {
        data_centroids(j, k) = 0.0;
      }
      for (int j = 0; j < coords.nrow(); j++) {
        coord_centroids(j, k) = 0.0;
      }
      cluster_counts[k] = 0;
      
      // Accumulate for this cluster
      for (int i = 0; i < cluster_ids.size(); i++) {
        if (cluster_ids[i] == (k + 1)) { // 1-based indexing
          cluster_counts[k]++;
          for (int j = 0; j < data.nrow(); j++) {
            data_centroids(j, k) += data(j, i);
          }
          for (int j = 0; j < coords.nrow(); j++) {
            coord_centroids(j, k) += coords(j, i);
          }
        }
      }
      
      // Divide by count to get mean
      if (cluster_counts[k] > 0) {
        double inv_count = 1.0 / cluster_counts[k];
        for (int j = 0; j < data.nrow(); j++) {
          data_centroids(j, k) *= inv_count;
        }
        for (int j = 0; j < coords.nrow(); j++) {
          coord_centroids(j, k) *= inv_count;
        }
      }
    }
  }
};

// [[Rcpp::export]]
List compute_centroids_parallel(IntegerVector cluster_ids,
                                NumericMatrix data,
                                NumericMatrix coords,
                                int n_clusters,
                                int grain_size = 10) {
  // First, find the actual unique cluster IDs and create a mapping
  std::set<int> unique_ids;
  for (int i = 0; i < cluster_ids.size(); i++) {
    unique_ids.insert(cluster_ids[i]);
  }
  
  // Create mapping from cluster ID to index
  std::map<int, int> id_to_idx;
  int idx = 0;
  for (int id : unique_ids) {
    id_to_idx[id] = idx++;
  }
  
  int actual_n_clusters = unique_ids.size();
  
  // Preallocate output matrices based on actual number of clusters
  NumericMatrix data_centroids(data.nrow(), actual_n_clusters);
  NumericMatrix coord_centroids(coords.nrow(), actual_n_clusters);
  IntegerVector cluster_counts(actual_n_clusters);
  
  // Sequential computation (safer for now to avoid parallel issues with mapping)
  for (int k = 0; k < actual_n_clusters; k++) {
    // Zero out this cluster's accumulators
    for (int j = 0; j < data.nrow(); j++) {
      data_centroids(j, k) = 0.0;
    }
    for (int j = 0; j < coords.nrow(); j++) {
      coord_centroids(j, k) = 0.0;
    }
    cluster_counts[k] = 0;
  }
  
  // Accumulate for each cluster
  for (int i = 0; i < cluster_ids.size(); i++) {
    int cluster_id = cluster_ids[i];
    int cluster_idx = id_to_idx[cluster_id];
    
    cluster_counts[cluster_idx]++;
    for (int j = 0; j < data.nrow(); j++) {
      data_centroids(j, cluster_idx) += data(j, i);
    }
    for (int j = 0; j < coords.nrow(); j++) {
      coord_centroids(j, cluster_idx) += coords(j, i);
    }
  }
  
  // Divide by count to get mean
  for (int k = 0; k < actual_n_clusters; k++) {
    if (cluster_counts[k] > 0) {
      double inv_count = 1.0 / cluster_counts[k];
      for (int j = 0; j < data.nrow(); j++) {
        data_centroids(j, k) *= inv_count;
      }
      for (int j = 0; j < coords.nrow(); j++) {
        coord_centroids(j, k) *= inv_count;
      }
    }
  }
  
  return List::create(
    Named("centers") = data_centroids,
    Named("coord_centers") = coord_centroids,
    Named("counts") = cluster_counts
  );
}