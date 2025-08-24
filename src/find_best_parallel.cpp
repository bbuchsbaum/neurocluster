// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <unordered_set>
#include <cmath>
#include <atomic>

using namespace Rcpp;
using namespace RcppParallel;

// Sequential heat kernel functions (same as before)
inline double heat_kernel_impl(const double* x1, const double* x2, int len, double sigma) {
  double dist_sq = 0.0;
  for (int i = 0; i < len; i++) {
    double diff = x1[i] - x2[i];
    dist_sq += diff * diff;
  }
  // Fixed: Use dist_sq directly in Gaussian kernel, not sqrt(dist_sq)
  return std::exp(-dist_sq / (2.0 * sigma * sigma));
}

inline double normalized_heat_kernel_impl(const double* x1, const double* x2, int len, double sigma) {
  double dist_sq = 0.0;
  for (int i = 0; i < len; i++) {
    double diff = x1[i] - x2[i];
    dist_sq += diff * diff;
  }
  double norm_dist = dist_sq / (2.0 * len);
  return std::exp(-norm_dist / (2.0 * sigma * sigma));
}

// Parallel worker for best_candidate computation
struct BestCandidateWorker : public Worker {
  // Inputs
  const List& candidates;
  const RVector<int> curclus;
  const RMatrix<double> coords;
  const RMatrix<double> data_centroids;
  const RMatrix<double> coord_centroids;
  const RMatrix<double> data;
  const double sigma1;
  const double sigma2;
  const double alpha;
  
  // Output
  RVector<int> out;
  std::atomic<int>& nswitches;
  
  // Constructor
  BestCandidateWorker(const List& candidates_,
                      const IntegerVector& curclus_,
                      const NumericMatrix& coords_,
                      const NumericMatrix& data_centroids_,
                      const NumericMatrix& coord_centroids_,
                      const NumericMatrix& data_,
                      double sigma1_,
                      double sigma2_,
                      double alpha_,
                      IntegerVector& out_,
                      std::atomic<int>& nswitches_)
    : candidates(candidates_),
      curclus(curclus_),
      coords(coords_),
      data_centroids(data_centroids_),
      coord_centroids(coord_centroids_),
      data(data_),
      sigma1(sigma1_),
      sigma2(sigma2_),
      alpha(alpha_),
      out(out_),
      nswitches(nswitches_) {}
  
  // Parallel operator
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      IntegerVector cand = candidates[i];
      
      if (cand.size() <= 1) {
        out[i] = curclus[i];
        continue;
      }
      
      // Find best candidate for this voxel
      double best_score = -1.0;
      int best_cluster = curclus[i];
      
      for (int j = 0; j < cand.size(); j++) {
        int cid = cand[j] - 1;
        
        // Compute data similarity
        // Note: matrices are transposed when passed from R, so columns are rows here
        double c1 = normalized_heat_kernel_impl(
          data.column(i).begin(), 
          data_centroids.column(cid).begin(),
          data.nrow(), 
          sigma1
        );
        
        // Compute spatial similarity
        double c2 = heat_kernel_impl(
          coords.column(i).begin(),
          coord_centroids.column(cid).begin(),
          coords.nrow(),
          sigma2
        );
        
        double score = alpha * c1 + (1.0 - alpha) * c2;
        
        if (score > best_score) {
          best_score = score;
          best_cluster = cand[j];
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
IntegerVector best_candidate_parallel(List candidates,
                                      IntegerVector curclus,
                                      NumericMatrix coords,
                                      NumericMatrix data_centroids,
                                      NumericMatrix coord_centroids,
                                      NumericMatrix data,
                                      double sigma1,
                                      double sigma2,
                                      double alpha,
                                      int grain_size = 100) {
  int n = candidates.size();
  IntegerVector out(n);
  std::atomic<int> nswitches(0);
  
  // Create and run parallel worker
  BestCandidateWorker worker(candidates, curclus, coords, 
                             data_centroids, coord_centroids,
                             data, sigma1, sigma2, alpha,
                             out, nswitches);
  
  // Use grain size for load balancing
  // Smaller grain = better load balancing but more overhead
  parallelFor(0, n, worker, grain_size);
  
  out.attr("nswitches") = nswitches.load();
  return out;
}

// Keep the sequential version for comparison/fallback
// [[Rcpp::export]]
IntegerVector best_candidate_sequential(List candidates,
                                        IntegerVector curclus,
                                        NumericMatrix coords,
                                        NumericMatrix data_centroids,
                                        NumericMatrix coord_centroids,
                                        NumericMatrix data,
                                        double sigma1,
                                        double sigma2,
                                        double alpha) {
  int n = candidates.size();
  IntegerVector out(n);
  int nswitches = 0;
  
  for (int i = 0; i < n; i++) {
    IntegerVector cand = candidates[i];
    if (cand.size() <= 1) {
      out[i] = curclus[i];
      continue;
    }
    
    NumericVector score(cand.size());
    for (int j = 0; j < cand.size(); j++) {
      int cid = cand[j] - 1;
      
      // Get column vectors
      NumericVector data_i = data(_, i);
      NumericVector data_cent = data_centroids(_, cid);
      double c1 = normalized_heat_kernel_impl(
        &data_i[0], &data_cent[0], data_i.size(), sigma1
      );
      
      NumericVector coord_i = coords(_, i);
      NumericVector coord_cent = coord_centroids(_, cid);
      double c2 = heat_kernel_impl(
        &coord_i[0], &coord_cent[0], coord_i.size(), sigma2
      );
      
      score[j] = alpha * c1 + (1.0 - alpha) * c2;
    }
    
    int best_j = which_max(score);
    out[i] = cand[best_j];
    if (out[i] != curclus[i]) {
      nswitches++;
    }
  }
  
  out.attr("nswitches") = nswitches;
  return out;
}

// Export the original functions to maintain compatibility
// [[Rcpp::export]]
double heat_kernel(NumericVector x1, NumericVector x2, double sigma) {
  return heat_kernel_impl(&x1[0], &x2[0], x1.size(), sigma);
}

// [[Rcpp::export]]
double normalized_heat_kernel(NumericVector x1, NumericVector x2, double sigma) {
  return normalized_heat_kernel_impl(&x1[0], &x2[0], x1.size(), sigma);
}

// No need to re-export these functions - they're already in find_best.cpp
// We'll just use the parallel versions with different names