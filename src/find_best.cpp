// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <cmath> // std::sqrt, std::exp

using namespace Rcpp;

/*
 These C++ functions support the supervoxel clustering approach:
 - heat_kernel: Gaussian kernel on Euclidean distance
 - normalized_heat_kernel: normalized by dimension
 - compute_scores: compute combined heat-kernel scores for each voxel
 - best_candidate: picks the best cluster among a set of candidate cluster labels
 - find_candidates: determines which cluster labels to consider
 based on local neighbors where distance < dthresh
 */

// Internal implementation - not exported to avoid duplicate symbols
double heat_kernel_internal(NumericVector x1, NumericVector x2, double sigma) {
  double dist_sq = 0.0;
  int len = x1.size();
  for (int i = 0; i < len; i++) {
    double diff = x1[i] - x2[i];
    dist_sq += diff * diff;
  }
  double dist = std::sqrt(dist_sq);
  return std::exp(-dist / (2.0 * sigma * sigma));
}

// Internal implementation - not exported to avoid duplicate symbols
double normalized_heat_kernel_internal(NumericVector x1, NumericVector x2, double sigma) {
  double dist_sq = 0.0;
  int len = x1.size();
  for (int i = 0; i < len; i++) {
    double diff = x1[i] - x2[i];
    dist_sq += diff * diff;
  }
  double norm_dist = dist_sq / (2.0 * len);
  return std::exp(-norm_dist / (2.0 * sigma * sigma));
}

// [[Rcpp::export]]
NumericVector compute_scores(IntegerVector curclus,
                             NumericMatrix coords,
                             NumericMatrix data_centroids,
                             NumericMatrix coord_centroids,
                             NumericMatrix data,
                             double sigma1,
                             double sigma2) {
  int n = curclus.size();
  NumericVector out(n);

  for (int i = 0; i < n; i++) {
    int cluster_idx = curclus[i] - 1;
    NumericVector voxelData = data(_, i);
    NumericVector clusterData = data_centroids(_, cluster_idx);
    double c1 = normalized_heat_kernel_internal(voxelData, clusterData, sigma1);

    NumericVector voxelCoord = coords(_, i);
    NumericVector clusterCoord = coord_centroids(_, cluster_idx);
    double c2 = heat_kernel_internal(voxelCoord, clusterCoord, sigma2);

    out[i] = c1 + c2;
  }

  return out;
}

// [[Rcpp::export]]
IntegerVector best_candidate(List candidates,
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
      double c1 = normalized_heat_kernel_internal(data(_, i), data_centroids(_, cid), sigma1);
      double c2 = heat_kernel_internal(coords(_, i), coord_centroids(_, cid), sigma2);
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


// [[Rcpp::export]]
List find_candidates(IntegerMatrix nn_index,
                     NumericMatrix nn_dist,
                     IntegerVector curclus,
                     double dthresh) {
  int n = nn_index.nrow();
  List out(n);

  for(int i = 0; i < n; ++i) {

    // CHANGED:
    // Instead of ConstRow, just read the row as a vector:
    NumericVector D = nn_dist.row(i);
    IntegerVector ind = nn_index.row(i);

    std::unordered_set<int> cand;
    cand.insert(curclus[i]); // always include the voxel's own cluster

    for (int j = 0; j < D.size(); ++j) {
      if (D[j] < dthresh) {
        int neighbor_label = curclus[ind[j]];
        cand.insert(neighbor_label);
      }
    }

    // CHANGED:
    // Rcpp IntegerVector doesn't support .reserve().
    // Use a temporary std::vector<int> and then wrap it:
    std::vector<int> tmp;
    tmp.reserve(cand.size());
    for (auto c : cand) {
      tmp.push_back(c);
    }
    IntegerVector outCand = wrap(tmp);

    out[i] = outCand;
  }

  return out;
}
