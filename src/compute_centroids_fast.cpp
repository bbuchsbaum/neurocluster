// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace RcppParallel;

struct CentroidReducer : public Worker {
  const RMatrix<double> data;   // D x n
  const RMatrix<double> coords; // 3 x n
  const RVector<int> cluster_ids; // n
  const int K;
  const int D;

  // partial sums
  std::vector<double> data_sums;   // D*K
  std::vector<double> coord_sums;  // 3*K
  std::vector<int> counts;         // K

  CentroidReducer(const NumericMatrix &data_,
                  const NumericMatrix &coords_,
                  const IntegerVector &cluster_ids_,
                  int K_)
    : data(data_), coords(coords_), cluster_ids(cluster_ids_), K(K_), D(data_.nrow()),
      data_sums((size_t)data_.nrow() * K, 0.0),
      coord_sums(3u*K, 0.0),
      counts(K, 0) {}

  // helper because MSVC doesn't like capture of this->
  inline int D_() const { return D; }

  // Split constructor
  CentroidReducer(const CentroidReducer &other, Split)
    : data(other.data), coords(other.coords), cluster_ids(other.cluster_ids),
      K(other.K), D(other.D),
      data_sums((size_t)other.D*other.K, 0.0),
      coord_sums(3u*other.K, 0.0),
      counts(other.K, 0) {}

  void operator()(std::size_t begin, std::size_t end) {
    const int Dloc = D;
    for (std::size_t i = begin; i < end; ++i) {
      int cid = cluster_ids[i];
      if (cid < 0 || cid >= K) continue;
      counts[cid] += 1;
      // sum features
      const double *x = &data(0, (int)i);
      double *dst_data = &data_sums[(size_t)cid * Dloc];
      for (int j = 0; j < Dloc; ++j) {
        dst_data[j] += x[j];
      }
      // sum coords
      const double vx = coords(0, (int)i);
      const double vy = coords(1, (int)i);
      const double vz = coords(2, (int)i);
      double *dst_coord = &coord_sums[(size_t)cid * 3u];
      dst_coord[0] += vx;
      dst_coord[1] += vy;
      dst_coord[2] += vz;
    }
  }

  void join(const CentroidReducer &rhs) {
    // sum partials
    for (size_t t = 0; t < data_sums.size(); ++t) data_sums[t] += rhs.data_sums[t];
    for (size_t t = 0; t < coord_sums.size(); ++t) coord_sums[t] += rhs.coord_sums[t];
    for (int k = 0; k < K; ++k) counts[k] += rhs.counts[k];
  }
};

// [[Rcpp::export]]
List compute_centroids_parallel_fast(IntegerVector cluster_ids,
                                     NumericMatrix data,
                                     NumericMatrix coords,
                                     int n_clusters) {
  const int K = n_clusters;
  if (K <= 0) stop("n_clusters must be > 0");
  if (coords.nrow() != 3) stop("coords must be 3 x n");
  if (data.ncol() != coords.ncol()) stop("data and coords must have the same ncol");

  CentroidReducer red(data, coords, cluster_ids, K);
  parallelReduce(0, (size_t)coords.ncol(), red);

  // Build outputs
  NumericMatrix data_centroids(data.nrow(), K);
  NumericMatrix coord_centroids(3, K);
  IntegerVector counts(K);
  for (int k = 0; k < K; ++k) {
    int cnt = red.counts[k];
    counts[k] = cnt;
    if (cnt > 0) {
      double inv = 1.0 / (double)cnt;
      // features
      const double *src_data = &red.data_sums[(size_t)k * data.nrow()];
      for (int j = 0; j < data.nrow(); ++j) {
        data_centroids(j, k) = src_data[j] * inv;
      }
      // coords
      const double *src_coord = &red.coord_sums[(size_t)k * 3u];
      coord_centroids(0, k) = src_coord[0] * inv;
      coord_centroids(1, k) = src_coord[1] * inv;
      coord_centroids(2, k) = src_coord[2] * inv;
    } else {
      // leave zeros; caller may reinitialize empty clusters if desired
    }
  }

  return List::create(
    Named("centers") = data_centroids,
    Named("coord_centers") = coord_centroids,
    Named("counts") = counts
  );
}
