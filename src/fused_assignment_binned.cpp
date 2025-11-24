// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <stdint.h>

using namespace Rcpp;
using namespace RcppParallel;

// Simple hash for 3D grid index (if we needed unordered_map).
// In this implementation we use a dense 3D grid flattened to 1D, so we don't
// need a hash map; keeping this around for reference.
// static inline uint64_t pack3i(int x, int y, int z) {
//   return ( (uint64_t)(uint32_t)x << 42 ) ^ ( (uint64_t)(uint32_t)y << 21 ) ^ (uint64_t)(uint32_t)z;
// }

struct GridIndex {
  // Grid step size in world units (voxel space)
  double step;
  // Bounding box minima
  double minx, miny, minz;
  // Grid dimensions
  int nx, ny, nz;
  // Flattened bins: each cell holds a small vector of centroid IDs
  std::vector< std::vector<int> > bins;

  GridIndex() : step(1.0), minx(0), miny(0), minz(0), nx(0), ny(0), nz(0) {}

  inline int clampi(int v, int lo, int hi) const {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
  }

  inline int index_of(int ix, int iy, int iz) const {
    return (iz * ny + iy) * nx + ix;
  }

  inline void coord_to_cell(double x, double y, double z, int &ix, int &iy, int &iz) const {
    ix = (int)std::floor((x - minx) / step);
    iy = (int)std::floor((y - miny) / step);
    iz = (int)std::floor((z - minz) / step);
    ix = clampi(ix, 0, nx - 1);
    iy = clampi(iy, 0, ny - 1);
    iz = clampi(iz, 0, nz - 1);
  }
};

// Build a regular 3D grid over the voxel coordinate bounding box and bin centroids
// into cells. The grid step is chosen so there are ~K cells (heuristic), which
// keeps ~1-2 centroids per cell on average.
static GridIndex build_grid(const NumericMatrix &coords,          // 3 x n
                            const NumericMatrix &coord_centroids  // 3 x K
                            ) {
  GridIndex g;
  const int n = coords.ncol();
  const int K = coord_centroids.ncol();
  if (n == 0 || K == 0) return g;

  // Bounding box of voxel coords
  double minx =  std::numeric_limits<double>::infinity();
  double miny =  std::numeric_limits<double>::infinity();
  double minz =  std::numeric_limits<double>::infinity();
  double maxx = -std::numeric_limits<double>::infinity();
  double maxy = -std::numeric_limits<double>::infinity();
  double maxz = -std::numeric_limits<double>::infinity();

  // 3 rows: x,y,z
  for (int i = 0; i < n; ++i) {
    double x = coords(0, i);
    double y = coords(1, i);
    double z = coords(2, i);
    if (x < minx) minx = x; if (x > maxx) maxx = x;
    if (y < miny) miny = y; if (y > maxy) maxy = y;
    if (z < minz) minz = z; if (z > maxz) maxz = z;
  }

  // Guard against degenerate ranges
  const double eps = 1e-9;
  if (maxx - minx < eps) { maxx = minx + 1.0; }
  if (maxy - miny < eps) { maxy = miny + 1.0; }
  if (maxz - minz < eps) { maxz = minz + 1.0; }

  const double rangex = maxx - minx;
  const double rangey = maxy - miny;
  const double rangez = maxz - minz;
  const double volume = rangex * rangey * rangez;

  // Approximate centroid spacing S ≈ (volume / K)^(1/3)
  // Use a small safety factor (0.9) so neighborhood search (±1 cells)
  // roughly corresponds to a radius of ~2S in world units.
  const double S = std::cbrt(std::max(volume / (double)K, 1e-9)) * 0.9;

  // Choose grid step equal to S and compute integer grid sizes
  g.step = S;
  g.minx = minx; g.miny = miny; g.minz = minz;
  g.nx = std::max(1, (int)std::floor(rangex / S) + 1);
  g.ny = std::max(1, (int)std::floor(rangey / S) + 1);
  g.nz = std::max(1, (int)std::floor(rangez / S) + 1);

  // Avoid pathological memory if ranges are tiny
  const int64_t total_cells = (int64_t)g.nx * (int64_t)g.ny * (int64_t)g.nz;
  if (total_cells > 10000000LL) {
    // Fallback: coarsen grid
    double factor = std::cbrt((double)total_cells / 10000000.0);
    g.nx = std::max(1, (int)std::floor(g.nx / factor));
    g.ny = std::max(1, (int)std::floor(g.ny / factor));
    g.nz = std::max(1, (int)std::floor(g.nz / factor));
  }

  g.bins.clear();
  g.bins.resize((size_t)g.nx * g.ny * g.nz);

  // Bin centroids
  for (int k = 0; k < K; ++k) {
    double cx = coord_centroids(0, k);
    double cy = coord_centroids(1, k);
    double cz = coord_centroids(2, k);
    int ix, iy, iz;
    g.coord_to_cell(cx, cy, cz, ix, iy, iz);
    const int idx = g.index_of(ix, iy, iz);
    g.bins[(size_t)idx].push_back(k);
  }

  return g;
}

struct BinnedAssignWorker : public Worker {
  // Inputs (views)
  const RMatrix<int> nn_index;    // n x knn
  const RMatrix<double> nn_dist;  // n x knn
  const RVector<int> curclus;     // n
  const RMatrix<double> coords;   // 3 x n
  const RMatrix<double> data_centroids;   // D x K
  const RMatrix<double> coord_centroids;  // 3 x K
  const RMatrix<double> data;     // D x n
  const double dthresh;
  const double inv2sigma1;
  const double inv2sigma2;
  const double alpha;
  const double one_minus_alpha;
  const int n_clusters;

  // Spatial binning (shared, read-only)
  const GridIndex &grid;
  const int bin_expand;       // how many neighbor cells to visit (±bin_expand)
  const double window_radius2; // squared world radius to accept a centroid candidate

  // Output
  RVector<int> out;

  // Thread-local scratch
  // We keep these mutable and re-used across voxels to avoid allocations.
  mutable std::vector<int> cand;
  mutable std::vector<int> seen_stamp;
  mutable int cur_stamp;

  BinnedAssignWorker(const IntegerMatrix &nn_index_,
                     const NumericMatrix &nn_dist_,
                     const IntegerVector &curclus_,
                     const NumericMatrix &coords_,
                     const NumericMatrix &data_centroids_,
                     const NumericMatrix &coord_centroids_,
                     const NumericMatrix &data_,
                     double dthresh_,
                     double sigma1_,
                     double sigma2_,
                     double alpha_,
                     const GridIndex &grid_,
                     int bin_expand_,
                     double window_radius_,
                     IntegerVector &out_)
    : nn_index(nn_index_), nn_dist(nn_dist_), curclus(curclus_),
      coords(coords_), data_centroids(data_centroids_),
      coord_centroids(coord_centroids_), data(data_),
      dthresh(dthresh_),
      inv2sigma1(1.0 / (2.0 * sigma1_ * sigma1_)),
      inv2sigma2(1.0 / (2.0 * sigma2_ * sigma2_)),
      alpha(alpha_),
      one_minus_alpha(1.0 - alpha_),
      n_clusters(coord_centroids_.ncol()),
      grid(grid_),
      bin_expand(bin_expand_),
      window_radius2(window_radius_ * window_radius_),
      out(out_),
      cand(),
      seen_stamp(n_clusters, 0),
      cur_stamp(1)
  {
    cand.reserve(64);
  }

  inline void reset_seen_if_needed() const {
    // Prevent overflow of cur_stamp; reset stamps if necessary.
    if (cur_stamp == std::numeric_limits<int>::max()) {
      std::fill(seen_stamp.begin(), seen_stamp.end(), 0);
      cur_stamp = 1;
    } else {
      ++cur_stamp;
    }
  }

  // Gather candidate clusters from spatial bins around the voxel coordinate.
  inline void gather_bin_candidates(double vx, double vy, double vz) const {
    cand.clear();
    int ix, iy, iz;
    grid.coord_to_cell(vx, vy, vz, ix, iy, iz);

    const int x0 = std::max(0, ix - bin_expand);
    const int x1 = std::min(grid.nx - 1, ix + bin_expand);
    const int y0 = std::max(0, iy - bin_expand);
    const int y1 = std::min(grid.ny - 1, iy + bin_expand);
    const int z0 = std::max(0, iz - bin_expand);
    const int z1 = std::min(grid.nz - 1, iz + bin_expand);

    for (int cz = z0; cz <= z1; ++cz) {
      for (int cy = y0; cy <= y1; ++cy) {
        for (int cx = x0; cx <= x1; ++cx) {
          const int cell_idx = grid.index_of(cx, cy, cz);
          const std::vector<int> &vec = grid.bins[(size_t)cell_idx];
          for (size_t t = 0; t < vec.size(); ++t) {
            const int cid = vec[t];
            if (seen_stamp[cid] == cur_stamp) continue; // already recorded for this voxel
            // Optional exact radius filter (avoid far centroids even if in neighbor cells)
            double dx = coord_centroids(0, cid) - vx;
            double dy = coord_centroids(1, cid) - vy;
            double dz = coord_centroids(2, cid) - vz;
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 > window_radius2) continue;
            seen_stamp[cid] = cur_stamp;
            cand.push_back(cid);
          }
        }
      }
    }
  }

  // Optionally intersect with neighbor-assigned clusters to preserve locality
  inline void intersect_with_neighbor_clusters(int voxel_index) const {
    const int knn = nn_index.ncol();
    // Only prune when there are materially more clusters than immediate neighbors.
    if (n_clusters <= knn) return;
    if (cand.empty()) return; // nothing to intersect; fall back later if needed
    // Mark which candidates are kept (we'll rebuild cand in-place)
    int write_pos = 0;
    // Mark neighbor clusters from nn graph (within dthresh)
    // We'll use a small stamp array keyed by cluster id to avoid a set.
    // To avoid an extra array, we reuse seen_stamp with a new stamp value.
    int nb_stamp_val = cur_stamp + 1;
    if (nb_stamp_val == std::numeric_limits<int>::max()) nb_stamp_val = 1; // extremely unlikely
    // Temporarily mark nothing
    // Now mark neighbor clusters
    for (int j = 0; j < knn; ++j) {
      int nb = nn_index(voxel_index, j);
      if (nb < 0) continue;
      double d = nn_dist(voxel_index, j);
      if (dthresh > 0 && d > dthresh) continue;
      int nb_cid = curclus[nb];
      if (nb_cid >= 0 && nb_cid < n_clusters) {
        seen_stamp[nb_cid] = nb_stamp_val;
      }
    }
    // Keep only those cand IDs that were marked by neighbor clusters
    for (size_t u = 0; u < cand.size(); ++u) {
      int cid = cand[u];
      if (seen_stamp[cid] == nb_stamp_val) {
        cand[write_pos++] = cid;
      }
    }
    cand.resize(write_pos);
  }

  inline int choose_best_for_voxel(int i) const {
    // If cand empty, optionally fall back to neighbor clusters
    if (cand.empty()) {
      // Gather neighbor cluster IDs as candidates
      const int knn = nn_index.ncol();
      for (int j = 0; j < knn; ++j) {
        int nb = nn_index(i, j);
        if (nb < 0) continue;
        double d = nn_dist(i, j);
        if (dthresh > 0 && d > dthresh) continue;
        int cid = curclus[nb];
        if (cid >= 0 && cid < n_clusters) {
          if (seen_stamp[cid] != cur_stamp) {
            seen_stamp[cid] = cur_stamp;
            cand.push_back(cid);
          }
        }
      }
      // Still empty? As a last resort, consider current assignment only.
      if (cand.empty()) {
        int cur = curclus[i];
        if (cur >= 0 && cur < n_clusters) cand.push_back(cur);
      }
    }

    const int D = data.nrow();
    const double *x = &data(0, i);
    const double vx = coords(0, i);
    const double vy = coords(1, i);
    const double vz = coords(2, i);

    double best_score = -std::numeric_limits<double>::infinity();
    int best_cluster = curclus[i];

    for (size_t t = 0; t < cand.size(); ++t) {
      const int cid = cand[t];
      const double *c = &data_centroids(0, cid);

      // data distance squared
      double dd2 = 0.0;
      for (int j = 0; j < D; ++j) {
        double diff = x[j] - c[j];
        dd2 += diff * diff;
      }

      // spatial distance squared
      double dx = vx - coord_centroids(0, cid);
      double dy = vy - coord_centroids(1, cid);
      double dz = vz - coord_centroids(2, cid);
      double sd2 = dx*dx + dy*dy + dz*dz;

      // combined score (exact)
      double score = alpha * std::exp(-dd2 * inv2sigma1) +
                     one_minus_alpha * std::exp(-sd2 * inv2sigma2);
      if (score > best_score) {
        best_score = score;
        best_cluster = cid;
      }
    }
    return best_cluster;
  }

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      reset_seen_if_needed();

      const double vx = coords(0, i);
      const double vy = coords(1, i);
      const double vz = coords(2, i);

      // 1) Spatial bin candidates around (vx,vy,vz)
      gather_bin_candidates(vx, vy, vz);

      // 2) Intersect with neighbor clusters (within dthresh) to preserve locality.
      //    If intersection becomes empty, we'll fallback to neighbors inside choose_best_for_voxel.
      intersect_with_neighbor_clusters((int)i);

      // 3) Score only the remaining few candidates
      const int best_c = choose_best_for_voxel((int)i);
      out[i] = best_c;
    }
  }
};

// [[Rcpp::export]]
IntegerVector fused_assignment_parallel_binned(IntegerMatrix nn_index,
                                               NumericMatrix nn_dist,
                                               IntegerVector curclus,
                                               NumericMatrix coords,           // 3 x n (columns are voxels)
                                               NumericMatrix data_centroids,   // D x K
                                               NumericMatrix coord_centroids,  // 3 x K
                                               NumericMatrix data,             // D x n
                                               double dthresh,
                                               double sigma1,
                                               double sigma2,
                                               double alpha,
                                               int grain_size = 2048,
                                               double window_factor = 2.0,
                                               int bin_expand = 1) {
  const int n = coords.ncol();
  if (n == 0) return IntegerVector();

  // Basic validation
  if (data.ncol() != n || nn_index.nrow() != n || nn_dist.nrow() != n) {
    stop("Dimension mismatch: ncols(data|coords) and nrows(nn_index|nn_dist) must match n voxels.");
  }
  if (data_centroids.nrow() != data.nrow()) stop("data_centroids rows must match data rows (feature dimension).");
  if (coord_centroids.nrow() != coords.nrow()) stop("coord_centroids rows must match coords rows (usually 3).");
  if (alpha < 0.0 || alpha > 1.0) stop("alpha must be in [0,1].");
  if (sigma1 <= 0.0 || sigma2 <= 0.0) stop("sigmas must be positive.");
  if (bin_expand < 0) stop("bin_expand must be >= 0.");
  if (window_factor <= 0.0) window_factor = 2.0;

  // Build centroid grid once per iteration
  GridIndex grid = build_grid(coords, coord_centroids);
  // Spatial window radius ~ window_factor * S
  const double window_radius = window_factor * grid.step;

  IntegerVector out(n);

  BinnedAssignWorker worker(nn_index, nn_dist, curclus, coords,
                            data_centroids, coord_centroids, data,
                            dthresh, sigma1, sigma2, alpha,
                            grid, bin_expand, window_radius, out);

  if (grain_size <= 0) {
    // Heuristic: aim for ~ (#threads * 8) chunks
    int threads = 8; // RcppParallel picks actual threads; this is just to size the chunks
    grain_size = std::max(1024, n / (threads * 8));
  }

  parallelFor(0, (size_t)n, worker, (size_t)grain_size);
  return out;
}

// A sequential version using the same binning; helpful if you want to keep
// parallel=FALSE but still benefit from a bounded candidate set.
// [[Rcpp::export]]
IntegerVector fused_assignment_binned(IntegerMatrix nn_index,
                                      NumericMatrix nn_dist,
                                      IntegerVector curclus,
                                      NumericMatrix coords,           // 3 x n (columns are voxels)
                                      NumericMatrix data_centroids,   // D x K
                                      NumericMatrix coord_centroids,  // 3 x K
                                      NumericMatrix data,             // D x n
                                      double dthresh,
                                      double sigma1,
                                      double sigma2,
                                      double alpha,
                                      double window_factor = 2.0,
                                      int bin_expand = 1) {
  const int n = coords.ncol();
  if (n == 0) return IntegerVector();

  GridIndex grid = build_grid(coords, coord_centroids);
  const double window_radius = window_factor * grid.step;

  const int K = coord_centroids.ncol();
  const int D = data.nrow();

  const double inv2sigma1 = 1.0 / (2.0 * sigma1 * sigma1);
  const double inv2sigma2 = 1.0 / (2.0 * sigma2 * sigma2);
  const double one_minus_alpha = 1.0 - alpha;

  IntegerVector out(n);

  std::vector<int> seen_stamp((size_t)K, 0);
  int cur_stamp = 1;
  std::vector<int> cand; cand.reserve(64);

  auto reset_seen_if_needed = [&](){
    if (cur_stamp == std::numeric_limits<int>::max()) {
      std::fill(seen_stamp.begin(), seen_stamp.end(), 0);
      cur_stamp = 1;
    } else {
      ++cur_stamp;
    }
  };

  auto gather_bin_candidates = [&](double vx, double vy, double vz){
    cand.clear();
    int ix, iy, iz;
    grid.coord_to_cell(vx, vy, vz, ix, iy, iz);

    const int x0 = std::max(0, ix - bin_expand);
    const int x1 = std::min(grid.nx - 1, ix + bin_expand);
    const int y0 = std::max(0, iy - bin_expand);
    const int y1 = std::min(grid.ny - 1, iy + bin_expand);
    const int z0 = std::max(0, iz - bin_expand);
    const int z1 = std::min(grid.nz - 1, iz + bin_expand);

    const double wr2 = window_radius * window_radius;

    for (int cz = z0; cz <= z1; ++cz) {
      for (int cy = y0; cy <= y1; ++cy) {
        for (int cx = x0; cx <= x1; ++cx) {
          const int cell_idx = grid.index_of(cx, cy, cz);
          const std::vector<int> &vec = grid.bins[(size_t)cell_idx];
          for (size_t t = 0; t < vec.size(); ++t) {
            const int cid = vec[t];
            if (seen_stamp[cid] == cur_stamp) continue;
            double dx = coord_centroids(0, cid) - vx;
            double dy = coord_centroids(1, cid) - vy;
            double dz = coord_centroids(2, cid) - vz;
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 > wr2) continue;
            seen_stamp[cid] = cur_stamp;
            cand.push_back(cid);
          }
        }
      }
    }
  };

  auto intersect_with_neighbor_clusters = [&](int i){
    const int knn = nn_index.ncol();
    // Skip pruning when cluster count is small relative to neighbor count.
    if (K <= knn) return;
    if (cand.empty()) return;
    int nb_stamp_val = cur_stamp + 1;
    if (nb_stamp_val == std::numeric_limits<int>::max()) nb_stamp_val = 1;
    for (int j = 0; j < knn; ++j) {
      int nb = nn_index(i, j);
      if (nb < 0) continue;
      double d = nn_dist(i, j);
      if (dthresh > 0 && d > dthresh) continue;
      int nb_cid = curclus[nb];
      if (nb_cid >= 0 && nb_cid < K) seen_stamp[nb_cid] = nb_stamp_val;
    }
    int w = 0;
    for (size_t u = 0; u < cand.size(); ++u) {
      int cid = cand[u];
      if (seen_stamp[cid] == nb_stamp_val) cand[w++] = cid;
    }
    cand.resize(w);
  };

  for (int i = 0; i < n; ++i) {
    reset_seen_if_needed();

    const double vx = coords(0, i);
    const double vy = coords(1, i);
    const double vz = coords(2, i);

    gather_bin_candidates(vx, vy, vz);
    intersect_with_neighbor_clusters(i);

    if (cand.empty()) {
      // fallback to neighbor clusters
      const int knn = nn_index.ncol();
      for (int j = 0; j < knn; ++j) {
        int nb = nn_index(i, j);
        if (nb < 0) continue;
        double d = nn_dist(i, j);
        if (dthresh > 0 && d > dthresh) continue;
        int cid = curclus[nb];
        if (cid >= 0 && cid < K) {
          if (seen_stamp[cid] != cur_stamp) {
            seen_stamp[cid] = cur_stamp;
            cand.push_back(cid);
          }
        }
      }
      if (cand.empty()) {
        int cur = curclus[i];
        if (cur >= 0 && cur < K) cand.push_back(cur);
      }
    }

    const int Dloc = D;
    const double *x = &data(0, i);

    double best_score = -std::numeric_limits<double>::infinity();
    int best_cluster = curclus[i];

    for (size_t t = 0; t < cand.size(); ++t) {
      const int cid = cand[t];
      const double *c = &data_centroids(0, cid);

      double dd2 = 0.0;
      for (int j = 0; j < Dloc; ++j) {
        double diff = x[j] - c[j];
        dd2 += diff * diff;
      }

      double dx = vx - coord_centroids(0, cid);
      double dy = vy - coord_centroids(1, cid);
      double dz = vz - coord_centroids(2, cid);
      double sd2 = dx*dx + dy*dy + dz*dz;

      double score = alpha * std::exp(-dd2 * inv2sigma1) +
                     one_minus_alpha * std::exp(-sd2 * inv2sigma2);
      if (score > best_score) {
        best_score = score;
        best_cluster = cid;
      }
    }
    out[i] = best_cluster;
  }

  return out;
}
