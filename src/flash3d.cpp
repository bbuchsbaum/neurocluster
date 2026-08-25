// FLASH-3D: Fast Low-rank Approximate Superclusters for Hemodynamics (3D)
// Native 3D, fMRI-aware superclustering via DCT-rank temporal hashes + 3D jump-flood propagation.
// Author: (c) 2025
//
// Build: Rcpp + RcppParallel (OpenMP/TBB underneath).
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdint.h>
#include <string.h>
#include <functional>

using namespace Rcpp;
using namespace RcppParallel;

// ----- portable popcount64 --------------------------------------------------
static inline uint32_t popcount64(uint64_t x) {
#if defined(_MSC_VER)
  return (uint32_t)__popcnt64(x);
#elif defined(__GNUC__) || defined(__clang__)
  return (uint32_t)__builtin_popcountll((unsigned long long)x);
#else
  // Optimized SWAR (SIMD Within A Register) algorithm for portable performance
  x -= (x >> 1) & 0x5555555555555555ULL;
  x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
  x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
  return (uint32_t)((x * 0x0101010101010101ULL) >> 56);
#endif
}

// ----- fixed comparator table (deterministic) -------------------------------
static void make_comparators(int M, int B, std::vector<uint8_t> &ci, std::vector<uint8_t> &cj) {
  ci.resize(B); cj.resize(B);
  // A one-frame structural image has only the DC coefficient. Its temporal
  // hash is therefore identically zero; avoid the impossible i != j loop.
  if (M <= 1) {
    std::fill(ci.begin(), ci.end(), 0);
    std::fill(cj.begin(), cj.end(), 0);
    return;
  }
  // simple LCG seeded constant for determinism
  uint64_t s = 1469598103934665603ull ^ (uint64_t)M << 32 ^ (uint64_t)B;
  auto rnd = [&](){
    s ^= s << 13; s ^= s >> 7; s ^= s << 17;
    return s;
  };
  for (int b = 0; b < B; ++b) {
    int i, j;
    do {
      i = (int)(rnd() % (uint64_t)M);
      j = (int)(rnd() % (uint64_t)M);
    } while (j == i);
    ci[b] = (uint8_t)i;
    cj[b] = (uint8_t)j;
  }
}

// ----- precompute DCT cos table (T x M) -------------------------------------
static void build_dct_table(int T, int M, std::vector<double> &ctab) {
  ctab.resize((size_t)T * (size_t)M);
  const double pi = 3.14159265358979323846;
  for (int t = 0; t < T; ++t) {
    double x = (t + 0.5) * (pi / (double)T);
    for (int m = 0; m < M; ++m) {
      ctab[(size_t)t * M + m] = std::cos(x * (double)m);
    }
  }
}

// ----- hashing worker --------------------------------------------------------
namespace {

struct HashWorker : public Worker {
  const RMatrix<double> X; // T x Nmask
  const int T;
  const int M;  // # low DCT coeffs
  const int B;  // # bits (64 or 128)
  const std::vector<double> &ctab;
  const std::vector<uint8_t> &ci;
  const std::vector<uint8_t> &cj;
  RVector<double> mean_out; // optional tails
  RVector<double> var_out;
  std::vector<uint64_t> &hashA; // output hashes (64 or lower half of 128)
  std::vector<uint64_t> &hashB; // upper half when B==128, else unused

  HashWorker(const NumericMatrix &X_,
             int M_,
             int B_,
             const std::vector<double> &ctab_,
             const std::vector<uint8_t> &ci_,
             const std::vector<uint8_t> &cj_,
             NumericVector &mean_,
             NumericVector &var_,
             std::vector<uint64_t> &hashA_,
             std::vector<uint64_t> &hashB_)
    : X(X_), T(X_.nrow()), M(M_), B(B_), ctab(ctab_), ci(ci_), cj(cj_),
      mean_out(mean_), var_out(var_), hashA(hashA_), hashB(hashB_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const int N = X.ncol();
    const int Mloc = M;
    const int Bloc = B;
    for (std::size_t idx = begin; idx < end; ++idx) {
      if ((int)idx >= N) break;
      // zscore
      const double *x = &X(0, (int)idx);
      double mu = 0.0;
      for (int t = 0; t < T; ++t) mu += x[t];
      mu /= (double)T;
      double v = 0.0;
      for (int t = 0; t < T; ++t) { double d = x[t] - mu; v += d * d; }
      v = (v <= 0.0) ? 1.0 : std::sqrt(v / (double)T);
      mean_out[idx] = mu;
      var_out[idx]  = v * v;

      // DCT-M
      double c[32]; // support up to M<=32
      for (int m = 0; m < Mloc; ++m) {
        double acc = 0.0;
        for (int t = 0; t < T; ++t) {
          // (t*M + m)
          acc += ((x[t] - mu) / v) * ctab[(size_t)t * Mloc + m];
        }
        c[m] = acc;
      }

      // rank hash (comparators)
      uint64_t ha = 0ull, hb = 0ull;
      for (int b = 0; b < Bloc; ++b) {
        int bit = (c[ci[b]] > c[cj[b]]) ? 1 : 0;
        if (b < 64)        ha |= (uint64_t)bit << b;
        else /*b<128*/     hb |= (uint64_t)bit << (b - 64);
      }
      hashA[idx] = ha;
      if (Bloc > 64) hashB[idx] = hb;
    }
  }
};

// ----- simple blue-noise seeding on coarse grid -----------------------------
struct Seed {
  int gridIndex; // linear index in full grid
  int voxIdx;    // index in mask order (0..Nmask-1)
};

static void seed_blue_noise(const std::vector<int> &mask_lin, // Nmask, 0-based linear indices
                            int nx, int ny, int nz,
                            int K,
                            std::vector<Seed> &seeds) {
  const int64_t Nmask = (int64_t)mask_lin.size();
  if (Nmask == 0) return;

  const double S_est = std::cbrt((double)Nmask / (double)K);
  const int cx = std::max(1, (int)std::floor((double)nx / std::max(1.0, S_est)));
  const int cy = std::max(1, (int)std::floor((double)ny / std::max(1.0, S_est)));
  const int cz = std::max(1, (int)std::floor((double)nz / std::max(1.0, S_est)));

  // Coarse grid dims
  const int gx = cx, gy = cy, gz = cz;
  const int64_t gcells = (int64_t)gx * gy * gz;
  std::vector<int> coarse(gcells, -1);

  auto cell_of = [&](int idx)->int {
    int z = idx / (nx*ny);
    int rem = idx - z*(nx*ny);
    int y = rem / nx;
    int x = rem - y*nx;
    int ix = (int)((int64_t)x * gx / nx);
    int iy = (int)((int64_t)y * gy / ny);
    int iz = (int)((int64_t)z * gz / nz);
    if (ix >= gx) ix = gx-1; if (iy >= gy) iy = gy-1; if (iz >= gz) iz = gz-1;
    return (iz * gy + iy) * gx + ix;
  };

  // First pass: build cell membership lists to enable spatial-balanced selection
  std::vector<std::vector<int>> cell_members(gcells);
  for (size_t vi = 0; vi < mask_lin.size(); ++vi) {
    int idx = mask_lin[vi];
    int c = cell_of(idx);
    cell_members[c].push_back((int)vi);
  }

  // Select voxel closest to each cell's geometric center (eliminates z-ordering bias)
  std::vector<int> reps;
  reps.reserve(K*2);
  for (int64_t c = 0; c < gcells; ++c) {
    if (cell_members[c].empty()) continue;

    // Calculate coarse cell center in voxel coordinates
    int iz = (int)(c / (gx * gy));
    int rem = (int)(c - iz * (gx * gy));
    int iy = rem / gx;
    int ix = rem - iy * gx;
    double cx = ((double)ix + 0.5) * (double)nx / (double)gx;
    double cy = ((double)iy + 0.5) * (double)ny / (double)gy;
    double cz = ((double)iz + 0.5) * (double)nz / (double)gz;

    // Find voxel closest to cell center
    int best_vi = -1;
    double best_dist = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < cell_members[c].size(); ++i) {
      int vi = cell_members[c][i];
      int idx = mask_lin[vi];
      int zc = idx / (nx * ny);
      int r = idx - zc * (nx * ny);
      int yc = r / nx;
      int xc = r - yc * nx;
      double dx = (double)xc - cx;
      double dy = (double)yc - cy;
      double dz = (double)zc - cz;
      double d2 = dx*dx + dy*dy + dz*dz;
      if (d2 < best_dist) {
        best_dist = d2;
        best_vi = vi;
      }
    }
    if (best_vi >= 0) {
      reps.push_back(best_vi);
    }
  }

  // Adjust to exactly K seeds
  std::vector<int> chosen;
  chosen.reserve(K);
  if ((int)reps.size() >= K) {
    // simple striding selection for determinism
    double step = (double)reps.size() / (double)K;
    double acc = 0.0;
    for (int k = 0; k < K; ++k) {
      int pick = reps[(int)std::floor(acc)];
      chosen.push_back(pick);
      acc += step;
    }
  } else {
    std::vector<uint8_t> selected((size_t)Nmask, 0);
    for (size_t i = 0; i < reps.size(); ++i) {
      int pick = reps[i];
      if (!selected[(size_t)pick]) {
        chosen.push_back(pick);
        selected[(size_t)pick] = 1;
      }
    }
    // Deterministic top-up with distinct valid mask voxels. K <= Nmask is
    // validated at the native boundary, so this must reach exactly K.
    for (int pick = 0; pick < Nmask && (int)chosen.size() < K; ++pick) {
      if (!selected[(size_t)pick]) {
        chosen.push_back(pick);
        selected[(size_t)pick] = 1;
      }
    }
  }

  if ((int)chosen.size() != K) {
    stop("FLASH-3D internal error: could not create K distinct seeds");
  }

  seeds.clear();
  seeds.reserve(K);
  for (int k = 0; k < K; ++k) {
    int voxIdx = chosen[k];
    Seed s;
    s.voxIdx = voxIdx;
    s.gridIndex = mask_lin[(size_t)voxIdx];
    seeds.push_back(s);
  }
}

// ----- Parallel centroid reduction worker -----------------------------------
struct ParallelCentroidWorker : public Worker {
  const std::vector<int> &labels;
  const uint8_t *is_mask;
  const uint64_t *hashA, *hashB;
  const float *best;
  const int nx, ny;
  const int K;
  const int bits;

  // Accumulator structure for each cluster
  struct CentroidAccum {
    double sx, sy, sz;
    int64_t count;
    uint32_t bitcA[64];
    uint32_t bitcB[64];
    float max_err;
    int   max_idx;
    CentroidAccum() : sx(0), sy(0), sz(0), count(0), max_err(-std::numeric_limits<float>::infinity()), max_idx(-1) {
      memset(bitcA, 0, sizeof(bitcA));
      memset(bitcB, 0, sizeof(bitcB));
    }
  };

  std::vector<CentroidAccum> acc;  // Per-thread accumulator

  ParallelCentroidWorker(const std::vector<int> &lbl, const uint8_t *msk,
                         const uint64_t *ha, const uint64_t *hb, const float *best_,
                         int nx_, int ny_, int K_, int bits_)
    : labels(lbl), is_mask(msk), hashA(ha), hashB(hb), best(best_),
      nx(nx_), ny(ny_), K(K_), bits(bits_) {
    acc.resize(K);
  }

  // Split constructor for RcppParallel
  ParallelCentroidWorker(const ParallelCentroidWorker& other, Split)
    : labels(other.labels), is_mask(other.is_mask),
      hashA(other.hashA), hashB(other.hashB), best(other.best),
      nx(other.nx), ny(other.ny), K(other.K), bits(other.bits) {
    acc.resize(K);
  }

  void operator()(std::size_t begin, std::size_t end) {
    // OPTIMIZED: Use incremental coordinates like RelaxWorker
    int temp = (int)begin;
    int zc = temp / (nx * ny);
    temp %= (nx * ny);
    int yc = temp / nx;
    int xc = temp % nx;

    for (std::size_t i = begin; i < end; ++i) {
      int lab = labels[i];
      if (is_mask[i] && lab >= 0 && lab < K) {
        CentroidAccum &a = acc[lab];
        a.sx += xc;
        a.sy += yc;
        a.sz += zc;
        a.count++;

        float err = best ? best[i] : 0.0f;
        if (err > a.max_err) {
          a.max_err = err;
          a.max_idx = (int)i;
        }

        // Bit voting for hash majority
        uint64_t h = hashA[i];
        for (int b = 0; b < 64; ++b) {
          if ((h >> b) & 1) a.bitcA[b]++;
        }

        if (bits > 64) {
          h = hashB[i];
          for (int b = 0; b < 64; ++b) {
            if ((h >> b) & 1) a.bitcB[b]++;
          }
        }
      }

      // Increment coordinates
      xc++;
      if (xc == nx) {
        xc = 0;
        yc++;
        if (yc == ny) {
          yc = 0;
          zc++;
        }
      }
    }
  }

  // Join operation to merge thread-local accumulators
  void join(const ParallelCentroidWorker& rhs) {
    for (int k = 0; k < K; ++k) {
      acc[k].sx += rhs.acc[k].sx;
      acc[k].sy += rhs.acc[k].sy;
      acc[k].sz += rhs.acc[k].sz;
      acc[k].count += rhs.acc[k].count;

      for (int b = 0; b < 64; ++b) {
        acc[k].bitcA[b] += rhs.acc[k].bitcA[b];
      }

      if (bits > 64) {
        for (int b = 0; b < 64; ++b) {
          acc[k].bitcB[b] += rhs.acc[k].bitcB[b];
        }
      }

      if (rhs.acc[k].max_err > acc[k].max_err) {
        acc[k].max_err = rhs.acc[k].max_err;
        acc[k].max_idx = rhs.acc[k].max_idx;
      }
    }
  }
};

// ----- Exact feature center computation worker ------------------------------
struct ExactFeatureWorker : public Worker {
  const RMatrix<double> ts_mat;  // T x Nmask time series matrix
  const int *mask_lin;            // Mapping: mask index -> grid index
  const int *grid_labels;         // Grid labels (full grid size)
  const int Nmask, T, K;

  // Thread-local accumulators: K x T sums + K counts
  std::vector<double> local_sums;
  std::vector<int> local_counts;

  ExactFeatureWorker(const NumericMatrix &ts, const int *mlin, const int *glbl,
                     int Nm, int t, int k)
    : ts_mat(ts), mask_lin(mlin), grid_labels(glbl),
      Nmask(Nm), T(t), K(k) {
    local_sums.assign((size_t)K * T, 0.0);
    local_counts.assign(K, 0);
  }

  // Split constructor
  ExactFeatureWorker(const ExactFeatureWorker& other, Split)
    : ts_mat(other.ts_mat), mask_lin(other.mask_lin), grid_labels(other.grid_labels),
      Nmask(other.Nmask), T(other.T), K(other.K) {
    local_sums.assign((size_t)K * T, 0.0);
    local_counts.assign(K, 0);
  }

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      int grid_idx = mask_lin[i];
      int lbl = grid_labels[grid_idx];

      if (lbl >= 0 && lbl < K) {
        local_counts[lbl]++;
        // Access column i of ts_mat (T x Nmask, column-major)
        double *dest = &local_sums[(size_t)lbl * T];
        for (int t = 0; t < T; ++t) {
          dest[t] += ts_mat(t, (int)i);
        }
      }
    }
  }

  // Join operation
  void join(const ExactFeatureWorker& rhs) {
    for (int k = 0; k < K; ++k) {
      local_counts[k] += rhs.local_counts[k];
      for (int t = 0; t < T; ++t) {
        local_sums[(size_t)k * T + t] += rhs.local_sums[(size_t)k * T + t];
      }
    }
  }
};

// Reassign one distinct high-error voxel to every empty cluster. Donors must
// retain at least one voxel, so the returned labels contain exactly K groups.
static std::vector<int> reseed_empty_clusters(
    std::vector<int> &labels,
    std::vector<int> &carry,
    std::vector<float> &best,
    const std::vector<int> &mask_lin,
    int K) {
  std::vector<int64_t> counts(K, 0);
  for (size_t pos = 0; pos < mask_lin.size(); ++pos) {
    int label = labels[(size_t)mask_lin[pos]];
    if (label >= 0 && label < K) counts[(size_t)label]++;
  }

  std::vector<int> candidates = mask_lin;
  std::stable_sort(candidates.begin(), candidates.end(), [&](int lhs, int rhs) {
    float lhs_score = best[(size_t)lhs];
    float rhs_score = best[(size_t)rhs];
    if (lhs_score != rhs_score) return lhs_score > rhs_score;
    return lhs < rhs;
  });

  std::vector<uint8_t> used(best.size(), 0);
  std::vector<int> reseed_idx(K, -1);
  for (int k = 0; k < K; ++k) {
    if (counts[(size_t)k] != 0) continue;
    for (size_t candidate_pos = 0; candidate_pos < candidates.size(); ++candidate_pos) {
      int candidate = candidates[candidate_pos];
      int donor = labels[(size_t)candidate];
      if (!used[(size_t)candidate] && donor >= 0 && donor < K &&
          counts[(size_t)donor] > 1) {
        counts[(size_t)donor]--;
        counts[(size_t)k] = 1;
        labels[(size_t)candidate] = k;
        carry[(size_t)candidate] = k;
        best[(size_t)candidate] = 0.0f;
        used[(size_t)candidate] = 1;
        reseed_idx[(size_t)k] = candidate;
        break;
      }
    }
    if (reseed_idx[(size_t)k] < 0) {
      stop("FLASH-3D could not find a distinct valid donor for an empty cluster");
    }
  }
  return reseed_idx;
}

} // namespace

// Testable boundary for the same distinct-reseed operation used by the core.
// Labels are one-based in R and every vector position is treated as masked.
// [[Rcpp::export]]
List flash3d_reseed_empty_cpp(IntegerVector labels, NumericVector best, int K) {
  if (labels.size() != best.size() || labels.size() < K || K < 1) {
    stop("labels and best must have equal length with K between 1 and N");
  }
  const int n = labels.size();
  std::vector<int> native_labels(n), carry(n), mask_lin(n);
  std::vector<float> native_best(n);
  for (int i = 0; i < n; ++i) {
    if (IntegerVector::is_na(labels[i]) || labels[i] < 1 || labels[i] > K ||
        !R_FINITE(best[i])) {
      stop("labels must be in 1..K and best must be finite");
    }
    native_labels[(size_t)i] = labels[i] - 1;
    carry[(size_t)i] = native_labels[(size_t)i];
    native_best[(size_t)i] = (float)best[i];
    mask_lin[(size_t)i] = i;
  }
  std::vector<int> reseeds = reseed_empty_clusters(
    native_labels, carry, native_best, mask_lin, K
  );
  IntegerVector output_labels(n);
  IntegerVector output_reseeds;
  for (int i = 0; i < n; ++i) output_labels[i] = native_labels[(size_t)i] + 1;
  for (int k = 0; k < K; ++k) {
    if (reseeds[(size_t)k] >= 0) output_reseeds.push_back(reseeds[(size_t)k] + 1);
  }
  return List::create(
    Named("labels") = output_labels,
    Named("reseed_indices") = output_reseeds
  );
}

// ----- FLASH-3D clustering core ---------------------------------------------

// [[Rcpp::export]]
List flash3d_supervoxels_cpp(NumericMatrix ts,     // T x Nmask (time series for masked voxels)
                                      IntegerVector mask_lin0, // Nmask linear indices (1-based from R)
                                      IntegerVector dims,      // c(nx,ny,nz)
                                      int K,
                                      NumericVector lambda,    // c(lambda_s, lambda_t, lambda_g)
                                      int rounds = 2,
                                      int bits = 64,
                                      int dctM = 12,
                                      NumericVector vox_scale = NumericVector::create(1.0,1.0,1.0),
                                      Nullable<NumericVector> barrier_opt = R_NilValue,
                                      bool verbose = false,
                                      bool diagnostics = false) {

  if (dims.size() != 3) stop("dims must contain three finite integers");
  for (int axis = 0; axis < 3; ++axis) {
    if (IntegerVector::is_na(dims[axis])) stop("dims must contain three finite integers");
  }
  if (dims[0] <= 0 || dims[1] <= 0 || dims[2] <= 0) stop("dims must be positive");
  if (bits != 64 && bits != 128) stop("bits must be 64 or 128");
  if (dctM < 4 || dctM > 32) stop("dctM must be in [4,32]");
  if (rounds < 1) stop("rounds must be positive");
  if (lambda.size() < 2) stop("lambda must be at least length 2: c(lambda_s, lambda_t, [lambda_g])");
  const double lambda_s0 = lambda[0];
  const double lambda_t  = lambda[1];
  const double lambda_g  = (lambda.size() >= 3) ? lambda[2] : 0.0;
  if (!R_FINITE(lambda_s0) || !R_FINITE(lambda_t) || !R_FINITE(lambda_g) ||
      lambda_s0 < 0.0 || lambda_t < 0.0 || lambda_g < 0.0) {
    stop("lambda weights must be finite and non-negative");
  }

  const int nx = dims[0], ny = dims[1], nz = dims[2];
  const size_t Ngrid = (size_t)nx * ny * nz;
  const int Nmask = ts.ncol();
  const int T = ts.nrow();
  const int dctM_eff = std::min(std::min(dctM, T), 32);
  if (Nmask < 1 || T < 1) stop("ts must have at least one row and one column");
  if (K < 1 || K > Nmask) stop("K must be between 1 and the number of masked voxels");
  if (Nmask != mask_lin0.size()) stop("ts ncol must equal length(mask_lin)");
  if ((int)vox_scale.size() != 3) stop("vox_scale must be length 3");
  for (int axis = 0; axis < 3; ++axis) {
    if (!R_FINITE(vox_scale[axis]) || vox_scale[axis] <= 0.0) {
      stop("vox_scale values must be finite and positive");
    }
  }
  for (int col = 0; col < Nmask; ++col) {
    for (int row = 0; row < T; ++row) {
      if (!R_FINITE(ts(row, col))) stop("ts values must be finite");
    }
  }

  // Convert 1-based to 0-based linear indices
  std::vector<int> mask_lin(Nmask);
  std::vector<uint8_t> seen_mask(Ngrid, 0);
  for (int i = 0; i < Nmask; ++i) {
    if (IntegerVector::is_na(mask_lin0[i])) stop("mask_lin contains NA");
    mask_lin[i] = mask_lin0[i] - 1;
    if (mask_lin[i] < 0 || (size_t)mask_lin[i] >= Ngrid) {
      stop("mask_lin index is out of bounds");
    }
    if (seen_mask[(size_t)mask_lin[i]]) stop("mask_lin indices must be distinct");
    seen_mask[(size_t)mask_lin[i]] = 1;
  }

  // Optional barrier
  std::vector<float> barrier;
  bool use_barrier = false;
  if (barrier_opt.isNotNull()) {
    NumericVector b(barrier_opt);
    if ((int)b.size() != (int)Ngrid) stop("barrier must be length nx*ny*nz (full grid)");
    barrier.resize(Ngrid);
    for (size_t i = 0; i < Ngrid; ++i) {
      if (!R_FINITE(b[i]) || b[i] < 0.0) {
        stop("barrier values must be finite and non-negative");
      }
      barrier[i] = (float)b[i];
    }
    use_barrier = (lambda_g > 0.0);
  } else if (lambda_g > 0.0) {
    stop("positive lambda_g requires a barrier");
  }

  // Build mask bitmap
  std::vector<uint8_t> is_mask(Ngrid, 0);
  for (int i = 0; i < Nmask; ++i) {
    int idx = mask_lin[i];
    if (idx >= 0 && (size_t)idx < Ngrid) is_mask[(size_t)idx] = 1;
  }

  // Estimate spacing S
  const double S_est = std::cbrt((double)Nmask / (double)K);
  const double S2 = std::max(1.0, S_est * S_est);

  // Precompute DCT table & comparator pairs
  std::vector<double> ctab;
  build_dct_table(T, dctM_eff, ctab);
  std::vector<uint8_t> ci, cj;
  make_comparators(dctM_eff, bits, ci, cj);

  // Hash time series
  std::vector<uint64_t> hashA((size_t)Nmask, 0ull), hashB;
  if (bits == 128) hashB.assign((size_t)Nmask, 0ull);
  NumericVector mean_tail(Nmask), var_tail(Nmask);
  {
    HashWorker hw(ts, dctM_eff, bits, ctab, ci, cj, mean_tail, var_tail, hashA, hashB);
    parallelFor(0, (size_t)Nmask, hw, (size_t)2048);
  }

  // Place hashes into full grid arrays
  std::vector<uint64_t> gridHashA(Ngrid, 0ull), gridHashB;
  if (bits == 128) gridHashB.assign(Ngrid, 0ull);
  for (int i = 0; i < Nmask; ++i) {
    int idx = mask_lin[i];
    gridHashA[(size_t)idx] = hashA[(size_t)i];
    if (bits == 128) gridHashB[(size_t)idx] = hashB[(size_t)i];
  }

  // Seed K sites
  std::vector<Seed> seeds;
  seed_blue_noise(mask_lin, nx, ny, nz, K, seeds);

  // Site state
  struct Site { float x,y,z; uint64_t ha, hb; int seed_idx; };
  std::vector<Site> sites((size_t)K);
  // Map from masked voxel index to coordinate
  auto idx_to_xyz = [&](int idx, float &x, float &y, float &z){
    int zc = idx / (nx*ny);
    int rem = idx - zc*(nx*ny);
    int yc = rem / nx;
    int xc = rem - yc*nx;
    x = (float)xc; y = (float)yc; z = (float)zc;
  };
  for (int k = 0; k < K; ++k) {
    int gridIdx = seeds[(size_t)k].gridIndex;
    float x,y,z;
    idx_to_xyz(gridIdx, x,y,z);
    sites[(size_t)k] = Site{ x, y, z, gridHashA[(size_t)gridIdx], (bits==128)?gridHashB[(size_t)gridIdx]:0ull, gridIdx };
  }


  // ---- Relaxation worker (thread-safe, embeds scoring logic; no std::function) ----
  struct RelaxWorker : public RcppParallel::Worker {
    // grid and dims
    const uint8_t *is_mask;
    const uint64_t *gridHashA;
    const uint64_t *gridHashB;
    const float *barrier;
    const bool has_hashB;
    const bool use_bar;
    const int nx, ny, nz;
    const double vx_sx, vx_sy, vx_sz; // voxel scaling (dx,dy,dz)

    // params
    const double lambda_s_cur;
    const double lambda_t;
    const double lambda_g;
    const double S2;

    // site state (read-only)
    const std::vector<Site> &sites;

    // double-buffered fields
    const std::vector<int>   &label_read;
    const std::vector<float> &best_read;
    const std::vector<int>   &carry_read;
          std::vector<int>   &label_write;
          std::vector<float> &best_write;
          std::vector<int>   &carry_write;

    const int step; // stride

    RelaxWorker(const uint8_t *is_mask_,
                const uint64_t *gridHashA_,
                const uint64_t *gridHashB_,
                const float *barrier_,
                bool has_hashB_,
                bool use_bar_,
                int nx_, int ny_, int nz_,
                double vx_sx_, double vx_sy_, double vx_sz_,
                double lambda_s_cur_, double lambda_t_,
                double lambda_g_, double S2_,
                const std::vector<Site> &sites_,
                const std::vector<int>   &label_read_,
                const std::vector<float> &best_read_,
                const std::vector<int>   &carry_read_,
                      std::vector<int>   &label_write_,
                      std::vector<float> &best_write_,
                      std::vector<int>   &carry_write_,
                int step_)
      : is_mask(is_mask_), gridHashA(gridHashA_), gridHashB(gridHashB_), barrier(barrier_),
        has_hashB(has_hashB_), use_bar(use_bar_),
        nx(nx_), ny(ny_), nz(nz_),
        vx_sx(vx_sx_), vx_sy(vx_sy_), vx_sz(vx_sz_),
        lambda_s_cur(lambda_s_cur_), lambda_t(lambda_t_), lambda_g(lambda_g_), S2(S2_),
        sites(sites_),
        label_read(label_read_), best_read(best_read_), carry_read(carry_read_),
        label_write(label_write_), best_write(best_write_), carry_write(carry_write_),
        step(step_) {}

    inline uint32_t hamming_dist(size_t vIdx, const Site &s) const {
      uint64_t xa = gridHashA[vIdx] ^ s.ha;
      uint32_t hd = (uint32_t)popcount64(xa);
      if (has_hashB) {
        uint64_t xb = gridHashB[vIdx] ^ s.hb;
        hd += (uint32_t)popcount64(xb);
      }
      return hd;
    }

    // OPTIMIZED: Pass coordinates directly to avoid expensive div/mod operations
    inline float score_vox(int vIdx, int xc, int yc, int zc, int siteId) const {
      if (siteId < 0) return std::numeric_limits<float>::infinity();
      const Site &s = sites[(size_t)siteId];

      double dx = ((double)xc - (double)s.x) * vx_sx;
      double dy = ((double)yc - (double)s.y) * vx_sy;
      double dz = ((double)zc - (double)s.z) * vx_sz;
      double sd2 = dx*dx + dy*dy + dz*dz;

      uint32_t hd = hamming_dist((size_t)vIdx, s);

      double sc = lambda_s_cur * (sd2 / S2) +
        lambda_t * ((double)hd / (double)(has_hashB ? 128 : 64));
      if (use_bar) {
        // Candidate-dependent boundary-field penalty. The previous additive
        // barrier[vIdx] term was constant across sites and could never change
        // an assignment. A site now pays for crossing to a different barrier
        // level, symmetrically in either direction.
        sc += lambda_g * std::fabs(
          (double)barrier[(size_t)vIdx] -
          (double)barrier[(size_t)s.seed_idx]
        );
      }
      return (float)sc;
    }

    void operator()(std::size_t begin, std::size_t end) {
      // OPTIMIZED: Compute coordinates once for 'begin', then increment
      int temp = (int)begin;
      int zc = temp / (nx * ny);
      temp %= (nx * ny);
      int yc = temp / nx;
      int xc = temp % nx;

      for (size_t idx = begin; idx < end; ++idx) {
        if (!is_mask[idx]) {
          best_write[idx] = best_read[idx];
          label_write[idx] = label_read[idx];
          carry_write[idx] = carry_read[idx];
        } else {
          float b = best_read[idx];
          int lab = label_read[idx];
          int car = carry_read[idx];

          // Check 26 neighbors using current coordinates
          for (int dz = -1; dz <= 1; ++dz) {
            int zn = zc + dz * step;
            if (zn < 0 || zn >= nz) continue;

            for (int dy = -1; dy <= 1; ++dy) {
              int yn = yc + dy * step;
              if (yn < 0 || yn >= ny) continue;

              for (int dx = -1; dx <= 1; ++dx) {
                if (dx==0 && dy==0 && dz==0) continue;
                int xn = xc + dx * step;
                if (xn < 0 || xn >= nx) continue;

                size_t nb = (size_t)(zn * (nx*ny) + yn * nx + xn);
                if (!is_mask[nb]) continue;

                int cand = carry_read[nb];
                if (cand < 0) continue;

                // Pass current coordinates to avoid div/mod in score_vox
                float sc = score_vox((int)idx, xc, yc, zc, cand);
                if (sc < b) { b = sc; lab = cand; car = cand; }
              }
            }
          }
          best_write[idx]  = b;
          label_write[idx] = lab;
          carry_write[idx] = car;
        }

        // OPTIMIZED: Increment coordinates for next voxel
        xc++;
        if (xc == nx) {
          xc = 0;
          yc++;
          if (yc == ny) {
            yc = 0;
            zc++;
          }
        }
      }
    }
  };

    // Working grids (double buffered)
  const float inf = std::numeric_limits<float>::infinity();
  std::vector<int>   label(Ngrid, -1), label_next(Ngrid, -1);
  std::vector<float> best (Ngrid,  inf), best_next (Ngrid, inf);
  std::vector<int>   carry(Ngrid, -1), carry_next(Ngrid, -1); // site id carried by cell

  double lambda_s_current = lambda_s0;
  IntegerVector round_reseed_counts(rounds);
  std::vector<int> reseed_history;
  std::vector<Site> final_assignment_sites;

  // Helper lambda for power-of-2 ceiling
  auto ceil_pow2 = [](int v)->int {
    int p = 1; while (p < v) p <<= 1; return p;
  };

  // ---- outer rounds ----
  for (int r = 0; r < rounds; ++r) {
    if (verbose) Rcpp::Rcout << "[FLASH-3D] Round " << (r+1) << " / " << rounds << " [JFA]\n";

    // Anneal spatial weight slightly upward per round (encourages compactness to finish)
    lambda_s_current = lambda_s0 * (1.0 + 0.5 * (double)r);
    if (r == rounds - 1) final_assignment_sites = sites;

    // reset grids and seed from current medoids
    std::fill(best.begin(),  best.end(),  inf);
    std::fill(label.begin(), label.end(), -1);
    std::fill(carry.begin(), carry.end(), -1);

    for (int k = 0; k < K; ++k) {
      int g = sites[(size_t)k].seed_idx;
      if (g < 0 || (size_t)g >= Ngrid) continue;
      label[(size_t)g] = k;
      carry[(size_t)g] = k;
      best [(size_t)g] = 0.0f;
    }

    // mirror into next buffers before first swap
    best_next  = best;
    label_next = label;
    carry_next = carry;

    // jump-flood steps: from large to 1
    int maxdim = std::max(nx, std::max(ny, nz));
    for (int s = ceil_pow2(maxdim); s >= 1; s >>= 1) {
      // swap buffers
      std::swap(best, best_next);
      std::swap(label, label_next);
      std::swap(carry, carry_next);

      {
        RelaxWorker worker(
          is_mask.data(),
          gridHashA.data(),
          (bits==128 ? gridHashB.data() : nullptr),
          (use_barrier ? barrier.data() : nullptr),
          (bits==128),
          use_barrier,
          nx, ny, nz,
          (double)vox_scale[0], (double)vox_scale[1], (double)vox_scale[2],
          lambda_s_current, lambda_t, lambda_g, S2,
          sites,
          label, best, carry,
          label_next, best_next, carry_next,
          s
        );
        parallelFor(0, (size_t)Ngrid, worker, (size_t)65536);
      }
    }

    // ensure we use the latest buffers after the last JFA step
    std::swap(best, best_next);
    std::swap(label, label_next);
    std::swap(carry, carry_next);

    // ----- recenter sites (OPTIMIZED: parallel reduction) -----
    ParallelCentroidWorker centroid_worker(
      label,
      is_mask.data(),
      gridHashA.data(),
      (bits == 128) ? gridHashB.data() : nullptr,
      best.data(),
      nx, ny, K, bits
    );
    parallelReduce(0, Ngrid, centroid_worker);

    // compute medoids (closest voxel to centroid per cluster)
    std::vector<int> medoid_idx(K, -1);
    std::vector<double> medoid_best(K, std::numeric_limits<double>::infinity());
    for (int m = 0; m < Nmask; ++m) {
      int idx = mask_lin[m];
      int lab = label[(size_t)idx];
      if (lab < 0 || lab >= K) continue;
      double cnt = (double)centroid_worker.acc[(size_t)lab].count;
      if (cnt <= 0) continue;
      double cx = centroid_worker.acc[(size_t)lab].sx / cnt;
      double cy = centroid_worker.acc[(size_t)lab].sy / cnt;
      double cz = centroid_worker.acc[(size_t)lab].sz / cnt;
      int zc = idx / (nx*ny);
      int rem = idx - zc*(nx*ny);
      int yc = rem / nx;
      int xc = rem - yc*nx;
      double dx = (double)xc - cx;
      double dy = (double)yc - cy;
      double dz = (double)zc - cz;
      double d2 = dx*dx + dy*dy + dz*dz;
      if (d2 < medoid_best[(size_t)lab]) {
        medoid_best[(size_t)lab] = d2;
        medoid_idx[(size_t)lab]  = idx;
      }
    }

    // Reseed every empty cluster from a different high-error voxel. Donors are
    // required to retain at least one voxel, so the reassignment preserves K
    // immediately, including after the final round.
    std::vector<int> reseed_idx = reseed_empty_clusters(
      label, carry, best, mask_lin, K
    );
    int reseeded = 0;
    for (int k = 0; k < K; ++k) {
      if (reseed_idx[(size_t)k] >= 0) {
        reseed_history.push_back(reseed_idx[(size_t)k] + 1);
        reseeded++;
      }
    }
    round_reseed_counts[r] = reseeded;

    for (int k = 0; k < K; ++k) {
      uint64_t ha = 0ull, hb = 0ull;
      if (reseed_idx[(size_t)k] >= 0) {
        int replacement = reseed_idx[(size_t)k];
        float x,y,z; idx_to_xyz(replacement, x,y,z);
        ha = gridHashA[(size_t)replacement];
        if (bits==128) hb = gridHashB[(size_t)replacement];
        sites[(size_t)k] = Site{ x, y, z, ha, hb, replacement };
        continue;
      }

      int seed_idx = medoid_idx[(size_t)k];
      if (seed_idx < 0) seed_idx = sites[(size_t)k].seed_idx; // fallback
      float x,y,z; idx_to_xyz(seed_idx, x,y,z);

      for (int b=0;b<64;b++) {
        if (centroid_worker.acc[(size_t)k].bitcA[b] * 2 >= centroid_worker.acc[(size_t)k].count) {
          ha |= (1ull<<b);
        }
      }
      if (bits==128) {
        for (int b=0;b<64;b++) {
          if (centroid_worker.acc[(size_t)k].bitcB[b] * 2 >= centroid_worker.acc[(size_t)k].count) {
            hb |= (1ull<<b);
          }
        }
      }
      sites[(size_t)k] = Site{ x, y, z, ha, hb, seed_idx };
    }
  }

  // Pull labels for mask voxels (1-based cluster ids to be friendly in R)
  IntegerVector lab_mask(Nmask);
  std::vector<double> coord_sum_x(K, 0.0), coord_sum_y(K, 0.0), coord_sum_z(K, 0.0);
  std::vector<int> final_counts(K, 0);
  for (int i = 0; i < Nmask; ++i) {
    int idx = mask_lin[i];
    int lab = label[(size_t)idx];
    if (lab < 0 || lab >= K) stop("FLASH-3D left a masked voxel unassigned");
    lab_mask[i] = lab + 1;
    int zc = idx / (nx * ny);
    int rem = idx - zc * (nx * ny);
    int yc = rem / nx;
    int xc = rem - yc * nx;
    coord_sum_x[(size_t)lab] += xc;
    coord_sum_y[(size_t)lab] += yc;
    coord_sum_z[(size_t)lab] += zc;
    final_counts[(size_t)lab]++;
  }

  // ----- OPTIMIZED: Compute exact feature centers in C++ (Phase 3) -----
  if (verbose) Rcpp::Rcout << "[FLASH-3D] Computing final feature centers in C++...\n";

  ExactFeatureWorker feature_worker(ts, mask_lin.data(), label.data(), Nmask, T, K);
  parallelReduce(0, (size_t)Nmask, feature_worker);

  // Build output matrices
  NumericMatrix centers(K, T);
  NumericMatrix coords(K, 3);

  for (int k = 0; k < K; ++k) {
    double c = (double)feature_worker.local_counts[k];
    if (c > 0 && final_counts[(size_t)k] > 0) {
      // Feature space centers (average time series)
      for (int t = 0; t < T; ++t) {
        centers(k, t) = feature_worker.local_sums[(size_t)k * T + t] / c;
      }
    } else {
      stop("FLASH-3D internal error: final cluster is empty after reseeding");
    }

    // Exact zero-based voxel-coordinate mean of the final assignments. The R
    // boundary independently converts final labels to physical millimeters.
    coords(k, 0) = coord_sum_x[(size_t)k] / final_counts[(size_t)k];
    coords(k, 1) = coord_sum_y[(size_t)k] / final_counts[(size_t)k];
    coords(k, 2) = coord_sum_z[(size_t)k] / final_counts[(size_t)k];
  }

  List result = List::create(
    Named("labels") = lab_mask,
    Named("centers") = centers,
    Named("coords") = coords,
    Named("K") = K,
    Named("reseed_count") = reseed_history.size(),
    Named("round_reseed_counts") = round_reseed_counts,
    Named("reseed_indices_grid") = wrap(reseed_history)
  );

  if (diagnostics) {
    NumericMatrix assignment_site_coords(K, 3);
    IntegerVector assignment_site_seed_grid(K);
    IntegerMatrix assignment_site_hash_bits(K, bits);
    IntegerMatrix voxel_hash_bits(Nmask, bits);
    for (int k = 0; k < K; ++k) {
      const Site &site = final_assignment_sites[(size_t)k];
      assignment_site_coords(k, 0) = site.x;
      assignment_site_coords(k, 1) = site.y;
      assignment_site_coords(k, 2) = site.z;
      assignment_site_seed_grid[k] = site.seed_idx + 1;
      for (int bit = 0; bit < bits; ++bit) {
        uint64_t word = bit < 64 ? site.ha : site.hb;
        assignment_site_hash_bits(k, bit) = (word >> (bit % 64)) & 1ull;
      }
    }
    for (int voxel = 0; voxel < Nmask; ++voxel) {
      for (int bit = 0; bit < bits; ++bit) {
        uint64_t word = bit < 64 ? hashA[(size_t)voxel] : hashB[(size_t)voxel];
        voxel_hash_bits(voxel, bit) = (word >> (bit % 64)) & 1ull;
      }
    }
    result["assignment_site_coords"] = assignment_site_coords;
    result["assignment_site_seed_grid"] = assignment_site_seed_grid;
    result["assignment_site_hash_bits"] = assignment_site_hash_bits;
    result["voxel_hash_bits"] = voxel_hash_bits;
    result["S2"] = S2;
    result["lambda_s_final"] = lambda_s_current;
  }
  return result;
}
