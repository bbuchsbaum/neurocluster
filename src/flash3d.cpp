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
  // portable fallback
  uint32_t c = 0;
  while (x) { x &= (x - 1); ++c; }
  return c;
#endif
}

// ----- fixed comparator table (deterministic) -------------------------------
static void make_comparators(int M, int B, std::vector<uint8_t> &ci, std::vector<uint8_t> &cj) {
  ci.resize(B); cj.resize(B);
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

  // First pass: pick one representative per non-empty coarse cell
  std::vector<int> reps;
  reps.reserve(K*2);
  for (size_t vi = 0; vi < mask_lin.size(); ++vi) {
    int idx = mask_lin[vi];
    int c = cell_of(idx);
    if (coarse[c] < 0) {
      coarse[c] = (int)vi;
      reps.push_back((int)vi);
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
    chosen = reps;
    // pad with evenly spaced additional mask voxels
    double step = (double)mask_lin.size() / (double)std::max(1, K - (int)reps.size());
    double acc = 0.0;
    while ((int)chosen.size() < K) {
      int pick = (int)std::floor(acc);
      chosen.push_back(pick % (int)mask_lin.size());
      acc += step;
    }
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

// ----- FLASH-3D clustering core ---------------------------------------------

// [[Rcpp::export]]
IntegerVector flash3d_supervoxels_cpp(NumericMatrix ts,     // T x Nmask (time series for masked voxels)
                                      IntegerVector mask_lin0, // Nmask linear indices (1-based from R)
                                      IntegerVector dims,      // c(nx,ny,nz)
                                      int K,
                                      NumericVector lambda,    // c(lambda_s, lambda_t, lambda_g)
                                      int rounds = 2,
                                      int bits = 64,
                                      int dctM = 12,
                                      NumericVector vox_scale = NumericVector::create(1.0,1.0,1.0),
                                      Nullable<NumericVector> barrier_opt = R_NilValue,
                                      bool verbose = false) {

  if (bits != 64 && bits != 128) stop("bits must be 64 or 128");
  if (dctM < 4 || dctM > 32) stop("dctM must be in [4,32]");
  if (lambda.size() < 2) stop("lambda must be at least length 2: c(lambda_s, lambda_t, [lambda_g])");
  const double lambda_s0 = lambda[0];
  const double lambda_t  = lambda[1];
  const double lambda_g  = (lambda.size() >= 3) ? lambda[2] : 0.0;

  const int nx = dims[0], ny = dims[1], nz = dims[2];
  const size_t Ngrid = (size_t)nx * ny * nz;
  const int Nmask = ts.ncol();
  const int T = ts.nrow();
  if (Nmask != mask_lin0.size()) stop("ts ncol must equal length(mask_lin)");
  if ((int)vox_scale.size() != 3) stop("vox_scale must be length 3");

  // Convert 1-based to 0-based linear indices
  std::vector<int> mask_lin(Nmask);
  for (int i = 0; i < Nmask; ++i) mask_lin[i] = mask_lin0[i] - 1;

  // Optional barrier
  std::vector<float> barrier;
  bool use_barrier = false;
  if (barrier_opt.isNotNull()) {
    NumericVector b(barrier_opt);
    if ((int)b.size() != (int)Ngrid) stop("barrier must be length nx*ny*nz (full grid)");
    barrier.resize(Ngrid);
    for (size_t i = 0; i < Ngrid; ++i) barrier[i] = (float)b[i];
    use_barrier = (lambda_g > 0.0);
  }

  // Build mask bitmap
  std::vector<uint8_t> is_mask(Ngrid, 0);
  for (int i = 0; i < Nmask; ++i) {
    int idx = mask_lin[i];
    if (idx >= 0 && (size_t)idx < Ngrid) is_mask[(size_t)idx] = 1;
  }

  // Estimate spacing S
  const double S_est = std::cbrt((double)Nmask / (double)K);
  const double S2 = S_est * S_est;

  // Precompute DCT table & comparator pairs
  std::vector<double> ctab;
  build_dct_table(T, dctM, ctab);
  std::vector<uint8_t> ci, cj;
  make_comparators(dctM, bits, ci, cj);

  // Hash time series
  std::vector<uint64_t> hashA((size_t)Nmask, 0ull), hashB;
  if (bits == 128) hashB.assign((size_t)Nmask, 0ull);
  NumericVector mean_tail(Nmask), var_tail(Nmask);
  {
    HashWorker hw(ts, dctM, bits, ctab, ci, cj, mean_tail, var_tail, hashA, hashB);
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
  struct Site { float x,y,z; uint64_t ha, hb; };
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
    sites[(size_t)k] = Site{ x, y, z, gridHashA[(size_t)gridIdx], (bits==128)?gridHashB[(size_t)gridIdx]:0ull };
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

    inline float score_vox(int vIdx, int siteId) const {
      if (siteId < 0) return std::numeric_limits<float>::infinity();
      const Site &s = sites[(size_t)siteId];
      // coords
      int zc = vIdx / (nx*ny);
      int rem = vIdx - zc*(nx*ny);
      int yc = rem / nx;
      int xc = rem - yc*nx;

      double dx = ((double)xc - (double)s.x) * vx_sx;
      double dy = ((double)yc - (double)s.y) * vx_sy;
      double dz = ((double)zc - (double)s.z) * vx_sz;
      double sd2 = dx*dx + dy*dy + dz*dz;

      uint32_t hd = hamming_dist((size_t)vIdx, s);

      double sc = lambda_s_cur * (sd2 / S2) + lambda_t * ((double)hd / (double)(has_hashB ? 128 : 64));
      if (use_bar) sc += lambda_g * barrier[(size_t)vIdx];
      return (float)sc;
    }

    void operator()(std::size_t begin, std::size_t end) {
      const size_t plane = (size_t)nx * (size_t)ny;
      for (size_t idx = begin; idx < end; ++idx) {
        if (!is_mask[idx]) {
          best_write[idx] = best_read[idx];
          label_write[idx] = label_read[idx];
          carry_write[idx] = carry_read[idx];
          continue;
        }
        float b = best_read[idx];
        int lab = label_read[idx];
        int car = carry_read[idx];

        int id = (int)idx;
        int zc = id / (nx*ny);
        int rem = id - zc*(nx*ny);
        int yc = rem / nx;
        int xc = rem - yc*nx;

        for (int dz = -1; dz <= 1; ++dz) {
          int zn = zc + dz * step; if (zn < 0 || zn >= nz) continue;
          for (int dy = -1; dy <= 1; ++dy) {
            int yn = yc + dy * step; if (yn < 0 || yn >= ny) continue;
            for (int dx = -1; dx <= 1; ++dx) {
              if (dx==0 && dy==0 && dz==0) continue;
              int xn = xc + dx * step; if (xn < 0 || xn >= nx) continue;
              size_t nb = (size_t)(zn * (nx*ny) + yn * nx + xn);
              if (!is_mask[nb]) continue;
              int cand = carry_read[nb];
              if (cand < 0) continue;
              float sc = score_vox((int)idx, cand);
              if (sc < b) { b = sc; lab = cand; car = cand; }
            }
          }
        }
        best_write[idx]  = b;
        label_write[idx] = lab;
        carry_write[idx] = car;
      }
    }
  };

    // Working grids (double buffered)
  const float inf = std::numeric_limits<float>::infinity();
  std::vector<int>   label(Ngrid, -1), label_next(Ngrid, -1);
  std::vector<float> best (Ngrid,  inf), best_next (Ngrid, inf);
  std::vector<int>   carry(Ngrid, -1), carry_next(Ngrid, -1); // site id carried by cell

  // Initialize with seeds
  for (int k = 0; k < K; ++k) {
    int g = seeds[(size_t)k].gridIndex;
    label[(size_t)g] = k;
    carry[(size_t)g] = k;
    best [(size_t)g] = 0.0f;
  }


  double lambda_s_current = lambda_s0;
  auto score_vox = [&](int vIdx, int siteId)->float {
    if (siteId < 0) return inf;
    const Site &s = sites[(size_t)siteId];
    // coords
    int zc = vIdx / (nx*ny);
    int rem = vIdx - zc*(nx*ny);
    int yc = rem / nx;
    int xc = rem - yc*nx;

    double dx = ((double)xc - (double)s.x) * (double)vox_scale[0];
    double dy = ((double)yc - (double)s.y) * (double)vox_scale[1];
    double dz = ((double)zc - (double)s.z) * (double)vox_scale[2];
    double sd2 = dx*dx + dy*dy + dz*dz;
    // hamming
    uint64_t ha = gridHashA[(size_t)vIdx] ^ s.ha;
    uint32_t hd = popcount64(ha);
    if (bits == 128) {
      uint64_t hb = gridHashB[(size_t)vIdx] ^ s.hb;
      hd += popcount64(hb);
    }

    double sc = lambda_s_current * (sd2 / S2) + lambda_t * ((double)hd / (double)bits);
    if (use_barrier) sc += lambda_g * barrier[(size_t)vIdx];
    return (float)sc;
  };

  // Create a worker class for parallel relaxation
  
  auto ceil_pow2 = [](int v)->int {
    int p = 1; while (p < v) p <<= 1; return p;
  };

  // ---- outer rounds ----
  for (int r = 0; r < rounds; ++r) {
    if (verbose) Rcpp::Rcout << "[FLASH-3D] Round " << (r+1) << " / " << rounds << " [JFA]\n";

    // Anneal spatial weight slightly upward per round (encourages compactness to finish)
    lambda_s_current = lambda_s0 * (1.0 + 0.5 * (double)r);

    // sync best/label/carry first
    std::copy(best.begin(), best.end(), best_next.begin());
    std::copy(label.begin(), label.end(), label_next.begin());
    std::copy(carry.begin(), carry.end(), carry_next.begin());

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

    // ----- recenter sites (parallel reduction) -----
    struct Accum {
      double sx, sy, sz;
      int64_t count;
      uint32_t bitcA[64];
      uint32_t bitcB[64];
      Accum() : sx(0),sy(0),sz(0),count(0) { memset(bitcA,0,sizeof(bitcA)); memset(bitcB,0,sizeof(bitcB)); }
    };
    std::vector<Accum> acc((size_t)K);

    // Sequential accumulation to avoid race conditions
    // (parallel reduction would need thread-local accumulators and merging)
    for (size_t i = 0; i < Ngrid; ++i) {
      int lab = label[i];
      if (lab < 0 || lab >= K) continue;
      if (!is_mask[i]) continue;
      // coords
      int id = (int)i;
      int zc = id / (nx*ny);
      int rem = id - zc*(nx*ny);
      int yc = rem / nx;
      int xc = rem - yc*nx;
      Accum &A = acc[(size_t)lab];
      A.sx += (double)xc; A.sy += (double)yc; A.sz += (double)zc; A.count++;
      uint64_t ha = gridHashA[i];
      for (int b=0;b<64;b++) if (ha & (1ull<<b)) A.bitcA[b]++;
      if (bits==128) {
        uint64_t hb = gridHashB[i];
        for (int b=0;b<64;b++) if (hb & (1ull<<b)) A.bitcB[b]++;
      }
    }

    // update sites (handle empties by stealing from biggest cluster center)
    // find biggest cluster
    int big_k = 0; int64_t big_c = 0;
    for (int k=0;k<K;k++) if (acc[(size_t)k].count > big_c) { big_c=acc[(size_t)k].count; big_k=k; }

    for (int k = 0; k < K; ++k) {
      if (acc[(size_t)k].count == 0) {
        // steal: perturb the biggest cluster center slightly
        Site s = sites[(size_t)big_k];
        s.x += 0.5f; s.y -= 0.5f; // tiny shift
        sites[(size_t)k] = s;
        continue;
      }
      double inv = 1.0 / (double)acc[(size_t)k].count;
      double x = acc[(size_t)k].sx * inv;
      double y = acc[(size_t)k].sy * inv;
      double z = acc[(size_t)k].sz * inv;
      uint64_t ha = 0ull, hb = 0ull;
      for (int b=0;b<64;b++) if (acc[(size_t)k].bitcA[b] * 2 >= acc[(size_t)k].count) ha |= (1ull<<b);
      if (bits==128) for (int b=0;b<64;b++) if (acc[(size_t)k].bitcB[b] * 2 >= acc[(size_t)k].count) hb |= (1ull<<b);
      sites[(size_t)k] = Site{ (float)x, (float)y, (float)z, ha, hb };
    }

    // reinitialize carry to current labels so next flood starts from assigned seeds
    for (size_t i=0;i<Ngrid;++i) carry[i] = label[i];
  }

  // Pull labels for mask voxels (1-based cluster ids to be friendly in R)
  IntegerVector lab_mask(Nmask);
  for (int i = 0; i < Nmask; ++i) {
    int idx = mask_lin[i];
    int lab = label[(size_t)idx];
    lab_mask[i] = (lab >= 0) ? (lab + 1) : NA_INTEGER;
  }
  return lab_mask;
}