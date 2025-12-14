// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Worker for normalizing each volume (removing mean offset per timepoint)
struct NormalizeVolumesWorker : public Worker {
  // Input matrix: voxels x timepoints
  const RMatrix<double> input;
  // Output matrix: voxels x timepoints
  RMatrix<double> output;
  // Volume means (one per timepoint)
  const RVector<double> vol_means;

  NormalizeVolumesWorker(const NumericMatrix& input_,
                         NumericMatrix& output_,
                         const NumericVector& vol_means_)
    : input(input_), output(output_), vol_means(vol_means_) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::size_t n_time = input.ncol();
    for (std::size_t i = begin; i < end; i++) {
      for (std::size_t t = 0; t < n_time; t++) {
        output(i, t) = input(i, t) - vol_means[t];
      }
    }
  }
};

//' Normalize Volumes by Removing Mean Offset (Fast C++ Implementation)
//'
//' Removes the mean offset from each volume/timepoint in a matrix.
//' This centers each volume to have zero mean across voxels.
//'
//' @param data Numeric matrix with voxels as rows and timepoints as columns
//' @return Numeric matrix with each column (volume) centered to mean zero
//'
//' @details
//' For each timepoint t, computes mean_t = mean(data[,t]) and subtracts
//' it from all voxels: output[,t] = data[,t] - mean_t
//'
//' Uses RcppParallel for fast parallel computation across voxels.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix normalize_volumes_cpp(NumericMatrix data) {
  int n_voxels = data.nrow();
  int n_time = data.ncol();

  // Compute mean for each volume/timepoint
  NumericVector vol_means(n_time);
  for (int t = 0; t < n_time; t++) {
    double sum = 0.0;
    for (int i = 0; i < n_voxels; i++) {
      sum += data(i, t);
    }
    vol_means[t] = sum / n_voxels;
  }

  // Create output matrix
  NumericMatrix output(n_voxels, n_time);

  // Parallel subtraction
  NormalizeVolumesWorker worker(data, output, vol_means);
  parallelFor(0, n_voxels, worker);

  return output;
}


// Worker for detrending each voxel's timeseries
struct DetrendWorker : public Worker {
  // Input matrix: voxels x timepoints
  const RMatrix<double> input;
  // Output matrix: voxels x timepoints
  RMatrix<double> output;
  // Precomputed values for linear regression
  const double t_mean;
  const double t_var;
  const int n_time;

  DetrendWorker(const NumericMatrix& input_,
                NumericMatrix& output_,
                double t_mean_,
                double t_var_,
                int n_time_)
    : input(input_), output(output_),
      t_mean(t_mean_), t_var(t_var_), n_time(n_time_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      // Compute mean of this voxel's timeseries
      double y_mean = 0.0;
      for (int t = 0; t < n_time; t++) {
        y_mean += input(i, t);
      }
      y_mean /= n_time;

      // Compute slope: sum((t - t_mean) * (y - y_mean)) / sum((t - t_mean)^2)
      double cov_ty = 0.0;
      for (int t = 0; t < n_time; t++) {
        cov_ty += (t - t_mean) * (input(i, t) - y_mean);
      }
      double slope = cov_ty / t_var;
      double intercept = y_mean - slope * t_mean;

      // Remove linear trend
      for (int t = 0; t < n_time; t++) {
        output(i, t) = input(i, t) - (intercept + slope * t);
      }
    }
  }
};

//' Detrend Voxel Timeseries (Fast C++ Implementation)
//'
//' Removes linear trend from each voxel's timeseries using fast
//' parallel linear regression.
//'
//' @param data Numeric matrix with voxels as rows and timepoints as columns
//' @return Numeric matrix with linear trend removed from each row (voxel)
//'
//' @details
//' For each voxel, fits a linear model y = intercept + slope * t and
//' subtracts the fitted values: residuals = y - fitted.
//'
//' Uses RcppParallel for fast parallel computation across voxels.
//' The linear regression is computed efficiently using precomputed
//' time statistics.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix detrend_time_cpp(NumericMatrix data) {
  int n_voxels = data.nrow();
  int n_time = data.ncol();

  // Precompute time statistics (0-indexed time)
  // t_mean = (n-1)/2, t_var = sum((t - t_mean)^2)
  double t_mean = (n_time - 1) / 2.0;
  double t_var = 0.0;
  for (int t = 0; t < n_time; t++) {
    double diff = t - t_mean;
    t_var += diff * diff;
  }

  // Create output matrix
  NumericMatrix output(n_voxels, n_time);

  // Parallel detrending
  DetrendWorker worker(data, output, t_mean, t_var, n_time);
  parallelFor(0, n_voxels, worker);

  return output;
}


// Worker for combined normalize + detrend
struct NormalizeDetrendWorker : public Worker {
  const RMatrix<double> input;
  RMatrix<double> output;
  const RVector<double> vol_means;
  const double t_mean;
  const double t_var;
  const int n_time;

  NormalizeDetrendWorker(const NumericMatrix& input_,
                         NumericMatrix& output_,
                         const NumericVector& vol_means_,
                         double t_mean_,
                         double t_var_,
                         int n_time_)
    : input(input_), output(output_), vol_means(vol_means_),
      t_mean(t_mean_), t_var(t_var_), n_time(n_time_) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> centered(n_time);

    for (std::size_t i = begin; i < end; i++) {
      // First: center by volume means
      double y_mean = 0.0;
      for (int t = 0; t < n_time; t++) {
        centered[t] = input(i, t) - vol_means[t];
        y_mean += centered[t];
      }
      y_mean /= n_time;

      // Second: compute and remove linear trend
      double cov_ty = 0.0;
      for (int t = 0; t < n_time; t++) {
        cov_ty += (t - t_mean) * (centered[t] - y_mean);
      }
      double slope = cov_ty / t_var;
      double intercept = y_mean - slope * t_mean;

      for (int t = 0; t < n_time; t++) {
        output(i, t) = centered[t] - (intercept + slope * t);
      }
    }
  }
};

//' Normalize Volumes and Detrend Timeseries (Fast C++ Implementation)
//'
//' Combined operation: first removes volume mean offsets, then removes
//' linear trend from each voxel's timeseries. More efficient than
//' calling both functions separately.
//'
//' @param data Numeric matrix with voxels as rows and timepoints as columns
//' @return Numeric matrix with volume offsets and linear trends removed
//'
//' @details
//' Performs two operations in a single parallel pass:
//' 1. Centers each volume to zero mean
//' 2. Removes linear trend from each voxel's (now centered) timeseries
//'
//' @export
// [[Rcpp::export]]
NumericMatrix normalize_detrend_cpp(NumericMatrix data) {
  int n_voxels = data.nrow();
  int n_time = data.ncol();

  // Compute volume means
  NumericVector vol_means(n_time);
  for (int t = 0; t < n_time; t++) {
    double sum = 0.0;
    for (int i = 0; i < n_voxels; i++) {
      sum += data(i, t);
    }
    vol_means[t] = sum / n_voxels;
  }

  // Precompute time statistics
  double t_mean = (n_time - 1) / 2.0;
  double t_var = 0.0;
  for (int t = 0; t < n_time; t++) {
    double diff = t - t_mean;
    t_var += diff * diff;
  }

  // Create output matrix
  NumericMatrix output(n_voxels, n_time);

  // Combined parallel operation
  NormalizeDetrendWorker worker(data, output, vol_means, t_mean, t_var, n_time);
  parallelFor(0, n_voxels, worker);

  return output;
}


// =============================================================================
// Flexible Detrending with Polynomial or DCT Basis
// =============================================================================

//' Generate DCT-II Basis Matrix
//'
//' Creates a discrete cosine transform (type II) basis matrix for detrending.
//'
//' @param n_time Number of timepoints
//' @param n_basis Number of DCT basis functions (including constant)
//' @return Numeric matrix (n_time x n_basis) with orthonormal DCT basis
//'
//' @details
//' DCT basis functions are excellent for removing low-frequency drift in fMRI.
//' The first basis is constant (mean), subsequent bases capture increasingly
//' higher frequency components.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix make_dct_basis(int n_time, int n_basis) {
  NumericMatrix basis(n_time, n_basis);

  for (int k = 0; k < n_basis; k++) {
    double norm = (k == 0) ? std::sqrt(1.0 / n_time) : std::sqrt(2.0 / n_time);
    for (int t = 0; t < n_time; t++) {
      basis(t, k) = norm * std::cos(M_PI * k * (t + 0.5) / n_time);
    }
  }

  return basis;
}


//' Generate Polynomial Basis Matrix
//'
//' Creates a polynomial basis matrix for detrending.
//'
//' @param n_time Number of timepoints
//' @param degree Polynomial degree (0 = constant, 1 = linear, 2 = quadratic, etc.)
//' @return Numeric matrix (n_time x (degree+1)) with orthonormalized polynomial basis
//'
//' @details
//' Uses Gram-Schmidt orthonormalization for numerical stability.
//' Degree 0 = mean removal, 1 = linear detrend, 2 = quadratic, etc.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix make_poly_basis(int n_time, int degree) {
  int n_basis = degree + 1;
  NumericMatrix basis(n_time, n_basis);

  // Normalized time from -1 to 1
  std::vector<double> t_norm(n_time);
  for (int t = 0; t < n_time; t++) {
    t_norm[t] = 2.0 * t / (n_time - 1) - 1.0;
  }

  // Build polynomial basis with Gram-Schmidt orthonormalization
  for (int d = 0; d <= degree; d++) {
    // Start with raw polynomial
    for (int t = 0; t < n_time; t++) {
      basis(t, d) = std::pow(t_norm[t], d);
    }

    // Orthogonalize against previous basis vectors
    for (int j = 0; j < d; j++) {
      double dot_prod = 0.0;
      for (int t = 0; t < n_time; t++) {
        dot_prod += basis(t, d) * basis(t, j);
      }
      for (int t = 0; t < n_time; t++) {
        basis(t, d) -= dot_prod * basis(t, j);
      }
    }

    // Normalize
    double norm = 0.0;
    for (int t = 0; t < n_time; t++) {
      norm += basis(t, d) * basis(t, d);
    }
    norm = std::sqrt(norm);
    if (norm > 1e-10) {
      for (int t = 0; t < n_time; t++) {
        basis(t, d) /= norm;
      }
    }
  }

  return basis;
}


// Worker for basis projection detrending
struct BasisDetrendWorker : public Worker {
  const RMatrix<double> input;
  RMatrix<double> output;
  const RMatrix<double> basis;      // n_time x n_basis
  const RMatrix<double> proj_mat;   // Precomputed I - B * B^T (n_time x n_time)
  const int n_time;
  const bool use_proj_mat;

  BasisDetrendWorker(const NumericMatrix& input_,
                     NumericMatrix& output_,
                     const NumericMatrix& basis_,
                     const NumericMatrix& proj_mat_,
                     int n_time_,
                     bool use_proj_mat_)
    : input(input_), output(output_), basis(basis_), proj_mat(proj_mat_),
      n_time(n_time_), use_proj_mat(use_proj_mat_) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n_basis = basis.ncol();
    std::vector<double> coeffs(n_basis);

    for (std::size_t i = begin; i < end; i++) {
      if (use_proj_mat) {
        // Use precomputed projection matrix: y_detrend = (I - B*B^T) * y
        for (int t = 0; t < n_time; t++) {
          double val = 0.0;
          for (int s = 0; s < n_time; s++) {
            val += proj_mat(t, s) * input(i, s);
          }
          output(i, t) = val;
        }
      } else {
        // Compute coefficients: beta = B^T * y (since B is orthonormal)
        for (int k = 0; k < n_basis; k++) {
          double coef = 0.0;
          for (int t = 0; t < n_time; t++) {
            coef += basis(t, k) * input(i, t);
          }
          coeffs[k] = coef;
        }

        // Remove fitted values: y - B * beta
        for (int t = 0; t < n_time; t++) {
          double fitted = 0.0;
          for (int k = 0; k < n_basis; k++) {
            fitted += basis(t, k) * coeffs[k];
          }
          output(i, t) = input(i, t) - fitted;
        }
      }
    }
  }
};


//' Flexible Detrending with Basis Functions (Fast C++ Implementation)
//'
//' Removes trends from each voxel's timeseries by projecting out a basis.
//' Supports polynomial or DCT (discrete cosine transform) basis functions.
//'
//' @param data Numeric matrix with voxels as rows and timepoints as columns
//' @param basis Orthonormal basis matrix (n_time x n_basis) to project out.
//'   Use make_dct_basis() or make_poly_basis() to generate.
//' @return Numeric matrix with basis components removed from each row
//'
//' @details
//' For orthonormal basis B, computes residuals as:
//' y_detrend = y - B * (B^T * y)
//'
//' This is equivalent to regressing out the basis and keeping residuals.
//'
//' @examples
//' \dontrun
//' # Remove linear + quadratic trend (polynomial degree 2)
//' basis <- make_poly_basis(n_time = 100, degree = 2)
//' data_detrend <- detrend_basis_cpp(data, basis)
//'
//' # Remove low-frequency drift with DCT (first 5 components)
//' basis <- make_dct_basis(n_time = 100, n_basis = 5)
//' data_detrend <- detrend_basis_cpp(data, basis)
//' }
//'
//' @export
// [[Rcpp::export]]
NumericMatrix detrend_basis_cpp(NumericMatrix data, NumericMatrix basis) {
  int n_voxels = data.nrow();
  int n_time = data.ncol();
  int n_basis = basis.ncol();

  if (basis.nrow() != n_time) {
    stop("Basis matrix rows must equal number of timepoints");
  }

  // Create output matrix
  NumericMatrix output(n_voxels, n_time);

  // For small n_basis relative to n_time, direct computation is faster
  // For larger n_basis, precompute projection matrix
  bool use_proj_mat = (n_basis > n_time / 4);

  NumericMatrix proj_mat(1, 1);  // Dummy if not used

  if (use_proj_mat) {
    // Precompute projection matrix: I - B * B^T
    proj_mat = NumericMatrix(n_time, n_time);
    for (int i = 0; i < n_time; i++) {
      for (int j = 0; j < n_time; j++) {
        double bb = 0.0;
        for (int k = 0; k < n_basis; k++) {
          bb += basis(i, k) * basis(j, k);
        }
        proj_mat(i, j) = (i == j ? 1.0 : 0.0) - bb;
      }
    }
  }

  // Parallel detrending
  BasisDetrendWorker worker(data, output, basis, proj_mat, n_time, use_proj_mat);
  parallelFor(0, n_voxels, worker);

  return output;
}


//' Polynomial Detrending (Convenience Function)
//'
//' Removes polynomial trend up to specified degree from each voxel.
//'
//' @param data Numeric matrix with voxels as rows and timepoints as columns
//' @param degree Polynomial degree (0 = mean, 1 = linear, 2 = quadratic, etc.)
//' @return Numeric matrix with polynomial trend removed
//'
//' @details
//' Convenience wrapper that generates polynomial basis and calls detrend_basis_cpp.
//' For degree=1, this is equivalent to detrend_time_cpp but slightly slower.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix detrend_poly_cpp(NumericMatrix data, int degree = 1) {
  int n_time = data.ncol();
  NumericMatrix basis = make_poly_basis(n_time, degree);
  return detrend_basis_cpp(data, basis);
}


//' DCT Detrending (Convenience Function)
//'
//' Removes low-frequency components using discrete cosine transform.
//'
//' @param data Numeric matrix with voxels as rows and timepoints as columns
//' @param n_basis Number of DCT basis functions to remove (including constant).
//'   Default 4 removes constant, linear, and first two cosine harmonics.
//' @return Numeric matrix with low-frequency components removed
//'
//' @details
//' DCT detrending is commonly used in fMRI preprocessing. The number of
//' basis functions controls the high-pass filter cutoff:
//' - n_basis = 1: mean removal only
//' - n_basis = 2: ~ linear detrend
//' - n_basis = 4-6: typical for fMRI (removes drift < ~0.01 Hz for TR=2s)
//'
//' Higher n_basis = more aggressive high-pass filtering.
//'
//' @export
// [[Rcpp::export]]
NumericMatrix detrend_dct_cpp(NumericMatrix data, int n_basis = 4) {
  int n_time = data.ncol();
  if (n_basis > n_time) {
    n_basis = n_time;
  }
  NumericMatrix basis = make_dct_basis(n_time, n_basis);
  return detrend_basis_cpp(data, basis);
}
