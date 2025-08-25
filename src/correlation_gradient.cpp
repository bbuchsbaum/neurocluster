#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

int index_3d_to_1d(int i, int j, int k, const IntegerVector& dims) {
  return i + j * dims[0] + k * dims[0] * dims[1];
}

// [[Rcpp::export]]
NumericVector correlation_gradient_cpp(const NumericVector& img_4d, const NumericVector& brain_mask) {
  // Safely get dimensions with error checking
  if (!img_4d.hasAttribute("dim")) {
    stop("img_4d must have dim attribute");
  }
  if (!brain_mask.hasAttribute("dim")) {
    stop("brain_mask must have dim attribute");
  }
  
  IntegerVector img_4d_dims = img_4d.attr("dim");
  IntegerVector brain_mask_dims = brain_mask.attr("dim");

  NumericVector gradient_3d(brain_mask_dims[0] * brain_mask_dims[1] * brain_mask_dims[2]);

  for (int i = 1; i < brain_mask_dims[0] - 1; ++i) {
    for (int j = 1; j < brain_mask_dims[1] - 1; ++j) {
      for (int k = 1; k < brain_mask_dims[2] - 1; ++k) {
        if (brain_mask(index_3d_to_1d(i, j, k, brain_mask_dims)) == 1) {
          double local_correlations[3][3][3];
          double local_gradient = 0.0;

          for (int i_n = 0; i_n < 3; ++i_n) {
            for (int j_n = 0; j_n < 3; ++j_n) {
              for (int k_n = 0; k_n < 3; ++k_n) {
                if (brain_mask(index_3d_to_1d(i - 1 + i_n, j - 1 + j_n, k - 1 + k_n, brain_mask_dims)) == 1) {
                  double correlation = 0.0;
                  double count = 0.0;

                  for (int t = 0; t < img_4d_dims[3]; ++t) {
                    correlation += img_4d(i - 1 + i_n + (j - 1 + j_n) * img_4d_dims[0] + (k - 1 + k_n) * img_4d_dims[0] * img_4d_dims[1] + t * img_4d_dims[0] * img_4d_dims[1] * img_4d_dims[2]) *
                      img_4d(i + j * img_4d_dims[0] + k * img_4d_dims[0] * img_4d_dims[1] + t * img_4d_dims[0] * img_4d_dims[1] * img_4d_dims[2]);
                    count += 1.0;
                  }

                  local_correlations[i_n][j_n][k_n] = correlation / count;
                } else {
                  local_correlations[i_n][j_n][k_n] = 0;
                }
              }
            }
          }

          for (int i_n = 0; i_n < 3; ++i_n) {
            for (int j_n = 0; j_n < 3; ++j_n) {
              for (int k_n = 0; k_n < 3; ++k_n) {
                if (brain_mask(index_3d_to_1d(i - 1 + i_n, j - 1 + j_n, k - 1 + k_n, brain_mask_dims)) == 1) {
                  local_gradient += pow(local_correlations[i_n][j_n][k_n] - local_correlations[1][1][1], 2);
                }
              }
            }
          }

          gradient_3d(index_3d_to_1d(i, j, k, brain_mask_dims)) = sqrt(local_gradient);
        }
      }
    }
  }

  gradient_3d.attr("dim") = brain_mask_dims;
  return gradient_3d;
}

RCPP_MODULE(correlation_gradient_module) {
  Rcpp::function("correlation_gradient_cpp", &correlation_gradient_cpp);
}

