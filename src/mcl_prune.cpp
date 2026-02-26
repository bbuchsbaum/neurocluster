#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List mcl_prune_sparse_cpp(const IntegerVector& p,
                                const IntegerVector& i,
                                const NumericVector& x,
                                const int ncol,
                                const int max_per_col,
                                const double min_value) {
  if (ncol < 0) {
    stop("mcl_prune_sparse_cpp: ncol must be non-negative");
  }
  if (max_per_col < 1) {
    stop("mcl_prune_sparse_cpp: max_per_col must be >= 1");
  }
  if (p.size() != (ncol + 1)) {
    stop("mcl_prune_sparse_cpp: invalid column pointer length");
  }
  if (i.size() != x.size()) {
    stop("mcl_prune_sparse_cpp: row index and value lengths differ");
  }

  const int nnz_in = x.size();
  std::vector<int> out_p(ncol + 1, 0);
  std::vector<int> out_i;
  std::vector<double> out_x;
  out_i.reserve(nnz_in);
  out_x.reserve(nnz_in);

  const bool use_threshold = R_finite(min_value) && (min_value > 0.0);
  std::vector<int> keep;

  for (int col = 0; col < ncol; ++col) {
    const int start = p[col];
    const int end = p[col + 1];

    if (start < 0 || end < start || end > nnz_in) {
      stop("mcl_prune_sparse_cpp: malformed sparse column pointers");
    }

    keep.clear();
    const int width = end - start;
    if (width > static_cast<int>(keep.capacity())) {
      keep.reserve(width);
    }

    if (!use_threshold) {
      keep.resize(width);
      std::iota(keep.begin(), keep.end(), start);
    } else {
      for (int idx = start; idx < end; ++idx) {
        if (x[idx] > min_value) {
          keep.push_back(idx);
        }
      }
    }

    if (keep.size() > static_cast<size_t>(max_per_col)) {
      auto by_value_desc = [&](const int a, const int b) {
        const double va = x[a];
        const double vb = x[b];
        if (va == vb) {
          return i[a] < i[b];
        }
        return va > vb;
      };

      std::nth_element(
        keep.begin(),
        keep.begin() + max_per_col,
        keep.end(),
        by_value_desc
      );
      keep.resize(max_per_col);
    }

    std::sort(
      keep.begin(),
      keep.end(),
      [&](const int a, const int b) {
        const int ra = i[a];
        const int rb = i[b];
        if (ra == rb) return a < b;
        return ra < rb;
      }
    );

    for (int idx : keep) {
      out_i.push_back(i[idx]);
      out_x.push_back(x[idx]);
    }

    out_p[col + 1] = static_cast<int>(out_i.size());
  }

  return List::create(
    _["p"] = wrap(out_p),
    _["i"] = wrap(out_i),
    _["x"] = wrap(out_x)
  );
}
