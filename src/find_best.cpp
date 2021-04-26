// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double heat_kernel(NumericVector x1, NumericVector x2, double sigma) {
  double dist = sqrt(sum(pow((x1-x2),2)));
  return exp(-dist / (2*(pow(sigma,2))));
}

// [[Rcpp::export]]
double normalized_heat_kernel(NumericVector x1, NumericVector x2, double sigma) {
  double norm_dist = sum(pow((x1-x2),2))/(2*x1.length());
  return exp(-norm_dist/(2*pow(sigma,2)));
}

// [[Rcpp::export]]
NumericVector compute_scores(IntegerVector curclus, NumericMatrix coords, NumericMatrix data_centroids,
                             NumericMatrix coord_centroids, NumericMatrix data, double sigma1, double sigma2) {
  int n = curclus.length();

  //List scores(n);
  NumericVector out(n);

  for (int i=0; i<n; i++) {
    double c1 = normalized_heat_kernel(data(_, i), data_centroids(_, curclus[i]-1), sigma1);
    double c2 = heat_kernel(coords(_, i), coord_centroids(_, curclus[i]-1), sigma2);
    out[i] = c1 + c2;
  }

  return out;
}


// [[Rcpp::export]]
IntegerVector best_candidate(List candidates, IntegerVector curclus, NumericMatrix coords, NumericMatrix data_centroids,
                             NumericMatrix coord_centroids, NumericMatrix data, double sigma1, double sigma2, double alpha) {
  int n = candidates.length();
  int nswitches = 0;

  //List scores(n);
  IntegerVector out(n);

  for (int i=0; i<n; i++) {

    IntegerVector cand = candidates[i];

    if (cand.length() <= 1) {
      out[i] = curclus[i];
    } else {
      NumericVector score(cand.length());

      for (int j=0; j<cand.length(); j++) {
        double c1 = normalized_heat_kernel(data(_, i), data_centroids(_, cand[j]-1), sigma1);
        double c2 = heat_kernel(coords(_, i), coord_centroids(_, cand[j]-1), sigma2);
        score[j] = alpha*c1 + (1-alpha)*c2;
      }

      //scores[i] = score;
      out[i] = cand[which_max(score)];

      if (out[i] != curclus[i]) {
        nswitches++;
      }
    }

  }

  out.attr("nswitches") = nswitches;
  //out.attr("scores") = scores;

  return out;

}


// [[Rcpp::export]]
List find_candidates(IntegerMatrix nn_index, NumericMatrix nn_dist, IntegerVector curclus, double dthresh) {
  int n = nn_index.nrow();
  List out(n);


  for(int i = 0; i < n; ++i) {
    NumericMatrix::Row D = nn_dist( i, _);
    IntegerMatrix::Row ind = nn_index(i, _);

    std::unordered_set<int> cand;
    cand.insert(curclus[i]);

    for (int j = 0; j < D.size(); ++j) {
      if (D[j] < dthresh) {
        cand.insert(curclus[ind[j]]);
      }
    }

    out[i] = cand;
  }

  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
library(FNN)
coords  <- cbind(1:100, 1:100, 1:100)
nn <- FNN::get.knn(coords, 5)
kres <- kmeans(coords, 25)
cclus <- kres$cluster
dthresh=8
ret <- find_candidates(nn$nn.index-1, nn$nn.dist, as.integer(cclus), dthresh)

x = heat_kernel(1:4, 4:1, .5)

*/
