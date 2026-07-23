#include "families.h"
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

namespace {

// Per-group means of X and y, and group sizes Ti - the building block for
// the Swamy-Arora variance-components estimator behind the Gaussian
// random-effects (quasi-within) panel model. Parallelized over groups.
struct GroupMeanWorker : public Worker {
  const RMatrix<double> X;
  const RVector<double> y;
  const RVector<int> group_start;
  const RVector<int> group_size;
  RMatrix<double> Xbar;
  RVector<double> ybar;

  GroupMeanWorker(const NumericMatrix& X, const NumericVector& y,
                   const IntegerVector& group_start, const IntegerVector& group_size,
                   NumericMatrix& Xbar, NumericVector& ybar)
    : X(X), y(y), group_start(group_start), group_size(group_size), Xbar(Xbar), ybar(ybar) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::size_t ncol = X.ncol();
    for (std::size_t g = begin; g < end; ++g) {
      int start = group_start[g];
      int size  = group_size[g];
      double ys = 0.0;
      for (int t = 0; t < size; ++t) ys += y[start + t];
      ybar[g] = ys / size;
      for (std::size_t j = 0; j < ncol; ++j) {
        double xs = 0.0;
        for (int t = 0; t < size; ++t) xs += X(start + t, j);
        Xbar(g, j) = xs / size;
      }
    }
  }
};

// Quasi-demeaning: row (i,t) becomes value_it - theta_i * mean_i, where
// theta_i is the group-specific GLS weight from the Swamy-Arora RE model.
struct QuasiDemeanWorker : public Worker {
  const RMatrix<double> X;
  const RVector<double> y;
  const RVector<int> group_start;
  const RVector<int> group_size;
  const RVector<double> theta;
  RMatrix<double> Xq;
  RVector<double> yq;

  QuasiDemeanWorker(const NumericMatrix& X, const NumericVector& y,
                     const IntegerVector& group_start, const IntegerVector& group_size,
                     const NumericVector& theta, NumericMatrix& Xq, NumericVector& yq)
    : X(X), y(y), group_start(group_start), group_size(group_size), theta(theta), Xq(Xq), yq(yq) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::size_t ncol = X.ncol();
    for (std::size_t g = begin; g < end; ++g) {
      int start = group_start[g];
      int size  = group_size[g];
      double th = theta[g];

      double ys = 0.0;
      for (int t = 0; t < size; ++t) ys += y[start + t];
      double ymean = ys / size;
      for (int t = 0; t < size; ++t) yq[start + t] = y[start + t] - th * ymean;

      for (std::size_t j = 0; j < ncol; ++j) {
        double xs = 0.0;
        for (int t = 0; t < size; ++t) xs += X(start + t, j);
        double xmean = xs / size;
        for (int t = 0; t < size; ++t) Xq(start + t, j) = X(start + t, j) - th * xmean;
      }
    }
  }
};

} // namespace

// [[Rcpp::export]]
List group_means_cpp(NumericMatrix X, NumericVector y,
                      IntegerVector group_start, IntegerVector group_size) {
  int G = group_start.size(), k = X.ncol();
  NumericMatrix Xbar(G, k);
  NumericVector ybar(G);
  GroupMeanWorker worker(X, y, group_start, group_size, Xbar, ybar);
  parallelFor(0, G, worker);
  return List::create(Named("Xbar") = Xbar, Named("ybar") = ybar);
}

// [[Rcpp::export]]
List quasi_demean_cpp(NumericMatrix X, NumericVector y,
                       IntegerVector group_start, IntegerVector group_size,
                       NumericVector theta) {
  int n = X.nrow(), k = X.ncol();
  NumericMatrix Xq(n, k);
  NumericVector yq(n);
  QuasiDemeanWorker worker(X, y, group_start, group_size, theta, Xq, yq);
  parallelFor(0, group_start.size(), worker);
  return List::create(Named("X") = Xq, Named("y") = yq);
}
