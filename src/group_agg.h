#ifndef PANGLM_GROUP_AGG_H
#define PANGLM_GROUP_AGG_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>

namespace panglm {

// Per-group aggregation of a Poisson linear predictor: Li[g] = sum(lit)
// within group g, lXi(g, :) = colSum(lit * X) within group g. Groups are
// contiguous row ranges [group_start[g], group_start[g] + group_size[g]).
// Parallelized over groups (panel units).
struct GroupAggWorker : public RcppParallel::Worker {
  const RcppParallel::RMatrix<double> X;
  const RcppParallel::RVector<double> lit;
  const RcppParallel::RVector<int> group_start;
  const RcppParallel::RVector<int> group_size;
  RcppParallel::RVector<double> Li;
  RcppParallel::RMatrix<double> lXi;

  GroupAggWorker(const Rcpp::NumericMatrix& X, const Rcpp::NumericVector& lit,
                 const Rcpp::IntegerVector& group_start, const Rcpp::IntegerVector& group_size,
                 Rcpp::NumericVector& Li, Rcpp::NumericMatrix& lXi)
    : X(X), lit(lit), group_start(group_start), group_size(group_size), Li(Li), lXi(lXi) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::size_t ncol = X.ncol();
    for (std::size_t g = begin; g < end; ++g) {
      int start = group_start[g];
      int size  = group_size[g];
      double s = 0.0;
      for (int t = 0; t < size; ++t) s += lit[start + t];
      Li[g] = s;
      for (std::size_t j = 0; j < ncol; ++j) {
        double sj = 0.0;
        for (int t = 0; t < size; ++t) sj += lit[start + t] * X(start + t, j);
        lXi(g, j) = sj;
      }
    }
  }
};

inline void compute_group_agg(const arma::mat& X, const arma::vec& lit,
                               const Rcpp::IntegerVector& group_start,
                               const Rcpp::IntegerVector& group_size,
                               arma::vec& Li, arma::mat& lXi) {
  int n = X.n_rows, k = X.n_cols, G = group_start.size();
  Rcpp::NumericMatrix Xr(n, k);
  for (int i = 0; i < n; ++i) for (int j = 0; j < k; ++j) Xr(i, j) = X(i, j);
  Rcpp::NumericVector litr(lit.begin(), lit.end());
  Rcpp::NumericVector Li_r(G);
  Rcpp::NumericMatrix lXi_r(G, k);

  GroupAggWorker worker(Xr, litr, group_start, group_size, Li_r, lXi_r);
  RcppParallel::parallelFor(0, G, worker);

  Li.set_size(G);
  lXi.set_size(G, k);
  for (int g = 0; g < G; ++g) {
    Li[g] = Li_r[g];
    for (int j = 0; j < k; ++j) lXi(g, j) = lXi_r(g, j);
  }
}

} // namespace panglm

#endif
