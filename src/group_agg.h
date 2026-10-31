#ifndef PANGLM_GROUP_AGG_H
#define PANGLM_GROUP_AGG_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "threads.h"

namespace panglm {

// Per-group aggregation of a Poisson linear predictor: Li[g] = sum(lit)
// within group g, lXi(g, :) = colSum(lit * X) within group g. Groups are
// contiguous row ranges [group_start[g], group_start[g] + group_size[g]).
// Parallelized over groups (panel units).
// The worker reads and writes plain column-major buffers, so the inputs and
// outputs are used in place without copies to and from R matrices.
struct GroupAggWorker : public RcppParallel::Worker {
  const double* X;
  const double* lit;
  const int* group_start;
  const int* group_size;
  double* Li;
  double* lXi;
  std::size_t n, G, ncol;

  GroupAggWorker(const double* X, const double* lit, const int* group_start,
                 const int* group_size, double* Li, double* lXi,
                 std::size_t n, std::size_t G, std::size_t ncol)
    : X(X), lit(lit), group_start(group_start), group_size(group_size), Li(Li),
      lXi(lXi), n(n), G(G), ncol(ncol) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t g = begin; g < end; ++g) {
      int start = group_start[g];
      int size  = group_size[g];
      double s = 0.0;
      for (int t = 0; t < size; ++t) s += lit[start + t];
      Li[g] = s;
      for (std::size_t j = 0; j < ncol; ++j) {
        const double* xj = X + j * n;
        double sj = 0.0;
        for (int t = 0; t < size; ++t) sj += lit[start + t] * xj[start + t];
        lXi[g + j * G] = sj;
      }
    }
  }
};

inline void compute_group_agg(const arma::mat& X, const arma::vec& lit,
                              const Rcpp::IntegerVector& group_start,
                              const Rcpp::IntegerVector& group_size,
                              arma::vec& Li, arma::mat& lXi) {
  std::size_t n = X.n_rows, k = X.n_cols, G = group_start.size();
  Li.set_size(G);
  lXi.set_size(G, k);
  GroupAggWorker worker(X.memptr(), lit.memptr(), group_start.begin(), group_size.begin(),
                        Li.memptr(), lXi.memptr(), n, G, k);
  RcppParallel::parallelFor(0, G, worker, panglm::grain(G));
}

} // namespace panglm

#endif
