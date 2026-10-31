#include "families.h"
#include <RcppParallel.h>
#include "threads.h"
using namespace Rcpp;
using namespace RcppParallel;

namespace {

// Group means are accumulated in one serial pass. A parallel reduction would
// give every thread its own n_groups x p accumulator and merge them, which
// costs far more than the pass itself when there are many groups (one per
// individual). The pass is O(n p) and memory bound; the demeaning step that
// follows runs in parallel across observations.

// Parallel elementwise demeaning: row i -= group_mean(code[i], .).
struct DemeanApplyWorker : public Worker {
  RMatrix<double> M;
  const RVector<int> code;
  const std::vector<double>& means; // n_groups x p, row-major
  std::size_t p;

  DemeanApplyWorker(NumericMatrix& M, const IntegerVector& code, const std::vector<double>& means)
    : M(M), code(code), means(means), p(M.ncol()) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = 0; j < p; ++j) {
      RMatrix<double>::Column col = M.column(j);
      for (std::size_t i = begin; i < end; ++i) col[i] -= means[(std::size_t) code[i] * p + j];
    }
  }
};

void demean_by(NumericMatrix& M, const IntegerVector& code, int n_groups) {
  const std::size_t n = M.nrow(), p = M.ncol();
  std::vector<double> means((std::size_t) n_groups * p, 0.0);
  std::vector<int> counts(n_groups, 0);
  const int* cd = code.begin();
  for (std::size_t i = 0; i < n; ++i) counts[cd[i]]++;
  for (std::size_t j = 0; j < p; ++j) {
    const double* col = &M[j * n];
    for (std::size_t i = 0; i < n; ++i) means[(std::size_t) cd[i] * p + j] += col[i];
  }
  for (int g = 0; g < n_groups; ++g) {
    double cnt = std::max(counts[g], 1);
    for (std::size_t j = 0; j < p; ++j) means[(std::size_t) g * p + j] /= cnt;
  }

  DemeanApplyWorker applier(M, code, means);
  parallelFor(0, n, applier, panglm::grain(n));
}

// Weighted analogue: group_mean = sum(w_i * M_i) / sum(w_i). Used by the
// IRLS-nested weighted two-way demeaning that non-gaussian twoways models
// (currently poisson) need; the IRLS weight changes every outer iteration,
// so this is re-run from scratch each time.
void demean_by_weighted(NumericMatrix& M, const NumericVector& w, const IntegerVector& code, int n_groups) {
  const std::size_t n = M.nrow(), p = M.ncol();
  std::vector<double> means((std::size_t) n_groups * p, 0.0);
  std::vector<double> wsums(n_groups, 0.0);
  const int* cd = code.begin();
  const double* wp = w.begin();
  for (std::size_t i = 0; i < n; ++i) wsums[cd[i]] += wp[i];
  for (std::size_t j = 0; j < p; ++j) {
    const double* col = &M[j * n];
    for (std::size_t i = 0; i < n; ++i) means[(std::size_t) cd[i] * p + j] += wp[i] * col[i];
  }
  for (int g = 0; g < n_groups; ++g) {
    double wg = std::max(wsums[g], 1e-12);
    for (std::size_t j = 0; j < p; ++j) means[(std::size_t) g * p + j] /= wg;
  }

  DemeanApplyWorker applier(M, code, means);
  parallelFor(0, n, applier, panglm::grain(n));
}

} // namespace

// Two-way (individual + time) demeaning via alternating projections
// (Gauss-Seidel), iterated until the largest change in one sweep is below
// `tol` or `maxit` sweeps are reached (`converged` reports which),
// parallelized across observations for both the group-sum
// reduction and the demeaning step. `id_code`/`time_code` are 0-based
// integer group codes, one per row (need not be contiguous/sorted).
//
// [[Rcpp::export]]
List twoway_demean_cpp(NumericMatrix M, IntegerVector id_code, IntegerVector time_code,
                       int n_id, int n_time, int maxit = 10000, double tol = 1e-10) {
  int n = M.nrow(), p = M.ncol();
  NumericMatrix cur = Rcpp::clone(M);
  NumericMatrix prev(n, p);

  int iter = 0;
  bool converged = false;
  for (iter = 0; iter < maxit; ++iter) {
    std::copy(cur.begin(), cur.end(), prev.begin());
    demean_by(cur, id_code, n_id);
    demean_by(cur, time_code, n_time);

    double maxdiff = 0.0;
    for (int i = 0; i < n * p; ++i) maxdiff = std::max(maxdiff, std::fabs(cur[i] - prev[i]));
    if (maxdiff < tol) { iter++; converged = true; break; }
  }

  return List::create(Named("M") = cur, Named("iterations") = iter,
                      Named("converged") = converged);
}

// Weighted two-way demeaning: the same alternating-projections algorithm,
// but each group mean is an IRLS-weighted mean. This is the inner step of
// the outer-IRLS / inner-FWL algorithm used for non-gaussian two-way FE
// (poisson, currently) -- see fit_within_twoways_poisson() in R.
//
// [[Rcpp::export]]
List twoway_demean_weighted_cpp(NumericMatrix M, NumericVector w,
                                IntegerVector id_code, IntegerVector time_code,
                                int n_id, int n_time, int maxit = 10000, double tol = 1e-10) {
  int n = M.nrow(), p = M.ncol();
  NumericMatrix cur = Rcpp::clone(M);
  NumericMatrix prev(n, p);

  int iter = 0;
  bool converged = false;
  for (iter = 0; iter < maxit; ++iter) {
    std::copy(cur.begin(), cur.end(), prev.begin());
    demean_by_weighted(cur, w, id_code, n_id);
    demean_by_weighted(cur, w, time_code, n_time);

    double maxdiff = 0.0;
    for (int i = 0; i < n * p; ++i) maxdiff = std::max(maxdiff, std::fabs(cur[i] - prev[i]));
    if (maxdiff < tol) { iter++; converged = true; break; }
  }

  return List::create(Named("M") = cur, Named("iterations") = iter,
                      Named("converged") = converged);
}
