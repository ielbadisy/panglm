#include "families.h"
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

namespace {

// Parallel group-sum reduction: for each of n_groups groups (indexed by a
// 0-based code, one per row), accumulate the column sums of M and the
// group counts. Standard RcppParallel reduce pattern - each split keeps
// its own local accumulator, combined via join().
struct GroupSumWorker : public Worker {
  const RMatrix<double> M;
  const RVector<int> code;
  std::size_t n_groups, p;
  std::vector<double> sums;
  std::vector<int> counts;

  GroupSumWorker(const NumericMatrix& M, const IntegerVector& code, std::size_t n_groups)
    : M(M), code(code), n_groups(n_groups), p(M.ncol()),
      sums(n_groups * M.ncol(), 0.0), counts(n_groups, 0) {}

  GroupSumWorker(const GroupSumWorker& other, Split)
    : M(other.M), code(other.code), n_groups(other.n_groups), p(other.p),
      sums(other.n_groups * other.p, 0.0), counts(other.n_groups, 0) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      int g = code[i];
      counts[g]++;
      for (std::size_t j = 0; j < p; ++j) sums[g * p + j] += M(i, j);
    }
  }

  void join(const GroupSumWorker& rhs) {
    for (std::size_t g = 0; g < n_groups; ++g) {
      counts[g] += rhs.counts[g];
      for (std::size_t j = 0; j < p; ++j) sums[g * p + j] += rhs.sums[g * p + j];
    }
  }
};

// Parallel elementwise demeaning: row i -= group_mean(code[i], .). The
// natural per-observation parallelism once group means are known.
struct DemeanApplyWorker : public Worker {
  RMatrix<double> M;
  const RVector<int> code;
  const std::vector<double>& means; // n_groups x p, row-major
  std::size_t p;

  DemeanApplyWorker(NumericMatrix& M, const IntegerVector& code, const std::vector<double>& means)
    : M(M), code(code), means(means), p(M.ncol()) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      int g = code[i];
      for (std::size_t j = 0; j < p; ++j) M(i, j) -= means[g * p + j];
    }
  }
};

void demean_by(NumericMatrix& M, const IntegerVector& code, int n_groups) {
  GroupSumWorker worker(M, code, n_groups);
  parallelReduce(0, M.nrow(), worker);

  std::vector<double> means(worker.sums.size());
  std::size_t p = M.ncol();
  for (int g = 0; g < n_groups; ++g) {
    double cnt = std::max(worker.counts[g], 1);
    for (std::size_t j = 0; j < p; ++j) means[g * p + j] = worker.sums[g * p + j] / cnt;
  }

  DemeanApplyWorker applier(M, code, means);
  parallelFor(0, M.nrow(), applier);
}

} // namespace

// Two-way (individual + time) demeaning via alternating projections
// (Gauss-Seidel), parallelized across observations for both the group-sum
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
  for (iter = 0; iter < maxit; ++iter) {
    std::copy(cur.begin(), cur.end(), prev.begin());
    demean_by(cur, id_code, n_id);
    demean_by(cur, time_code, n_time);

    double maxdiff = 0.0;
    for (int i = 0; i < n * p; ++i) maxdiff = std::max(maxdiff, std::fabs(cur[i] - prev[i]));
    if (maxdiff < tol) { iter++; break; }
  }

  return List::create(Named("M") = cur, Named("iterations") = iter);
}
