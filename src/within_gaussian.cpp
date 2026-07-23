#include "families.h"
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// Parallel per-group demeaning: for each panel unit i, subtract the
// group mean of y and each column of X. Groups are contiguous row ranges
// [group_start[g], group_start[g] + group_size[g]) - the caller (R side,
// via data.table) must have already sorted the data by (id, time).
struct DemeanWorker : public Worker {
  const RMatrix<double> X;
  const RVector<double> y;
  const RVector<int> group_start;
  const RVector<int> group_size;
  RMatrix<double> Xd;
  RVector<double> yd;

  DemeanWorker(const NumericMatrix& X, const NumericVector& y,
               const IntegerVector& group_start, const IntegerVector& group_size,
               NumericMatrix& Xd, NumericVector& yd)
    : X(X), y(y), group_start(group_start), group_size(group_size), Xd(Xd), yd(yd) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::size_t ncol = X.ncol();
    for (std::size_t g = begin; g < end; ++g) {
      int start = group_start[g];
      int size  = group_size[g];
      if (size <= 0) continue;

      double ymean = 0.0;
      for (int t = 0; t < size; ++t) ymean += y[start + t];
      ymean /= size;
      for (int t = 0; t < size; ++t) yd[start + t] = y[start + t] - ymean;

      for (std::size_t j = 0; j < ncol; ++j) {
        double xmean = 0.0;
        for (int t = 0; t < size; ++t) xmean += X(start + t, j);
        xmean /= size;
        for (int t = 0; t < size; ++t) Xd(start + t, j) = X(start + t, j) - xmean;
      }
    }
  }
};

// [[Rcpp::export]]
List within_demean_cpp(NumericMatrix X, NumericVector y,
                        IntegerVector group_start, IntegerVector group_size) {
  int n = X.nrow(), k = X.ncol();
  NumericMatrix Xd(n, k);
  NumericVector yd(n);

  DemeanWorker worker(X, y, group_start, group_size, Xd, yd);
  parallelFor(0, group_start.size(), worker);

  return List::create(Named("X") = Xd, Named("y") = yd);
}
