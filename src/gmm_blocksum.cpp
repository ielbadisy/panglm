#include <Rcpp.h>
using namespace Rcpp;

// Sums of each block of R consecutive rows of M (N x ncol result), adding
// the rows in order, as rowsum() does on a grouping with contiguous groups.
//
// [[Rcpp::export]]
NumericMatrix gmm_blocksum_cpp(const NumericMatrix& M, int R) {
  const int n = M.nrow(), p = M.ncol(), N = n / R;
  NumericMatrix out(N, p);
  const double* m = M.begin();
  double* o = out.begin();
  for (int j = 0; j < p; ++j) {
    const double* col = m + (size_t)j * n;
    double* oc = o + (size_t)j * N;
    for (int i = 0; i < N; ++i) {
      const double* b = col + (size_t)i * R;
      double s = 0.0;
      for (int r = 0; r < R; ++r) s += b[r];
      oc[i] = s;
    }
  }
  return out;
}
