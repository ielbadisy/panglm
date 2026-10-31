#include <Rcpp.h>
using namespace Rcpp;

// Sums of v by 0-based integer code (length n_codes), adding the
// observations in order, as rowsum() does.
//
// [[Rcpp::export]]
NumericVector code_sum_cpp(const NumericVector& v, const IntegerVector& code, int n_codes) {
  NumericVector out(n_codes);
  const int n = v.size();
  for (int i = 0; i < n; ++i) out[code[i]] += v[i];
  return out;
}
