#ifndef PANGLM_THREADS_H
#define PANGLM_THREADS_H

#include <Rcpp.h>
#include <cstdlib>
#include <cstring>
#include <algorithm>

// Thread control for the RcppParallel loops. The package builds RcppParallel
// with its tinythread backend, which ignores the numThreads argument of
// parallelFor()/parallelReduce() and starts one thread per chunk. The number
// of threads is therefore capped through the grain size: a range of length n
// split with grain ceil(n / k) gives at most k chunks.
namespace panglm {

// Maximum number of threads, or 0 for no cap (RcppParallel's default).
// Sources, combined by taking the smallest positive value:
//   options(panglm.threads = k)
//   OMP_THREAD_LIMIT (set to 2 by CRAN's check farm)
//   _R_CHECK_LIMIT_CORES_ (set by R CMD check --as-cran): cap at 2
inline int num_threads() {
  int k = 0;
  SEXP opt = Rf_GetOption1(Rf_install("panglm.threads"));
  if (!Rf_isNull(opt) && Rf_length(opt) >= 1 &&
      (Rf_isInteger(opt) || Rf_isReal(opt))) {
    double v = Rf_asReal(opt);
    if (R_FINITE(v) && v >= 1) k = (int) v;
  }
  const char* lim = std::getenv("OMP_THREAD_LIMIT");
  if (lim != NULL) {
    int v = std::atoi(lim);
    if (v > 0) k = (k > 0) ? std::min(k, v) : v;
  }
  const char* chk = std::getenv("_R_CHECK_LIMIT_CORES_");
  if (chk != NULL && std::strlen(chk) > 0 &&
      std::strcmp(chk, "false") != 0 && std::strcmp(chk, "FALSE") != 0) {
    k = (k > 0) ? std::min(k, 2) : 2;
  }
  return k;
}

// Grain size passed to parallelFor()/parallelReduce() so that a loop over
// n items uses at most num_threads() threads.
inline std::size_t grain(std::size_t n) {
  int k = num_threads();
  if (k <= 0 || n == 0) return 1;
  return std::max<std::size_t>(1, (n + (std::size_t) k - 1) / (std::size_t) k);
}

} // namespace panglm

#endif
