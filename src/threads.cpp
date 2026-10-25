#include "threads.h"

// Thread cap used by the parallel loops (0 = no cap); for tests.
// [[Rcpp::export]]
int panglm_num_threads_cpp() {
  return panglm::num_threads();
}
