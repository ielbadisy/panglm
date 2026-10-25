#' panglm: Generalized Linear Models for Panel Data
#'
#' @section Threads:
#' The within and random-effects estimators for Gaussian, binomial and
#' Poisson models compute per-group quantities in parallel with
#' RcppParallel. By default they use all available cores. Set
#' `options(panglm.threads = k)` to use at most `k` threads; the
#' environment variable `RCPP_PARALLEL_NUM_THREADS` is also honoured. When
#' `OMP_THREAD_LIMIT` is set, it caps the number of threads, and under
#' `R CMD check` (when `_R_CHECK_LIMIT_CORES_` is set) at most 2 threads
#' are used. Matrix algebra done in R goes through the BLAS library R is
#' linked to, whose threads are controlled outside panglm (for example with
#' `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` or `MKL_NUM_THREADS`).
#'
#' @keywords internal
#' @useDynLib panglm, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats model.frame model.matrix model.response pnorm printCoefmat nobs
"_PACKAGE"
