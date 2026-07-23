#' panglm: Fast Generalized Linear Models for Panel Data
#'
#' @keywords internal
#' @useDynLib panglm, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom data.table as.data.table setDT setkeyv
#' @importFrom stats model.frame model.matrix model.response pnorm printCoefmat nobs
"_PACKAGE"

utils::globalVariables(c(".id", ".time", ".SD", "..cols", "..xnames", ":="))
