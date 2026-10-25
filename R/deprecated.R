#' Deprecated functions in panglm
#'
#' These functions still work but will be removed in a future release.
#' Each one emits a deprecation warning once per session and calls its
#' replacement, so it returns exactly the same result.
#'
#' Exported function names in panglm are glued lowercase, without
#' underscores:
#'
#' \tabular{ll}{
#'   \strong{Deprecated} \tab \strong{Replacement} \cr
#'   \code{gaussian_family()} \tab \code{\link{gaussianfamily}()} \cr
#'   \code{binomial_family()} \tab \code{\link{binomialfamily}()} \cr
#'   \code{poisson_family()} \tab \code{\link{poissonfamily}()} \cr
#'   \code{negbin_family()} \tab \code{\link{negbinfamily}()} \cr
#'   \code{gamma_family()} \tab \code{\link{gammafamily}()} \cr
#'   \code{beta_family()} \tab \code{\link{betafamily}()} \cr
#'   \code{ordered_family()} \tab \code{\link{orderedfamily}()} \cr
#'   \code{tobit_family()} \tab \code{\link{tobitfamily}()} \cr
#'   \code{panglm_dispersiontest()} \tab \code{\link{pandispersiontest}()} \cr
#'   \code{panglm_hausman()} \tab \code{\link{panhausman}()} \cr
#'   \code{panglm_hurdle()} \tab \code{\link{panhurdle}()} \cr
#' }
#'
#' \code{panglm_hurdle(count_family = )} is \code{panhurdle(family = )}.
#' The objects returned keep their classes (\code{"panglm_dispersiontest"},
#' \code{"panglm_hausman"}, \code{"panglm_hurdle"}), so code that tests
#' or dispatches on them is unaffected.
#'
#' @param link,left,right passed to the replacement
#' @param object,fe,re,method,data,vcov,formula,index,maxit,tol passed to the
#'   replacement
#' @param count_family passed as \code{family} to \code{\link{panhurdle}()}
#' @return the value of the replacement function
#' @name panglm-deprecated
#' @keywords internal
NULL

# Names already warned about in this session.
.panglm_deprecated <- new.env(parent = emptyenv())

#' Warn once per session that a function or an argument is deprecated
#'
#' @param old,new names of the deprecated item and of its replacement
#' @param what `"function"` or `"argument"`
#' @keywords internal
#' @noRd
deprecate_once <- function(old, new, what = c("function", "argument")) {
  what <- match.arg(what)
  key <- paste(what, old)
  if (isTRUE(.panglm_deprecated[[key]])) return(invisible(FALSE))
  assign(key, TRUE, envir = .panglm_deprecated)
  msg <- if (what == "function") {
    sprintf("`%s` is deprecated; use `%s` instead.", old, new)
  } else {
    sprintf("argument %s is deprecated; use %s instead.", old, new)
  }
  warning(warningCondition(
    paste(msg, "This warning is shown once per session."),
    class = c("deprecatedWarning", "panglm_deprecated")
  ))
  invisible(TRUE)
}

#' Forget which deprecation warnings were shown (used by the tests)
#' @keywords internal
#' @noRd
deprecate_reset <- function() {
  rm(list = ls(.panglm_deprecated, all.names = TRUE), envir = .panglm_deprecated)
  invisible(NULL)
}

#' @rdname panglm-deprecated
#' @export
gaussian_family <- function(link = "identity") {
  deprecate_once("gaussian_family()", "gaussianfamily()")
  gaussianfamily(link = link)
}

#' @rdname panglm-deprecated
#' @export
binomial_family <- function(link = c("logit", "probit")) {
  deprecate_once("binomial_family()", "binomialfamily()")
  binomialfamily(link = link)
}

#' @rdname panglm-deprecated
#' @export
poisson_family <- function(link = "log") {
  deprecate_once("poisson_family()", "poissonfamily()")
  poissonfamily(link = link)
}

#' @rdname panglm-deprecated
#' @export
negbin_family <- function(link = "log") {
  deprecate_once("negbin_family()", "negbinfamily()")
  negbinfamily(link = link)
}

#' @rdname panglm-deprecated
#' @export
gamma_family <- function(link = "log") {
  deprecate_once("gamma_family()", "gammafamily()")
  gammafamily(link = link)
}

#' @rdname panglm-deprecated
#' @export
beta_family <- function(link = c("logit", "probit")) {
  deprecate_once("beta_family()", "betafamily()")
  betafamily(link = link)
}

#' @rdname panglm-deprecated
#' @export
ordered_family <- function(link = c("logit", "probit")) {
  deprecate_once("ordered_family()", "orderedfamily()")
  orderedfamily(link = link)
}

#' @rdname panglm-deprecated
#' @export
tobit_family <- function(left = 0, right = Inf) {
  deprecate_once("tobit_family()", "tobitfamily()")
  tobitfamily(left = left, right = right)
}

#' @rdname panglm-deprecated
#' @export
panglm_dispersiontest <- function(object) {
  deprecate_once("panglm_dispersiontest()", "pandispersiontest()")
  pandispersiontest(object)
}

#' @rdname panglm-deprecated
#' @export
panglm_hausman <- function(fe, re = NULL, method = c("classical", "mundlak"), data = NULL,
                           vcov = c("cluster", "HC1", "classical")) {
  deprecate_once("panglm_hausman()", "panhausman()")
  panhausman(fe, re = re, method = method, data = data, vcov = vcov)
}

#' @rdname panglm-deprecated
#' @export
panglm_hurdle <- function(formula, data, index,
                          count_family = c("poisson", "negbin"),
                          maxit = 100, tol = 1e-10) {
  deprecate_once("panglm_hurdle()", "panhurdle()")
  panhurdle(formula, data, index, family = count_family, maxit = maxit, tol = tol)
}
