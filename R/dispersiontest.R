#' Pearson chi-squared/df overdispersion test
#'
#' A quick diagnostic for count models: the Pearson goodness-of-fit
#' statistic `sum((y - mu)^2 / mu)` divided by its residual degrees of
#' freedom should be close to 1 under a correctly-specified Poisson model;
#' values well above 1 indicate overdispersion (evidence for a negative
#' binomial specification instead). Requires `fitted.values`, so it's
#' available for `model = "pooling"`/`"within"` (poisson/negbin) and
#' `model = "random", family = "gaussian"`; random-effects poisson/negbin/
#' binomial fits don't currently expose `fitted.values` (see the package
#' vignette) so the test isn't available there.
#'
#' @param object a `"panglm"` fit with `family` `"poisson"` or `"negbin"`
#' @return a list of class `"panglm_dispersiontest"` with the Pearson
#'   statistic, its degrees of freedom, the ratio (statistic/df), and a
#'   one-sided p-value against `chisq(df)`
#' @export
panglm_dispersiontest <- function(object) {
  if (!inherits(object, "panglm")) stop("'object' must be a panglm fit", call. = FALSE)
  if (!object$family$family %in% c("poisson", "negbin")) {
    stop("panglm_dispersiontest() only applies to family = 'poisson' or 'negbin'", call. = FALSE)
  }
  if (is.null(object$fitted.values)) {
    stop("fitted.values are not available for this fit (model = '", object$model,
         "'), so the Pearson dispersion statistic cannot be computed", call. = FALSE)
  }

  mu <- object$fitted.values
  y <- object$y
  pearson_chisq <- sum((y - mu)^2 / mu)
  df <- object$df.residual
  ratio <- pearson_chisq / df
  pval <- stats::pchisq(pearson_chisq, df, lower.tail = FALSE)

  out <- list(statistic = pearson_chisq, parameter = df, ratio = ratio, p.value = pval)
  class(out) <- "panglm_dispersiontest"
  out
}

#' @export
print.panglm_dispersiontest <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\n\tPearson chi-squared/df overdispersion test\n\n")
  cat("chisq =", format(x$statistic, digits = digits),
      ", df =", x$parameter,
      ", ratio =", format(x$ratio, digits = digits),
      ", p-value =", format.pval(x$p.value, digits = digits), "\n")
  cat("alternative hypothesis: overdispersion (ratio > 1)\n\n")
  invisible(x)
}
