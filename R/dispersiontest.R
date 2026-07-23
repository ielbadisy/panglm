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
#'   statistic, its degrees of freedom, the ratio (statistic/df), a
#'   one-sided p-value against `chisq(df)`, and `n_excluded`: the number of
#'   observations dropped from the statistic because their fitted value was
#'   exactly zero or unavailable (e.g. rows in a structurally-zero panel
#'   individual screened out of a fixed-effects negbin fit)
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

  # Structurally-zero groups under a closed-form fixed-effects fit (e.g. an
  # all-zero panel individual under fit_within_negbin_dummy()) can produce
  # mu == 0 exactly (or NA, if the row was screened out of the fit
  # entirely), which turns a single Pearson term into 0/0 = NaN and poisons
  # the whole statistic. Exclude those observations and report how many.
  usable <- !is.na(mu) & mu > .Machine$double.eps
  n_excluded <- sum(!usable)
  mu <- mu[usable]; y <- y[usable]

  pearson_chisq <- sum((y - mu)^2 / mu)
  df <- object$df.residual - n_excluded
  ratio <- pearson_chisq / df
  pval <- stats::pchisq(pearson_chisq, df, lower.tail = FALSE)

  out <- list(statistic = pearson_chisq, parameter = df, ratio = ratio, p.value = pval,
              n_excluded = n_excluded)
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
  cat("alternative hypothesis: overdispersion (ratio > 1)\n")
  if (!is.null(x$n_excluded) && x$n_excluded > 0) {
    cat(x$n_excluded, "observation(s) with a zero or unavailable fitted value excluded\n")
  }
  cat("\n")
  invisible(x)
}
