#' Tidy a panglm fit into a data frame of coefficients
#'
#' Implements the `generics::tidy()` generic so `panglm` results work with
#' `broom`/`modelsummary` table-building workflows.
#'
#' @param x a `"panglm"` fit
#' @param conf.int if `TRUE`, add `conf.low`/`conf.high` columns
#' @param conf.level confidence level for the interval
#' @param ... unused
#' @return a `data.frame` with one row per coefficient
#' @examples
#' data(copd)
#' fit <- panglm(fev1 ~ treatment + age + smoker + crp, data = copd,
#'               index = c("id", "visit"), model = "pooling", family = "gaussian")
#' if (requireNamespace("generics", quietly = TRUE)) {
#'   generics::tidy(fit, conf.int = TRUE)
#' }
#' @exportS3Method generics::tidy
tidy.panglm <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  coefs <- x$coefficients
  se <- sqrt(diag(x$vcov))
  zval <- coefs / se
  pval <- 2 * stats::pnorm(-abs(zval))

  out <- data.frame(
    term = names(coefs),
    estimate = as.numeric(coefs),
    std.error = as.numeric(se),
    statistic = as.numeric(zval),
    p.value = as.numeric(pval),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  if (conf.int) {
    crit <- stats::qnorm(1 - (1 - conf.level) / 2)
    out$conf.low <- out$estimate - crit * out$std.error
    out$conf.high <- out$estimate + crit * out$std.error
  }
  out
}

#' One-row model summary of a panglm fit
#'
#' Implements the `generics::glance()` generic.
#'
#' @param x a `"panglm"` fit
#' @param ... unused
#' @return a one-row `data.frame` including log-likelihood, AIC, BIC, the
#'   likelihood parameter count when a likelihood is available, and
#'   (Gaussian-family fits only) R-squared and adjusted R-squared
#' @examples
#' data(copd)
#' fit <- panglm(fev1 ~ treatment + age + smoker + crp, data = copd,
#'               index = c("id", "visit"), model = "pooling", family = "gaussian")
#' if (requireNamespace("generics", quietly = TRUE)) {
#'   generics::glance(fit)
#' }
#' @exportS3Method generics::glance
glance.panglm <- function(x, ...) {
  ll <- stats::logLik(x)
  has_likelihood <- length(ll) == 1L && is.finite(as.numeric(ll))

  ## R-squared/adjusted R-squared: defined here as the usual proportion of
  ## variance explained, 1 - RSS/TSS, computed on whatever scale `fitted()`
  ## and `y` are on (the linear predictor for gaussian identity-link models,
  ## i.e. the standard OLS R-squared; the same ratio for other families,
  ## which is a Gaussian-style pseudo-R-squared, not McFadden's or another
  ## likelihood-based pseudo-R-squared - reported as NA when TSS is 0 or the
  ## fit has no residual df left to adjust for).
  r2 <- adj_r2 <- NA_real_
  if (!is.null(x$fitted.values) && !is.null(x$y)) {
    yv <- x$y; fv <- x$fitted.values
    tss <- sum((yv - mean(yv))^2)
    if (is.finite(tss) && tss > 0) {
      rss <- sum((yv - fv)^2)
      r2 <- 1 - rss / tss
      k <- length(x$coefficients)
      n <- length(yv)
      if (n - k > 0) adj_r2 <- 1 - (1 - r2) * (n - 1) / (n - k)
    }
  }

  data.frame(
    model = x$model,
    effect = x$effect,
    family = x$family$family,
    link = x$family$link,
    vcov.type = x$vcov_type,
    logLik = if (has_likelihood) as.numeric(ll) else NA_real_,
    AIC = if (has_likelihood) stats::AIC(x) else NA_real_,
    BIC = if (has_likelihood) stats::BIC(x) else NA_real_,
    df.logLik = attr(ll, "df"),
    df.residual = x$df.residual,
    r.squared = r2,
    adj.r.squared = adj_r2,
    nobs = x$nobs,
    n.groups = x$n_groups,
    stringsAsFactors = FALSE
  )
}
