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
#' @return a one-row `data.frame`
#' @exportS3Method generics::glance
glance.panglm <- function(x, ...) {
  data.frame(
    model = x$model,
    effect = x$effect,
    family = x$family$family,
    link = x$family$link,
    vcov.type = x$vcov_type,
    logLik = if (is.null(x$loglik)) NA_real_ else x$loglik,
    df.residual = x$df.residual,
    nobs = x$nobs,
    n.groups = x$n_groups,
    stringsAsFactors = FALSE
  )
}
