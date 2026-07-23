#' Summarize a panglm fit
#'
#' Produces a classic coefficient table (Estimate, Std. Error, z value,
#' Pr(>|z|)) in the style of [summary.glm()].
#'
#' @param object a `"panglm"` object
#' @param ... unused
#' @export
summary.panglm <- function(object, ...) {
  coefs <- object$coefficients
  se <- sqrt(diag(object$vcov))
  zval <- coefs / se
  pval <- 2 * stats::pnorm(-abs(zval))

  coeftable <- cbind(Estimate = coefs, `Std. Error` = se,
                      `z value` = zval, `Pr(>|z|)` = pval)
  rownames(coeftable) <- names(coefs)

  out <- list(call = object$call, model = object$model, effect = object$effect, family = object$family,
              coefficients = coeftable, nobs = object$nobs, n_groups = object$n_groups,
              df.residual = object$df.residual, loglik = object$loglik,
              dispersion = object$dispersion, iterations = object$iterations,
              vcov_type = object$vcov_type)
  class(out) <- "summary.panglm"
  out
}

#' @export
print.summary.panglm <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nFamily:", x$family$family, " Link:", x$family$link,
      " Model:", x$model, " Effect:", x$effect, "\n")
  cat("N =", x$nobs, " groups =", x$n_groups,
      " vcov:", x$vcov_type, "\n\n")
  cat("Coefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = TRUE, has.Pvalue = TRUE)
  cat("\n")
  if (!is.na(x$loglik)) cat("Log-likelihood:", format(x$loglik, digits = digits), "\n")
  if (!is.null(x$dispersion)) cat("Dispersion parameter:", format(x$dispersion, digits = digits), "\n")
  cat("Residual df:", x$df.residual, " Iterations:", x$iterations, "\n\n")
  invisible(x)
}
