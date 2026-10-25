#' Summarize a panglm fit
#'
#' Produces a classic coefficient table (Estimate, Std. Error, z value,
#' Pr(>|z|)) in the style of [summary.glm()].
#'
#' @param object a `"panglm"` object
#' @param ... unused
#' @return an object of class `"summary.panglm"`, a list with components
#'   `call` (the matched call), `model`, `effect`, and `family` (as stored on
#'   the fit), `coefficients` (a numeric matrix with one row per parameter
#'   and columns `Estimate`, `Std. Error`, `z value`, `Pr(>|z|)`), `nobs`,
#'   `n_groups`, `df.residual`, `loglik`, `dispersion` (`NULL` when not
#'   applicable, as for negative binomial fits), `iterations`, and `vcov_type` (which covariance estimator
#'   produced the reported standard errors). It has an associated
#'   `print.summary.panglm()` method and is normally printed rather than
#'   used programmatically. Also included are `sigma` (random-intercept
#'   standard deviation) and `shape` (gamma shape) when estimated,
#'   `varcomp` (covariance matrix of the random intercept and slopes) for
#'   random-slope fits, for `family = "negbin"` the NB2 `negbin_theta`
#'   (variance `mu + mu^2 / theta`) of pooled, fixed-effects and random-slope
#'   fits, with its standard error `negbin_theta_se` (pooled and
#'   fixed-effects fits, from the observed information at the fitted means,
#'   as in `MASS::glm.nb()`), or `negbin_ab`, the beta mixing parameters
#'   `a` and `b` of the random-intercept model, and, for
#'   `cre = TRUE` fits, `cre_test`: the Wald (Mundlak) test that the
#'   individual-mean coefficients are jointly zero, computed with the
#'   covariance stored on the fit.
#' @family panglm methods
#' @examples
#' data(copd)
#' fit <- panglm(fev1 ~ treatment + age + smoker + crp, data = copd,
#'               index = c("id", "visit"), model = "pooling", family = "gaussian")
#' summary(fit)
#' @export
summary.panglm <- function(object, ...) {
  coefs <- object$coefficients
  se <- sqrt(diag(object$vcov))
  zval <- coefs / se
  pval <- 2 * stats::pnorm(-abs(zval))

  coeftable <- cbind(Estimate = coefs, `Std. Error` = se,
                      `z value` = zval, `Pr(>|z|)` = pval)
  rownames(coeftable) <- names(coefs)

  nb <- negbin_parameters(object)
  out <- list(call = object$call, model = object$model, effect = object$effect, family = object$family,
              coefficients = coeftable, nobs = object$nobs, n_groups = object$n_groups,
              df.residual = object$df.residual, loglik = object$loglik,
              dispersion = if (is.null(nb)) object$dispersion, iterations = object$iterations,
              vcov_type = object$vcov_type,
              sigma = if (is.null(object$varcomp)) object$sigma, shape = object$shape,
              negbin_theta = nb$theta, negbin_theta_se = nb$se,
              negbin_ab = nb$ab,
              varcomp = object$varcomp,
              weighted = !is.null(object$weights), cre_test = cre_wald_test(object))
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
  if (!is.null(x$sigma) && identical(x$model, "random")) {
    cat("Random-intercept SD:", format(x$sigma, digits = digits), "\n")
  }
  if (!is.null(x$varcomp)) {
    sds <- sqrt(diag(x$varcomp))
    cat("Random effects (SD):", paste0(names(sds), " ", format(sds, digits = digits),
                                       collapse = ", "), "\n")
    if (length(sds) > 1L) {
      cors <- stats::cov2cor(x$varcomp)
      pairs <- which(lower.tri(cors), arr.ind = TRUE)
      cat("Random-effects correlation:",
          paste0(rownames(cors)[pairs[, 1]], " ~ ", colnames(cors)[pairs[, 2]], " ",
                 format(cors[pairs], digits = digits), collapse = ", "), "\n")
    }
  }
  if (!is.null(x$shape)) cat("Gamma shape:", format(x$shape, digits = digits), "\n")
  if (!is.null(x$negbin_theta)) {
    cat("Negative binomial theta (NB2):", format(x$negbin_theta, digits = digits),
        if (!is.null(x$negbin_theta_se) && is.finite(x$negbin_theta_se)) {
          paste0("(SE ", format(x$negbin_theta_se, digits = digits), ")")
        }, "\n")
  }
  if (!is.null(x$negbin_ab)) {
    cat("Beta mixing parameters (Hausman, Hall and Griliches): a =",
        format(x$negbin_ab[["a"]], digits = digits), " b =",
        format(x$negbin_ab[["b"]], digits = digits), "\n")
  }
  if (isTRUE(x$weighted)) cat("Observation weights: yes (weighted pseudo-likelihood)\n")
  if (!is.null(x$cre_test)) {
    cat("Mundlak test (individual means jointly zero): chi2(", x$cre_test$df, ") = ",
        format(x$cre_test$statistic, digits = digits), ", p = ",
        format.pval(x$cre_test$p.value, digits = digits), "\n", sep = "")
  }
  # the residual df is not defined for every estimator (e.g. the conditional
  # logit), so it is only printed when it is a number
  df_line <- c(if (length(x$df.residual) == 1L && is.finite(x$df.residual)) {
    paste("Residual df:", x$df.residual)
  }, if (length(x$iterations) == 1L && !is.na(x$iterations)) {
    paste("Iterations:", x$iterations)
  })
  if (length(df_line)) cat(paste(df_line, collapse = "  "), "\n", sep = "")
  cat("\n")
  invisible(x)
}

#' Negative binomial parameters of a panglm fit
#'
#' The NB2 `theta` of pooled and fixed-effects fits (stored as `theta` by
#' the estimators) and of random-slope fits (`negbin_theta`), with the
#' standard error of `theta` from its observed information at the fitted
#' means, as `MASS::glm.nb()` reports it; the beta mixing parameters `a`
#' and `b` of the random-intercept model. `NULL` for other families.
#'
#' @keywords internal
#' @noRd
negbin_parameters <- function(object) {
  if (!identical(object$family$family, "negbin")) return(NULL)
  if (!is.null(object$negbin_theta)) return(list(theta = object$negbin_theta))
  if (identical(object$model, "random")) {
    if (is.null(object$negbin_a)) return(list())
    return(list(ab = c(a = object$negbin_a, b = object$negbin_b)))
  }
  theta <- object$theta
  if (length(theta) != 1L || !is.finite(theta)) return(list())
  list(theta = theta, se = negbin_theta_se(object$y, object$fitted.values, theta,
                                           object$weights))
}

negbin_theta_se <- function(y, mu, theta, weights = NULL) {
  w <- if (is.null(weights)) rep(1, length(y)) else weights
  if (length(mu) != length(y)) return(NA_real_)
  ok <- is.finite(mu) & is.finite(y) & w > 0
  y <- y[ok]; mu <- mu[ok]; w <- w[ok]
  info <- sum(w * (-trigamma(theta + y) + trigamma(theta) - 1 / theta +
                     2 / (mu + theta) - (y + theta) / (mu + theta)^2))
  if (is.finite(info) && info > 0) 1 / sqrt(info) else NA_real_
}

#' Wald test that the correlated random effects means are jointly zero
#'
#' @keywords internal
#' @noRd
cre_wald_test <- function(object) {
  columns <- object$cre_columns %||% character(0)
  if (!length(columns)) return(NULL)
  b <- object$coefficients[columns]
  V <- object$vcov[columns, columns, drop = FALSE]
  statistic <- tryCatch(as.numeric(t(b) %*% solve(V, b)), error = function(e) NA_real_)
  list(statistic = statistic, df = length(columns),
       p.value = stats::pchisq(statistic, length(columns), lower.tail = FALSE))
}
