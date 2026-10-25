# Further standard methods, so that the model classes share coef(),
# vcov(), confint() and nobs() (and logLik() where the fit maximises a
# likelihood): confint() for twopart(), logLik() for eventstudy(), and the
# coefficient methods of the hurdle model.

wald_confint <- function(b, V, parm, level) {
  tidy_check_level(level)
  if (missing(parm) || is.null(parm)) {
    parm <- names(b)
  } else if (is.numeric(parm)) {
    parm <- names(b)[parm]
  }
  if (anyNA(parm) || !all(parm %in% names(b))) {
    stop("'parm' must name or index coefficients of the fit", call. = FALSE)
  }
  se <- sqrt(diag(V))
  names(se) <- names(b)
  crit <- stats::qnorm(1 - (1 - level) / 2)
  out <- cbind(b[parm] - crit * se[parm], b[parm] + crit * se[parm])
  dimnames(out) <- list(parm, sprintf("%.1f %%", c((1 - level) / 2, 1 - (1 - level) / 2) * 100))
  out
}

#' Confidence intervals for a two-part model
#'
#' Wald intervals for the coefficients of both parts of a [twopart()] fit,
#' named as in `coef()` (`"binary:..."`, `"positive:..."`), from `vcov()`
#' (block diagonal for independent parts, the joint covariance for
#' `correlated = TRUE`).
#'
#' @param object a `"twopart"` fit
#' @param parm names or positions of the coefficients; all by default
#' @param level confidence level
#' @param ... unused
#' @return a two-column matrix
#' @family twopart methods
#' @examples
#' data(copd)
#' copd$cost <- ifelse(copd$exacerbations > 0, exp(6 + 0.1 * copd$crp), 0)
#' fit <- twopart(cost ~ treatment + crp, data = copd, index = c("id", "visit"))
#' confint(fit)
#' @export
confint.twopart <- function(object, parm, level = 0.95, ...) {
  wald_confint(stats::coef(object), stats::vcov(object), if (missing(parm)) NULL else parm,
               level)
}

#' @rdname eventstudy-methods
#' @export
logLik.panglm_eventstudy <- function(object, ...) stats::logLik(object$fit)

#' Coefficient methods for the hurdle model
#'
#' `coef()`, `vcov()`, `confint()`, `nobs()` and `summary()` for hurdle fits
#' (`panhurdle()`).
#' The coefficients of the zero part (conditional logit of `1(y > 0)`) and
#' of the count part (zero-truncated count model on the positive
#' observations) are named `"zero:..."` and `"count:..."`. The two parts are
#' estimated from separate factors of the hurdle likelihood, so `vcov()` is
#' block diagonal. No `logLik()` method is provided: the zero part maximises
#' a conditional likelihood and the count part an unconditional likelihood
#' with individual intercepts, so their sum is not the log-likelihood of one
#' model.
#'
#' @param object a `"panglm_hurdle"` fit
#' @param parm names or positions of the coefficients; all by default
#' @param level confidence level
#' @param x a `"summary.panglm_hurdle"` object
#' @param digits significant digits to print
#' @param ... unused
#' @return `coef()`: a named vector; `vcov()`: a matrix; `confint()`: a
#'   two-column matrix; `nobs()`: the number of observations; `summary()`:
#'   an object of class `"summary.panglm_hurdle"` with one coefficient table
#'   per part (`zero`, `count`: estimate, standard error, z value and
#'   p-value, from the inverse Hessian of each part) and the fit
#'   information.
#' @family panhurdle methods
#' @examples
#' data(copd)
#' fit <- panhurdle(exacerbations ~ crp, data = copd, index = c("id", "visit"))
#' coef(fit)
#' confint(fit)
#' summary(fit)
#' @name hurdle-methods
NULL

#' @rdname hurdle-methods
#' @export
coef.panglm_hurdle <- function(object, ...) {
  bz <- object$zero$coefficients
  bc <- object$count$coefficients
  stats::setNames(c(bz, bc), c(paste0("zero:", names(bz)), paste0("count:", names(bc))))
}

#' @rdname hurdle-methods
#' @export
vcov.panglm_hurdle <- function(object, ...) {
  Vz <- as.matrix(object$zero$vcov)
  Vc <- as.matrix(object$count$vcov)
  nm <- names(stats::coef(object))
  kz <- nrow(Vz); kc <- nrow(Vc)
  V <- matrix(0, kz + kc, kz + kc, dimnames = list(nm, nm))
  V[seq_len(kz), seq_len(kz)] <- Vz
  V[kz + seq_len(kc), kz + seq_len(kc)] <- Vc
  V
}

#' @rdname hurdle-methods
#' @export
nobs.panglm_hurdle <- function(object, ...) object$n

#' @rdname hurdle-methods
#' @export
confint.panglm_hurdle <- function(object, parm, level = 0.95, ...) {
  wald_confint(stats::coef(object), stats::vcov(object), if (missing(parm)) NULL else parm,
               level)
}

#' @rdname hurdle-methods
#' @export
summary.panglm_hurdle <- function(object, ...) {
  tab <- function(fit) {
    b <- fit$coefficients
    se <- sqrt(diag(as.matrix(fit$vcov)))
    z <- b / se
    cbind(Estimate = b, `Std. Error` = se, `z value` = z,
          `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  }
  out <- list(formula = object$formula, zero = tab(object$zero), count = tab(object$count),
              count_family = object$count_family, theta = object$count$theta,
              nobs = object$n, n_positive = object$n_positive,
              n_groups_zero = object$zero$n_used_groups,
              n_dropped_zero = object$zero$n_dropped_groups,
              n_groups_count = object$count$n_used_groups)
  class(out) <- "summary.panglm_hurdle"
  out
}

#' @rdname hurdle-methods
#' @export
print.summary.panglm_hurdle <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  family_label <- if (identical(x$count_family, "negbin")) "NB2" else "Poisson"
  cat("\nFormula: ", paste(deparse(x$formula), collapse = " "), "\n", sep = "")
  cat("\nEstimator: fixed-effects hurdle model (conditional logit and zero-truncated ",
      family_label, ")\n", sep = "")
  cat("N =", x$nobs, " positive =", x$n_positive, " vcov: model-based (inverse Hessian)\n")
  cat("\nZero part, 1(y > 0), conditional logit (groups used = ",
      or_na(x$n_groups_zero, NA_integer_), "):\n", sep = "")
  stats::printCoefmat(x$zero, digits = digits, signif.stars = TRUE, has.Pvalue = TRUE)
  cat("\nCount part, y | y > 0, zero-truncated ", family_label, " (groups used = ",
      or_na(x$n_groups_count, NA_integer_), "):\n", sep = "")
  stats::printCoefmat(x$count, digits = digits, signif.stars = TRUE, has.Pvalue = TRUE)
  if (!is.null(x$theta)) cat("Theta:", format(x$theta, digits = digits), "\n")
  cat("\n")
  invisible(x)
}

#' @rdname tidy-methods
#' @exportS3Method generics::tidy
tidy.panglm_hurdle <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  b <- stats::coef(x)
  tidy_wald(names(b), b, tidy_vcov_se(stats::vcov(x), names(b)), conf.int, conf.level,
            extra = data.frame(component = sub(":.*$", "", names(b)),
                               stringsAsFactors = FALSE))
}

#' @rdname tidy-methods
#' @exportS3Method generics::glance
glance.panglm_hurdle <- function(x, ...) {
  data.frame(count.family = x$count_family, nobs = x$n, n.positive = x$n_positive,
             n.groups.count = or_na(x$count$n_used_groups, NA_integer_),
             stringsAsFactors = FALSE)
}
