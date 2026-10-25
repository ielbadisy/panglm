# vcov(), nobs(), confint() and summary() for eventstudy() fits. The event
# coefficients are labelled by event time, as in coef(); the reference
# period is a coefficient fixed at zero (zero variance) and an event time
# that the data do not identify has a missing estimate and variance.

es_terms <- function(events) {
  ifelse(events < 0, paste0("event_lead", -events), paste0("event_lag", events))
}

#' Methods for event-study fits
#'
#' `vcov()`, `nobs()`, `confint()` and `summary()` for the fits of
#' [eventstudy()]. The coefficients are labelled by event time, as in
#' `coef()`. The reference period enters as a coefficient fixed at zero, with
#' zero variance and covariance and a degenerate interval at zero; an event
#' time that the data do not identify has missing values throughout.
#' Intervals and tests are Wald (normal) based, from the variance of the
#' underlying two-way fixed-effects fit (`object$fit`).
#'
#' @param object an `"panglm_eventstudy"` fit from [eventstudy()]
#' @param parm event times (as numbers or names) or positions of the
#'   coefficients; all by default
#' @param level confidence level
#' @param x a `"summary.panglm_eventstudy"` object
#' @param digits significant digits to print
#' @param signif.stars passed to [stats::printCoefmat()]
#' @param ... unused
#' @return `vcov()`: a square matrix with one row per event time; `nobs()`:
#'   the number of observations of the fit; `confint()`: a two-column matrix;
#'   `summary()`: an object of class `"summary.panglm_eventstudy"` with the
#'   coefficient table `coefficients` (estimate, standard error, z value,
#'   p-value), the intervals `conf.int`, the pre-trend test `pretrend` and
#'   the fit description.
#' @family eventstudy methods
#' @examples
#' set.seed(1)
#' d <- expand.grid(id = 1:200, year = 2015:2022)
#' d$first <- ifelse(d$id <= 100, 2019, NA)
#' d$visits <- 3 + rnorm(200)[d$id] + 0.1 * (d$year - 2015) +
#'   ifelse(!is.na(d$first) & d$year >= d$first, 0.5, 0) + rnorm(nrow(d))
#' es <- eventstudy(visits ~ 1, data = d, index = c("id", "year"),
#'                  treat = "first", leads = 3, lags = 2)
#' summary(es)
#' confint(es)
#' vcov(es)
#' @name eventstudy-methods
NULL

#' @rdname eventstudy-methods
#' @export
vcov.panglm_eventstudy <- function(object, ...) {
  ev <- object$estimates$event
  nm <- as.character(ev)
  terms <- es_terms(ev)
  V0 <- stats::vcov(object$fit)
  V <- matrix(NA_real_, length(ev), length(ev), dimnames = list(nm, nm))
  kept <- terms %in% rownames(V0)
  V[kept, kept] <- V0[terms[kept], terms[kept], drop = FALSE]
  ref <- ev == object$ref
  V[ref, ] <- 0
  V[, ref] <- 0
  V[ref, !is.finite(object$estimates$estimate)] <- NA_real_
  V[!is.finite(object$estimates$estimate), ref] <- NA_real_
  V
}

#' @rdname eventstudy-methods
#' @export
nobs.panglm_eventstudy <- function(object, ...) object$fit$nobs

#' @rdname eventstudy-methods
#' @export
confint.panglm_eventstudy <- function(object, parm, level = 0.95, ...) {
  tidy_check_level(level)
  b <- stats::coef(object)
  se <- object$estimates$se
  names(se) <- names(b)
  if (missing(parm)) {
    parm <- names(b)
  } else if (is.numeric(parm) && !all(as.character(parm) %in% names(b))) {
    parm <- names(b)[parm]
  } else {
    parm <- as.character(parm)
  }
  if (anyNA(parm) || !all(parm %in% names(b))) {
    stop("'parm' must name event times of the fit", call. = FALSE)
  }
  crit <- stats::qnorm(1 - (1 - level) / 2)
  out <- cbind(b[parm] - crit * se[parm], b[parm] + crit * se[parm])
  dimnames(out) <- list(parm, sprintf("%.1f %%", c((1 - level) / 2, 1 - (1 - level) / 2) * 100))
  out
}

#' @rdname eventstudy-methods
#' @export
summary.panglm_eventstudy <- function(object, level = object$level, ...) {
  tidy_check_level(level)
  b <- stats::coef(object)
  se <- object$estimates$se
  ok <- is.finite(se) & se > 0
  z <- ifelse(ok, b / se, NA_real_)
  tab <- cbind(Estimate = b, `Std. Error` = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  rownames(tab) <- names(b)
  out <- list(coefficients = tab, conf.int = stats::confint(object, level = level),
              pretrend = object$pretrend, ref = object$ref, leads = object$leads,
              lags = object$lags, bin = object$bin, level = level,
              family = object$fit$family$family, vcov_type = object$fit$vcov_type,
              nobs = object$fit$nobs, n_groups = object$fit$n_groups,
              call = object$fit$call)
  class(out) <- "summary.panglm_eventstudy"
  out
}

#' @rdname eventstudy-methods
#' @export
print.summary.panglm_eventstudy <- function(x, digits = max(3, getOption("digits") - 3),
                                            signif.stars = TRUE, ...) {
  cat("\nEvent study (two-way fixed effects, ", x$family, ")\n", sep = "")
  cat("N =", x$nobs, " groups =", x$n_groups, " vcov:", x$vcov_type, "\n")
  cat("Reference event time:", x$ref, if (x$bin) "(end points binned)" else "", "\n")
  cat("\nEvent-time coefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
                      na.print = "", P.values = TRUE, has.Pvalue = TRUE)
  if (!is.null(x$pretrend)) {
    cat("\nPre-trend test: Chi-squared = ", format(x$pretrend$statistic, digits = digits),
        ", df = ", x$pretrend$df, ", p-value = ",
        format.pval(x$pretrend$p.value, digits = digits), "\n", sep = "")
  }
  invisible(x)
}
