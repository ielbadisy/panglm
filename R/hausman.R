#' Hausman specification test (fixed vs. random effects)
#'
#' Tests the null that the random-effects estimator is consistent (i.e.
#' individual effects are uncorrelated with the regressors) against the
#' alternative that only the fixed-effects estimator is consistent, using
#' the classic Hausman (1978) statistic
#' `t(b_fe - b_re) %*% solve(V_fe - V_re) %*% (b_fe - b_re)`,
#' computed over the coefficients common to both models (the random-effects
#' model's intercept is dropped, since the within/FE model has none).
#'
#' @param fe a `"panglm"` fit with `model = "within"`
#' @param re a `"panglm"` fit with `model = "random"`, same formula/data/family
#' @return a list with class `"panglm_hausman"`: `statistic`, `parameter`
#'   (df), `p.value`, and the compared coefficient vectors
#' @examples
#' data(copd)
#' fit_fe <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
#'                   model = "within", family = "gaussian")
#' fit_re <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
#'                   model = "random", family = "gaussian")
#' panglm_hausman(fit_fe, fit_re)
#' @export
panglm_hausman <- function(fe, re) {
  if (!inherits(fe, "panglm") || fe$model != "within") {
    stop("'fe' must be a panglm fit with model = 'within'", call. = FALSE)
  }
  if (!inherits(re, "panglm") || re$model != "random") {
    stop("'re' must be a panglm fit with model = 'random'", call. = FALSE)
  }
  if (fe$family$family != re$family$family) {
    stop("'fe' and 're' must use the same family", call. = FALSE)
  }

  common <- intersect(names(fe$coefficients), names(re$coefficients))
  common <- setdiff(common, "(Intercept)")
  if (length(common) == 0) stop("no coefficients in common between 'fe' and 're'", call. = FALSE)

  b_fe <- fe$coefficients[common]
  b_re <- re$coefficients[common]
  V_fe <- fe$vcov[common, common, drop = FALSE]
  V_re <- re$vcov[common, common, drop = FALSE]

  diff <- b_fe - b_re
  V_diff <- V_fe - V_re
  V_diff_inv <- tryCatch(solve(V_diff), error = function(e) MASS_ginv(V_diff))

  stat <- as.numeric(t(diff) %*% V_diff_inv %*% diff)
  df <- length(common)
  pval <- stats::pchisq(stat, df, lower.tail = FALSE)

  out <- list(statistic = stat, parameter = df, p.value = pval,
              coef_fe = b_fe, coef_re = b_re, common = common)
  class(out) <- "panglm_hausman"
  out
}

MASS_ginv <- function(x) {
  s <- svd(x)
  tol <- max(dim(x)) * max(s$d) * .Machine$double.eps
  pos <- s$d > tol
  if (!any(pos)) return(matrix(0, nrow(x), ncol(x)))
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

#' @export
print.panglm_hausman <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\n\tHausman specification test (fixed vs. random effects)\n\n")
  cat("chisq =", format(x$statistic, digits = digits),
      ", df =", x$parameter,
      ", p-value =", format.pval(x$p.value, digits = digits), "\n")
  cat("alternative hypothesis: one model is inconsistent\n\n")
  invisible(x)
}
