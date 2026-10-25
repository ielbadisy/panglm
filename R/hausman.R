#' Hausman specification test (fixed vs. random effects)
#'
#' Tests the null that the random-effects estimator is consistent (i.e.
#' individual effects are uncorrelated with the regressors) against the
#' alternative that only the fixed-effects estimator is consistent.
#'
#' `method = "classical"` computes the Hausman (1978) statistic
#' `t(b_fe - b_re) %*% solve(V_fe - V_re) %*% (b_fe - b_re)` over the
#' coefficients common to both models (the random-effects intercept is
#' dropped), with a chi-square reference. It always uses the model-based
#' (classical) covariances, also when `fe` was fit with a robust `vcov`: the
#' contrast variance `V_fe - V_re` relies on the random-effects estimator
#' being efficient under the null, which holds only for model-based
#' covariances and independent, homoskedastic idiosyncratic errors.
#'
#' `method = "mundlak"` is the regression-based (correlated random effects)
#' version of the same test (Mundlak 1978; Wooldridge 2010, section 10.7.3),
#' valid under heteroskedasticity and serial correlation. The pooled model
#' `y ~ x + xbar` is fit with the individual means `xbar` of the time-varying
#' regressors of `fe`, and the null `coef(xbar) = 0` is tested with a Wald
#' statistic using the `vcov` covariance. With `vcov = "cluster"` the
#' reference is `F(k, G - 1)` (`G` individuals), otherwise chi-square(k).
#' Gaussian family only.
#'
#' @param fe a `"panglm"` fit with `model = "within"`
#' @param re a `"panglm"` fit with `model = "random"`, same formula/data/family;
#'   required for `method = "classical"`, ignored for `"mundlak"`
#' @param method `"classical"` (default) or `"mundlak"`
#' @param data for `method = "mundlak"`: the data `fe` was fit on
#' @param vcov for `method = "mundlak"`: covariance of the pooled fit, one of
#'   `"cluster"` (default), `"HC1"`, `"classical"`
#' @return a list with class `"panglm_hausman"`: `statistic`, `parameter`
#'   (df; two values for an F reference), `p.value`, `method`,
#'   `distribution` (`"chisq"` or `"F"`), and the compared coefficient
#'   vectors (`coef_fe`, `coef_re` for `"classical"`; `coef_mean` for
#'   `"mundlak"`)
#' @family tests
#' @seealso [panglm()] with `cre = TRUE`, whose summary reports the Mundlak
#'   test.
#' @examples
#' data(copd)
#' fit_fe <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
#'                   model = "within", family = "gaussian")
#' fit_re <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
#'                   model = "random", family = "gaussian")
#' panhausman(fit_fe, fit_re)
#' panhausman(fit_fe, method = "mundlak", data = copd)
#' @export
panhausman <- function(fe, re = NULL, method = c("classical", "mundlak"), data = NULL,
                           vcov = c("cluster", "HC1", "classical")) {
  method <- match.arg(method)
  if (!inherits(fe, "panglm") || fe$model != "within") {
    stop("'fe' must be a panglm fit with model = 'within'", call. = FALSE)
  }
  if (method == "mundlak") return(hausman_mundlak(fe, data, match.arg(vcov)))
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
  V_fe <- vcov.panglm(fe, type = "classical")[common, common, drop = FALSE]
  V_re <- vcov.panglm(re, type = "classical")[common, common, drop = FALSE]

  diff <- b_fe - b_re
  V_diff <- V_fe - V_re
  V_diff_inv <- tryCatch(solve(V_diff), error = function(e) MASS_ginv(V_diff))

  stat <- as.numeric(t(diff) %*% V_diff_inv %*% diff)
  df <- length(common)
  pval <- stats::pchisq(stat, df, lower.tail = FALSE)

  out <- list(statistic = stat, parameter = df, p.value = pval, method = "classical",
              distribution = "chisq", coef_fe = b_fe, coef_re = b_re, common = common)
  class(out) <- "panglm_hausman"
  out
}

hausman_mundlak <- function(fe, data, vcov) {
  if (is.null(data)) stop("method = 'mundlak' needs the 'data' that 'fe' was fit on", call. = FALSE)
  if (fe$family$family != "gaussian") {
    stop("method = 'mundlak' is only available for the gaussian family", call. = FALSE)
  }
  if (!identical(fe$effect, "individual")) {
    stop("method = 'mundlak' is only available for effect = 'individual'", call. = FALSE)
  }
  index <- fe$index
  vars <- unique(c(all.vars(fe$formula), index))
  data <- as.data.frame(data)[stats::complete.cases(as.data.frame(data)[, vars, drop = FALSE]), , drop = FALSE]
  mf <- stats::model.frame(fe$formula, data = data)
  X <- stats::model.matrix(fe$formula, mf)
  slopes <- intersect(names(fe$coefficients), colnames(X))
  if (!length(slopes)) stop("no time-varying regressors in 'fe'", call. = FALSE)
  X <- X[, slopes, drop = FALSE]
  id <- data[[index[1]]]
  xbar <- apply(X, 2, function(v) stats::ave(v, id))
  xbar <- matrix(xbar, nrow = nrow(X), dimnames = list(NULL, colnames(X)))
  k <- length(slopes)
  xn <- paste0(".x", seq_len(k)); mn <- paste0(".xbar", seq_len(k))
  dd <- data.frame(.y = stats::model.response(mf), X, xbar, data[, index, drop = FALSE])
  names(dd) <- c(".y", xn, mn, index)
  fit <- panglm(stats::reformulate(c(xn, mn), ".y"), data = dd, index = index,
                model = "pooling", family = "gaussian", vcov = vcov)
  b <- fit$coefficients[mn]
  V <- fit$vcov[mn, mn, drop = FALSE]
  W <- as.numeric(t(b) %*% solve(V, b))
  names(b) <- slopes
  if (vcov == "cluster") {
    G <- fit$n_groups
    stat <- W / k; df <- c(k, G - 1)
    pval <- stats::pf(stat, k, G - 1, lower.tail = FALSE)
    dist <- "F"
  } else {
    stat <- W; df <- k
    pval <- stats::pchisq(W, k, lower.tail = FALSE)
    dist <- "chisq"
  }
  out <- list(statistic = stat, parameter = df, p.value = pval, method = "mundlak",
              distribution = dist, vcov_type = vcov, coef_mean = b, common = slopes)
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
  if (identical(x$method, "mundlak")) {
    cat("\n\tMundlak specification test (fixed vs. random effects), vcov:", x$vcov_type, "\n\n")
  } else {
    cat("\n\tHausman specification test (fixed vs. random effects)\n\n")
  }
  if (identical(x$distribution, "F")) {
    cat("F =", format(x$statistic, digits = digits),
        ", df =", paste(x$parameter, collapse = ", "),
        ", p-value =", format.pval(x$p.value, digits = digits), "\n")
  } else {
    cat("chisq =", format(x$statistic, digits = digits),
        ", df =", x$parameter,
        ", p-value =", format.pval(x$p.value, digits = digits), "\n")
  }
  cat("alternative hypothesis: one model is inconsistent\n\n")
  invisible(x)
}
