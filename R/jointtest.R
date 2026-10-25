#' Joint Wald test of coefficients
#'
#' Tests `H0: R beta = r` with the statistic
#' `(R b - r)' (R V R')^(-1) (R b - r)`, referred to a chi-squared
#' distribution with `nrow(R)` degrees of freedom. The coefficients are
#' chosen by name (`terms`), by a regular expression (`pattern`) or through
#' an explicit restriction matrix `R`. Typical uses are the joint test of
#' the spline columns of a regressor (does `x` matter at all), of the
#' Mundlak means of a correlated random effects fit, or of the leads of an
#' event study.
#'
#' @param object a `"panglm"` fit
#' @param terms character vector of coefficient names
#' @param pattern regular expression selecting coefficient names (used when
#'   `terms` and `R` are `NULL`)
#' @param R optional restriction matrix with one column per coefficient
#' @param r right-hand side of the restrictions (recycled), default 0
#' @param vcov `NULL` (the covariance stored on the fit) or a type passed to
#'   [vcov.panglm()], e.g. `"cluster"`
#' @return an object of class `"panglm_jointtest"`: a list with
#'   `statistic`, `df`, `p.value`, `terms` and `vcov_type`
#' @examples
#' data(copd)
#' fit <- panglm(fev1 ~ splines::ns(crp, 3), data = copd,
#'               index = c("id", "visit"), model = "within", vcov = "cluster")
#' jointtest(fit, pattern = "crp")
#' @export
jointtest <- function(object, terms = NULL, pattern = NULL, R = NULL, r = 0, vcov = NULL) {
  b <- stats::coef(object)
  V <- if (is.null(vcov)) stats::vcov(object) else stats::vcov(object, type = vcov)
  if (is.null(R)) {
    if (is.null(terms)) {
      if (is.null(pattern)) stop("give 'terms', 'pattern' or 'R'", call. = FALSE)
      terms <- grep(pattern, names(b), value = TRUE)
      if (!length(terms)) stop("no coefficient matches '", pattern, "'", call. = FALSE)
    }
    unknown <- setdiff(terms, names(b))
    if (length(unknown)) {
      stop("unknown coefficient(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    R <- matrix(0, length(terms), length(b))
    R[cbind(seq_along(terms), match(terms, names(b)))] <- 1
  } else {
    R <- as.matrix(R)
    if (ncol(R) != length(b)) stop("'R' must have one column per coefficient", call. = FALSE)
    terms <- NULL
  }
  r <- rep_len(r, nrow(R))
  diff <- as.numeric(R %*% b) - r
  middle <- R %*% V %*% t(R)
  statistic <- as.numeric(t(diff) %*% solve(middle, diff))
  df <- nrow(R)
  structure(list(statistic = statistic, df = df,
                 p.value = stats::pchisq(statistic, df, lower.tail = FALSE),
                 terms = terms,
                 vcov_type = if (is.null(vcov)) object$vcov_type %||% "classical" else vcov),
            class = "panglm_jointtest")
}

#' @export
print.panglm_jointtest <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nJoint Wald test (vcov: ", x$vcov_type, ")\n", sep = "")
  if (!is.null(x$hypothesis)) cat("H0:", x$hypothesis, "\n")
  if (!is.null(x$terms)) cat("H0: all of", paste(x$terms, collapse = ", "), "are zero\n")
  cat("Chi-squared = ", format(x$statistic, digits = digits), ", df = ", x$df,
      ", p-value = ", format.pval(x$p.value, digits = digits), "\n\n", sep = "")
  invisible(x)
}

#' Design matrix of a panglm fit
#'
#' Returns the model matrix used for the fit (`data = NULL`) or rebuilt for
#' new data with the stored terms, factor levels, contrasts and spline
#' knots (`splines::ns()`, `splines::bs()` and `poly()` keep the basis of
#' the fit). Its columns match `coef(object)`. With a regression spline,
#' `model.matrix(fit, grid) %*% coef(fit)` traces the fitted curve, and the
#' covariance `vcov(fit)` gives pointwise confidence bands; for
#' fixed-effects fits the curve is identified up to a constant, so center
#' it at a reference value.
#'
#' @param object a `"panglm"` fit
#' @param data optional new data containing the regressors
#' @param cre_means for correlated random effects fits with new data: a
#'   matrix of the `mean_*` columns, one row per row of `data`; when
#'   `NULL` they are recomputed as individual means within `data`, which
#'   then needs the individual index
#' @param ... unused
#' @return a numeric matrix with one column per coefficient
#' @examples
#' data(copd)
#' fit <- panglm(fev1 ~ splines::ns(crp, 3), data = copd,
#'               index = c("id", "visit"), model = "within")
#' grid <- data.frame(crp = seq(min(copd$crp), max(copd$crp), length.out = 5))
#' X <- model.matrix(fit, grid)
#' X0 <- model.matrix(fit, data.frame(crp = median(copd$crp)))
#' Xc <- sweep(X, 2, X0)            # centered at the median CRP
#' cbind(grid, effect = drop(Xc %*% coef(fit)),
#'       se = sqrt(rowSums((Xc %*% vcov(fit)) * Xc)))
#' @export
model.matrix.panglm <- function(object, data = NULL, cre_means = NULL, ...) {
  if (is.null(data)) {
    X <- object$X
    return(X[object$inverse_order %||% seq_len(nrow(X)), , drop = FALSE])
  }
  panglm_design(object, data, cre_means)
}
