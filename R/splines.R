#' Columns of the model matrix that belong to one variable
#'
#' Terms whose only variable is `variable` (e.g. `x`, `ns(x, 4)`,
#' `poly(x, 2)`). Terms that mix `variable` with other regressors
#' (interactions) are refused because the curve would then depend on them.
#'
#' @keywords internal
#' @noRd
variable_columns <- function(object, variable) {
  labels <- attr(object$terms, "term.labels")
  vars <- lapply(labels, function(l) all.vars(stats::reformulate(l)))
  own <- labels[vapply(vars, function(v) identical(v, variable), logical(1))]
  mixed <- labels[vapply(vars, function(v) variable %in% v && length(v) > 1L, logical(1))]
  if (length(mixed)) {
    stop("'", variable, "' enters interaction term(s) ", paste(mixed, collapse = ", "),
         "; the curve is not defined without fixing the other variables", call. = FALSE)
  }
  if (!length(own)) stop("no model term involves only '", variable, "'", call. = FALSE)
  cn <- names(object$coefficients)
  hit <- vapply(cn, function(n) {
    any(n == own | (startsWith(n, own) & grepl("^[0-9]+$", substring(n, nchar(own) + 1L))))
  }, logical(1))
  if (!any(hit)) stop("the columns of '", variable, "' were all dropped", call. = FALSE)
  cn[hit]
}

#' Basis of one variable's terms evaluated at given values
#'
#' @keywords internal
#' @noRd
variable_basis <- function(object, variable, values, columns) {
  nd <- object$data[rep(1L, length(values)), , drop = FALSE]
  rownames(nd) <- NULL
  nd[[variable]] <- values
  cre_means <- if (length(object$cre_columns)) {
    matrix(0, length(values), length(object$cre_columns))
  }
  stats::model.matrix(object, nd, cre_means = cre_means)[, columns, drop = FALSE]
}

#' Estimated effect curve of a continuous regressor
#'
#' For a regressor entered through a regression spline (`splines::ns()`,
#' `splines::bs()`), a polynomial or linearly, returns the change in the
#' linear predictor when the regressor moves from the reference value
#' `ref` to each of `values`, with delta-method pointwise confidence limits.
#' Centering at a reference is what makes the curve comparable across
#' models: in fixed-effects fits the level is absorbed by the individual
#' effects. The scale is the link scale: a difference in mean for the
#' identity link, a log rate ratio for Poisson and negative binomial models,
#' a log odds ratio for the logit link. Use `plot()` to draw it.
#'
#' @param object a `"panglm"` fit
#' @param variable name of the regressor
#' @param values values at which to evaluate the curve; default: 100 points
#'   between the 1st and 99th percentile of the regressor in the fit
#' @param ref reference value; default: the median of the regressor
#' @param vcov `NULL` (the covariance stored on the fit) or a type passed to
#'   [vcov.panglm()]
#' @param level confidence level of the pointwise band
#' @return a data frame of class `"panglm_effectcurve"` with columns
#'   `value`, `estimate`, `se`, `lower` and `upper`
#' @examples
#' data(copd)
#' fit <- panglm(fev1 ~ splines::ns(crp, 3), data = copd,
#'               index = c("id", "visit"), model = "within", vcov = "cluster")
#' curve_crp <- effectcurve(fit, "crp")
#' head(curve_crp)
#' plot(curve_crp)
#' @export
effectcurve <- function(object, variable, values = NULL, ref = NULL, vcov = NULL, level = 0.95) {
  x <- object$data[[variable]]
  if (is.null(x) || !is.numeric(x)) stop("'", variable, "' is not a numeric regressor of the fit",
                                         call. = FALSE)
  columns <- variable_columns(object, variable)
  values <- values %||% seq(stats::quantile(x, 0.01, na.rm = TRUE),
                            stats::quantile(x, 0.99, na.rm = TRUE), length.out = 100L)
  ref <- ref %||% stats::median(x, na.rm = TRUE)
  B <- variable_basis(object, variable, values, columns)
  B0 <- variable_basis(object, variable, ref, columns)
  D <- sweep(B, 2, B0[1, ])
  V <- if (is.null(vcov)) stats::vcov(object) else stats::vcov(object, type = vcov)
  V <- V[columns, columns, drop = FALSE]
  estimate <- drop(D %*% object$coefficients[columns])
  se <- sqrt(pmax(rowSums((D %*% V) * D), 0))
  z <- stats::qnorm(1 - (1 - level) / 2)
  out <- data.frame(value = values, estimate = estimate, se = se,
                    lower = estimate - z * se, upper = estimate + z * se)
  attr(out, "variable") <- variable
  attr(out, "ref") <- ref
  attr(out, "level") <- level
  class(out) <- c("panglm_effectcurve", "data.frame")
  out
}

#' Plot an effect curve
#'
#' @param x a `"panglm_effectcurve"` object from [effectcurve()]
#' @param xlab,ylab axis labels
#' @param band_col colour of the confidence band
#' @param ... passed to [graphics::plot()]
#' @return `x`, invisibly
#' @export
plot.panglm_effectcurve <- function(x, xlab = attr(x, "variable"),
                                    ylab = "Effect on the linear predictor",
                                    band_col = "grey85", ...) {
  graphics::plot(x$value, x$estimate, type = "n", xlab = xlab, ylab = ylab,
                 ylim = range(x$lower, x$upper), ...)
  graphics::polygon(c(x$value, rev(x$value)), c(x$lower, rev(x$upper)),
                    col = band_col, border = NA)
  graphics::lines(x$value, x$estimate, lwd = 2)
  graphics::abline(h = 0, lty = 3)
  graphics::abline(v = attr(x, "ref"), lty = 3)
  invisible(x)
}

#' Test whether the effect of a regressor is linear
#'
#' For a regressor entered through a spline or polynomial basis whose span
#' contains the linear function (true for `splines::ns()`, `splines::bs()`
#' and `poly()`), linearity means the basis coefficients are proportional
#' to the coefficients `a` that reproduce `x` itself (`x = c + B a`). The
#' test is the joint Wald test of the `q - 1` restrictions that the
#' coefficient vector is orthogonal to the complement of `a`, so it accepts
#' any covariance type, including cluster-robust ones. For a test of whether
#' the regressor matters at all, use [jointtest()] on its columns.
#'
#' @inheritParams effectcurve
#' @return an object of class `"panglm_jointtest"`
#' @examples
#' data(copd)
#' fit <- panglm(fev1 ~ splines::ns(crp, 3), data = copd,
#'               index = c("id", "visit"), model = "within", vcov = "cluster")
#' lineartest(fit, "crp")
#' @export
lineartest <- function(object, variable, vcov = NULL) {
  x <- object$data[[variable]]
  if (is.null(x) || !is.numeric(x)) stop("'", variable, "' is not a numeric regressor of the fit",
                                         call. = FALSE)
  columns <- variable_columns(object, variable)
  if (length(columns) < 2L) stop("'", variable, "' enters with a single column: nothing to test",
                                 call. = FALSE)
  B <- object$X[, columns, drop = FALSE]
  ls <- stats::lm.fit(cbind(1, B), x)
  if (max(abs(ls$residuals)) > 1e-6 * max(1, abs(x))) {
    stop("the basis of '", variable, "' does not contain the linear function", call. = FALSE)
  }
  a <- ls$coefficients[-1L]
  Q <- qr.Q(qr(matrix(a)), complete = TRUE)[, -1L, drop = FALSE]
  R <- matrix(0, ncol(Q), length(object$coefficients))
  R[, match(columns, names(object$coefficients))] <- t(Q)
  out <- jointtest(object, R = R, vcov = vcov)
  out$terms <- NULL
  out$hypothesis <- paste0("the effect of ", variable, " is linear")
  out
}
