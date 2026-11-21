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
  # only data columns count: a symbol such as the knots in ns(x, knots = K) is not a regressor
  vars <- lapply(labels, function(l) intersect(all.vars(stats::reformulate(l)), names(object$data)))
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

#' Estimated effect curve of a regressor
#'
#' Three views of how the expected outcome depends on one regressor, each
#' with delta-method pointwise confidence limits from the covariance of the
#' fit. Use `plot()` to draw them.
#'
#' * `scale = "link"` (default): for a regressor entered through a
#'   regression spline (`splines::ns()`, `splines::bs()`), a polynomial or
#'   linearly, the change in the linear predictor when the regressor moves
#'   from the reference value `ref` to each of `values`. Centering at a
#'   reference is what makes the curve comparable across models: in
#'   fixed-effects fits the level is absorbed by the individual effects. The
#'   scale is the link scale: a difference in mean for the identity link, a
#'   log rate ratio for Poisson and negative binomial models, a log odds
#'   ratio for the logit link.
#' * `scale = "response"`: the average adjusted prediction (partial
#'   dependence), that is the mean outcome when the regressor is set to each
#'   value for every observation of the fitted data, all other regressors at
#'   their observed values. The mean is the one used by [partialeffects()]:
#'   the GLM mean for pooled fits and the population-averaged mean for
#'   random-effects fits, with the individual means of correlated random
#'   effects fits held at their fitted values. For a factor, logical or 0/1
#'   regressor the curve has one point per level. Not available for
#'   fixed-effects (within) fits, where the level of the mean is not
#'   identified.
#' * `scale = "slope"`: the average marginal effect `d E[y] / dx` when the
#'   regressor is set to each value for every observation, the derivative of
#'   the response-scale curve (a central difference in the regressor).
#'   Averaged over the sample at the observed values it would be the APE of
#'   [partialeffects()].
#'
#' On the response and slope scales the regressor may enter interactions,
#' and `by` gives one curve per level of a second variable, set to that
#' level for every observation. Standard errors treat variance components,
#' dispersion and fixed-effect intercepts as fixed, as in
#' [partialeffects()].
#'
#' With `by` and `contrast = TRUE` the curve is the difference between each
#' level of `by` and its first level, with the delta-method band of the
#' difference itself (not two overlapping bands): for example the effect of
#' smoking on birth weight as a function of maternal age in a model with a
#' spline-by-smoking interaction. On the link scale this is the average
#' difference in the linear predictor (a log odds ratio or log rate ratio
#' curve); on the response and slope scales the difference of the average
#' predictions or slopes. `at` fixes other regressors at given values for
#' every observation, so that with all of them fixed the curve is the
#' prediction for one covariate profile (a typical patient).
#'
#' On the response scale (without `contrast`) the band is computed on the
#' link scale of the family and back-transformed, `g^-1(g(m) -/+ z se_g)`
#' with `se_g = se / g'^-1(g(m))`, so that it stays inside the range of the
#' mean (a probability band inside \[0, 1\], a rate band above 0); `se` is
#' still the delta-method standard error on the response scale. For a
#' single covariate profile of a pooled fit this is the usual
#' `predict(type = "link", se.fit = TRUE)` band. Identity-link and Tobit
#' fits keep the symmetric Wald band.
#'
#' For random-effects fits the response and slope curves are
#' population-averaged by default (`type = "marginal"`): the mean is
#' integrated over the random effects, as in [partialeffects()].
#' `type = "conditional"` gives instead the curve of a typical individual,
#' whose random effects are zero (`b_i = 0`): the inverse link of the fixed
#' part of the linear predictor (the Tobit censored mean with the
#' within-individual standard deviation). For a random-intercept logit fit
#' at a single covariate profile it is `plogis(x'beta)`. Pooled and
#' fixed-effects fits have no random effects, so both types agree.
#'
#' @param object a `"panglm"` fit
#' @param variable name of the regressor
#' @param values values at which to evaluate the curve; default: 100 points
#'   between the 1st and 99th percentile of the regressor in the fit, or
#'   every level of a discrete regressor
#' @param ref reference value of the link-scale curve; default: the median
#'   of the regressor
#' @param vcov `NULL` (the covariance stored on the fit) or a type passed to
#'   [vcov.panglm()]
#' @param level confidence level of the pointwise band
#' @param scale `"link"`, `"response"` or `"slope"`, see Details
#' @param by optional name of a discrete variable (factor, character,
#'   logical, or at most 10 distinct values) giving one curve per level;
#'   response and slope scales, or every scale with `contrast = TRUE`
#' @param contrast if `TRUE` (with `by`), the difference curve of each
#'   level of `by` minus the first level, with its own confidence band
#' @param at optional named list of other regressors, each set to the given
#'   value for every observation (response, slope and contrast curves)
#' @param type `"marginal"` (default: population-averaged) or
#'   `"conditional"` (a typical individual, random effects at zero) for
#'   response and slope curves of random-effects fits, see Details
#' @return a data frame of class `"panglm_effectcurve"` with columns
#'   `value`, `estimate`, `se`, `lower` and `upper` (and `by` when `by` is
#'   given; with `contrast = TRUE` its levels read `"level - first"`). The
#'   attributes record the variable, scale, reference value and
#'   the percentiles of the observed regressor used for the rug.
#' @family panglm methods
#' @section Validation:
#' The package tests compare the response and slope curves of pooled fits
#' with `marginaleffects::avg_predictions()` and
#' `marginaleffects::avg_slopes()` on counterfactual grids.
#' @examples
#' data(copd)
#' fit <- panglm(fev1 ~ splines::ns(crp, 3), data = copd,
#'               index = c("id", "visit"), model = "within", vcov = "cluster")
#' curve_crp <- effectcurve(fit, "crp")
#' head(curve_crp)
#' plot(curve_crp)
#'
#' pois <- panglm(exacerbations ~ splines::ns(crp, 3) + factor(smoker),
#'                data = copd, index = c("id", "visit"), model = "random",
#'                family = "poisson", cre = TRUE, vcov = "cluster")
#' plot(effectcurve(pois, "crp", scale = "response", by = "smoker"))
#' plot(effectcurve(pois, "crp", scale = "slope"))
#' plot(effectcurve(pois, "smoker", scale = "response"))
#' @export
effectcurve <- function(object, variable, values = NULL, ref = NULL, vcov = NULL, level = 0.95,
                        scale = c("link", "response", "slope"), by = NULL, contrast = FALSE,
                        at = NULL, type = c("marginal", "conditional")) {
  scale <- match.arg(scale)
  type <- match.arg(type)
  if (!inherits(object, "panglm")) stop("'object' must be a panglm fit", call. = FALSE)
  if (contrast && is.null(by)) stop("'contrast' needs 'by'", call. = FALSE)
  if (scale != "link" || contrast) {
    return(response_curve(object, variable, values, vcov, level, scale, by, contrast, at,
                          type))
  }
  if (!is.null(by)) {
    stop("'by' needs scale = \"response\" or \"slope\", or contrast = TRUE; on the link scale ",
         "the curve of a variable without interactions does not depend on the other regressors",
         call. = FALSE)
  }
  if (length(at)) stop("'at' needs scale = \"response\" or \"slope\", or contrast = TRUE",
                       call. = FALSE)
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
  attr(out, "scale") <- "link"
  attr(out, "discrete") <- FALSE
  attr(out, "rug") <- curve_rug(x)
  class(out) <- c("panglm_effectcurve", "data.frame")
  out
}

#' Plot an effect curve
#'
#' Draws the estimate as a line over a shaded pointwise confidence band, a
#' rug of the percentiles of the observed regressor, and dotted reference
#' lines: at zero and at the reference value for link-scale curves, at zero
#' for slope curves. Curves by group share one panel, one colour each, with
#' a legend. Curves of a discrete regressor are drawn as points with
#' confidence whiskers.
#'
#' @param x a `"panglm_effectcurve"` object from [effectcurve()]
#' @param xlab,ylab axis labels; `ylab` defaults to a label for the scale
#' @param band_col colour of the confidence band (single curve)
#' @param rug if `TRUE`, mark the percentiles of the observed regressor on
#'   the horizontal axis
#' @param col colour(s) of the curve(s); black for a single curve and the
#'   Okabe-Ito palette for curves by group by default
#' @param legend position of the legend for curves by group: `"auto"` (the
#'   corner of the plot least covered by the curves and bands), a keyword
#'   accepted by [graphics::legend()], or `FALSE` for none
#' @param main plot title; `NULL` (default) for none
#' @param mar margins in lines; by default compact ones (`c(4, 4, 1, 1)`,
#'   the top larger with a title) when the current margins are R's
#'   defaults, the current ones otherwise. Margins set here are restored on
#'   exit, so to add layers after the call set `par(mar = )` beforehand.
#' @param bylabel legend title for curves by group; the name of the `by`
#'   column by default, `NULL` for none
#' @param ... passed to [graphics::plot()]
#' @return `x`, invisibly
#' @family panglm methods
#' @examples
#' data(copd)
#' fit <- panglm(fev1 ~ splines::ns(crp, 3), data = copd,
#'               index = c("id", "visit"), model = "within")
#' plot(effectcurve(fit, "crp"), band_col = "lightblue")
#' @export
plot.panglm_effectcurve <- function(x, xlab = attr(x, "variable"), ylab = NULL,
                                    band_col = "grey85", rug = TRUE, col = NULL,
                                    legend = "auto", main = NULL, mar = NULL,
                                    bylabel = attr(x, "by"), ...) {
  scale <- attr(x, "scale") %||% "link"
  variable <- attr(x, "variable")
  contrast <- isTRUE(attr(x, "contrast"))
  ylab <- ylab %||% if (contrast) {
    switch(scale, link = "Difference in the linear predictor",
           response = "Difference in average predicted outcome",
           slope = paste("Difference in average marginal effect of", variable))
  } else if (identical(attr(x, "type"), "conditional")) {
    switch(scale, response = "Predicted outcome (random effects at 0)",
           slope = paste("Marginal effect of", variable, "(random effects at 0)"))
  } else switch(scale,
    link = "Effect on the linear predictor",
    response = "Average predicted outcome",
    slope = paste("Average marginal effect of", variable)
  )
  op <- panglm_compact_mar(mar, if (is.null(main)) 1 else 2.5)
  on.exit(if (!is.null(op)) graphics::par(op), add = TRUE)
  # a single difference curve is drawn like a single curve
  groups <- if (is.null(x$by) || (contrast && nlevels(x$by) == 1L)) NULL else levels(x$by)
  k <- max(1L, length(groups))
  col <- col %||% if (is.null(groups)) "black" else curve_palette(k)
  col <- rep_len(col, k)
  group_rows <- if (is.null(groups)) list(seq_len(nrow(x))) else split(seq_len(nrow(x)), x$by)
  ylim <- range(x$estimate, x$lower, x$upper, finite = TRUE)

  if (isTRUE(attr(x, "discrete"))) {
    lev <- levels(x$value)
    offset <- if (k > 1L) seq(-0.15, 0.15, length.out = k) else 0
    graphics::plot(c(0.5, length(lev) + 0.5), ylim, type = "n", xaxt = "n",
                   xlab = xlab, ylab = ylab, main = main, ...)
    graphics::axis(1, at = seq_along(lev), labels = lev)
    for (g in seq_len(k)) {
      r <- group_rows[[g]]
      at <- match(as.character(x$value[r]), lev) + offset[g]
      graphics::segments(at, x$lower[r], at, x$upper[r], col = col[g], lwd = 1.5)
      graphics::points(at, x$estimate[r], pch = 19, col = col[g])
    }
  } else {
    graphics::plot(range(x$value), ylim, type = "n", xlab = xlab, ylab = ylab, main = main,
                   ...)
    fill <- if (is.null(groups)) band_col else curve_translucent(col)
    for (g in seq_len(k)) {
      r <- group_rows[[g]]
      graphics::polygon(c(x$value[r], rev(x$value[r])), c(x$lower[r], rev(x$upper[r])),
                        col = fill[g], border = NA)
    }
    for (g in seq_len(k)) {
      r <- group_rows[[g]]
      graphics::lines(x$value[r], x$estimate[r], lwd = 2, col = col[g])
    }
    marks <- attr(x, "rug")
    marks <- marks[marks >= min(x$value) & marks <= max(x$value)]
    if (rug && length(marks)) graphics::rug(marks, ticksize = 0.02, col = "grey40")
    if (scale == "link" && !contrast) graphics::abline(v = attr(x, "ref"), lty = 3)
  }
  if (scale != "response" || contrast) graphics::abline(h = 0, lty = 3)
  if (!is.null(groups) && !isFALSE(legend)) {
    discrete <- isTRUE(attr(x, "discrete"))
    if (identical(legend, "auto")) {
      at <- if (discrete) as.integer(x$value) else x$value
      legend <- curve_legend_corner(rep(at, 3L), c(x$estimate, x$lower, x$upper))
    }
    graphics::legend(legend, legend = groups, col = col, bty = "n",
                     lwd = if (discrete) 1.5 else 2, pch = if (discrete) 19 else NA,
                     title = bylabel)
  }
  invisible(x)
}

# Compact single-panel margins when the user kept R's defaults; returns the old par (or NULL)
panglm_compact_mar <- function(mar = NULL, top = 1) {
  if (is.null(mar)) {
    if (!isTRUE(all.equal(graphics::par("mar"), c(5.1, 4.1, 4.1, 2.1)))) return(NULL)
    mar <- c(4, 4, top, 1)
  }
  graphics::par(mar = mar)
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
#' @return an object of class `"panglm_jointtest"`, as returned by
#'   [jointtest()].
#' @family tests
#' @section Validation:
#' The package tests check that the test rejects a nonlinear curve and keeps
#' its size under linearity, and the fixed-effects likelihood-ratio statistic
#' against the Wald statistic.
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
