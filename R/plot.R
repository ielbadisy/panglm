#' Plot a fitted panglm model
#'
#' Base-graphics coefficient and diagnostic plots of a [panglm()] fit.
#'
#' * `"coefficients"`: estimates with Wald confidence intervals (the
#'   dot-and-whisker layout shared by every panglm coefficient plot and by
#'   [coefplot()]); the left margin is sized to the longest term label.
#' * `"residuals"`: residuals against fitted values with a
#'   [stats::lowess()] smooth. Pearson residuals by default (see
#'   [residuals.panglm()] for the variance functions); for the ordered and
#'   Tobit families, where they are not defined, response residuals are
#'   drawn with a message. For log-link families (Poisson, negative
#'   binomial, log-link gamma) the fitted values are on a logarithmic axis,
#'   which spreads the many small means of count data.
#' * `"fitted"`: observed against fitted values with the 45-degree line.
#' * `"qq"`: normal Q-Q plot of the Pearson residuals with a reference
#'   line through the quartiles; for the Gaussian and gamma families the
#'   residuals are divided by the square root of the dispersion
#'   (standardized Pearson residuals). For discrete outcomes the plot is
#'   lumpy by construction and mainly shows skewness and outliers.
#' * `"calibration"` (binomial families, logit or probit): observed
#'   proportion against mean predicted probability in `bins` quantile bins
#'   of the predicted probability, with exact (Clopper-Pearson) binomial
#'   confidence intervals and the 45-degree line.
#' * `"rootogram"` (Poisson and negative binomial families): hanging
#'   rootogram. The curve is `sqrt(E_k)`, the expected frequency of count
#'   `k` summed over observations, `E_k = sum_i P(y_i = k)`; the bars of
#'   height `sqrt(O_k)` (observed frequency) hang from it, so that bars
#'   ending above zero mark counts predicted too often and bars ending below
#'   zero counts predicted too rarely. `P(y_i = k)` is the Poisson or NB2
#'   probability at the fitted mean for pooled and fixed-effects fits (the
#'   fixed-effects means include the estimated individual effects). For
#'   random-intercept fits it is the marginal probability of a single
#'   observation implied by the mixing distribution, which is exactly the
#'   distribution whose frequencies the rootogram compares: negative
#'   binomial with size `d` for the Poisson-gamma model, beta-negative
#'   binomial for the Hausman-Hall-Griliches negative binomial model.
#'   Random-slope fits have no closed-form marginal and are not supported.
#'
#' @param x a fitted `"panglm"` object
#' @param which one of `"coefficients"`, `"residuals"`, `"fitted"`, `"qq"`,
#'   `"calibration"` or `"rootogram"`
#' @param level confidence level for coefficient and calibration intervals
#' @param parm names or positions of the terms of the coefficient plot; by
#'   default every term except `"(Intercept)"`
#' @param type residuals of the `"residuals"` view, `"pearson"` (default)
#'   or `"response"`
#' @param bins number of quantile bins of the calibration plot
#' @param max largest count of the rootogram; by default the 99th
#'   percentile of the outcome
#' @param main optional plot title
#' @param xlab optional horizontal-axis label
#' @param ... additional arguments passed to [graphics::plot.default()]
#' @return `x`, invisibly
#' @family panglm methods
#' @examples
#' data(copd)
#' fit <- panglm(fev1 ~ crp + age + treatment, data = copd,
#'               index = c("id", "visit"), model = "random",
#'               family = "gaussian")
#' plot(fit)
#' plot(fit, parm = c("crp", "treatment"))
#' plot(fit, which = "residuals")
#' plot(fit, which = "qq")
#'
#' counts <- panglm(exacerbations ~ crp + treatment, data = copd,
#'                  index = c("id", "visit"), model = "pooling",
#'                  family = "negbin")
#' plot(counts, which = "rootogram")
#'
#' any_exacerbation <- panglm(I(exacerbations > 0) ~ crp + treatment,
#'                            data = copd, index = c("id", "visit"),
#'                            model = "pooling", family = "binomial")
#' plot(any_exacerbation, which = "calibration", bins = 5)
#' @export
plot.panglm <- function(x, which = c("coefficients", "residuals", "fitted", "qq",
                                     "calibration", "rootogram"),
                        level = 0.95, parm = NULL, type = c("pearson", "response"),
                        bins = 10, max = NULL, main = NULL, xlab = NULL, ...) {
  which <- match.arg(which)
  type <- match.arg(type)

  if (which == "coefficients") {
    tab <- .panglm_coef_matrices(
      list(list(est = x$coefficients, se = .panglm_coef_se(x$coefficients, x$vcov))),
      parm, level
    )
    .panglm_dotwhisker(tab, main = main %||% "Coefficient estimates",
                       xlab = xlab %||% "Estimate", ...)
    return(invisible(x))
  }
  if (is.null(x$fitted.values) || is.null(x$y)) {
    stop("diagnostic plots require stored fitted values", call. = FALSE)
  }
  switch(which,
    residuals = .panglm_plot_residuals(x, type, main, xlab, ...),
    fitted = .panglm_plot_fitted(x, main, xlab, ...),
    qq = .panglm_plot_qq(x, main, ...),
    calibration = .panglm_plot_calibration(x, bins, level, main, xlab, ...),
    rootogram = .panglm_plot_rootogram(x, max, main, xlab, ...)
  )
  invisible(x)
}

.panglm_scatter_col <- function(n) if (n > 2000) "#00000026" else if (n > 300) "#00000059" else "#000000A6"

.panglm_plot_residuals <- function(x, type, main, xlab, ...) {
  if (type == "pearson") {
    pr <- panglm_pearson(x, fallback = TRUE)
    res <- pr$residuals
    if (pr$type == "response") type <- "response"
  } else {
    res <- x$y - x$fitted.values
  }
  fitted <- x$fitted.values
  keep <- is.finite(fitted) & is.finite(res)
  log_axis <- identical(x$family$link, "log") || x$family$family %in% c("poisson", "negbin")
  if (log_axis) keep <- keep & fitted > 0
  if (!any(keep)) stop("no finite fitted values are available to plot", call. = FALSE)
  fitted <- fitted[keep]
  res <- res[keep]
  log_axis <- log_axis && diff(range(fitted)) > 0
  graphics::plot.default(
    fitted, res, log = if (log_axis) "x" else "",
    main = main %||% "Residuals versus fitted values",
    xlab = xlab %||% if (log_axis) "Fitted values (log scale)" else "Fitted values",
    ylab = if (type == "pearson") "Pearson residuals" else "Response residuals",
    pch = 16, cex = 0.6, col = .panglm_scatter_col(length(res)), ...
  )
  graphics::abline(h = 0, lty = 2, col = "grey40")
  if (length(unique(fitted)) > 3L) {
    sm <- stats::lowess(if (log_axis) log(fitted) else fitted, res)
    graphics::lines(if (log_axis) exp(sm$x) else sm$x, sm$y, col = "#D55E00", lwd = 2)
  }
}

.panglm_plot_fitted <- function(x, main, xlab, ...) {
  keep <- is.finite(x$fitted.values) & is.finite(x$y)
  if (!any(keep)) stop("no finite fitted values are available to plot", call. = FALSE)
  fitted <- x$fitted.values[keep]
  observed <- x$y[keep]
  lim <- range(c(fitted, observed), finite = TRUE)
  graphics::plot.default(
    fitted, observed,
    main = main %||% "Observed versus fitted values",
    xlab = xlab %||% "Fitted values", ylab = "Observed", xlim = lim, ylim = lim,
    pch = 16, cex = 0.6, col = .panglm_scatter_col(length(fitted)), ...
  )
  graphics::abline(a = 0, b = 1, lty = 2, col = "grey40")
}

.panglm_plot_qq <- function(x, main, ...) {
  pr <- panglm_pearson(x, fallback = TRUE)
  res <- pr$residuals
  lab <- if (pr$type == "pearson") "Pearson residuals" else "Response residuals"
  if (pr$type == "pearson" && x$family$family %in% c("gaussian", "gamma") &&
      length(x$dispersion) == 1L && is.finite(x$dispersion) && x$dispersion > 0) {
    res <- res / sqrt(x$dispersion)
    lab <- "Standardized Pearson residuals"
  }
  res <- res[is.finite(res)]
  if (!length(res)) stop("no finite residuals are available to plot", call. = FALSE)
  q <- stats::qqnorm(res, plot.it = FALSE)
  graphics::plot.default(q$x, q$y, main = main %||% "Normal Q-Q plot",
                         xlab = "Theoretical quantiles", ylab = lab,
                         pch = 16, cex = 0.6, col = .panglm_scatter_col(length(res)), ...)
  stats::qqline(res, col = "#D55E00", lwd = 2)
}

.panglm_plot_calibration <- function(x, bins, level, main, xlab, ...) {
  if (!identical(x$family$family, "binomial")) {
    stop("the calibration plot is available for binomial (logit or probit) fits",
         call. = FALSE)
  }
  tidy_check_level(level)
  cal <- panglm_calibration(x$y, x$fitted.values, bins, level)
  graphics::plot.default(
    NA, NA, xlim = c(0, 1), ylim = c(0, 1), asp = 1,
    main = main %||% "Calibration",
    xlab = xlab %||% "Mean predicted probability", ylab = "Observed proportion", ...
  )
  graphics::abline(a = 0, b = 1, lty = 2, col = "grey40")
  graphics::lines(cal$predicted, cal$observed, col = "#0072B2", lwd = 1.5)
  graphics::segments(cal$predicted, cal$lower, cal$predicted, cal$upper, col = "#0072B2",
                     lwd = 1.5)
  graphics::points(cal$predicted, cal$observed, pch = 19, col = "#0072B2")
  graphics::rug(cal$predicted, ticksize = 0.02, col = "grey40")
}

.panglm_plot_rootogram <- function(x, max, main, xlab, ...) {
  rg <- panglm_rootogram(x, max)
  k <- rg$count
  top <- sqrt(rg$expected)
  bottom <- top - sqrt(rg$observed)
  ylim <- range(c(0, top, bottom), finite = TRUE)
  fam <- x$family$family
  graphics::plot.default(
    NA, NA, xlim = range(k) + c(-0.6, 0.6), ylim = ylim,
    main = main %||% paste0("Hanging rootogram (",
                            if (fam == "poisson") "Poisson" else "negative binomial", ")"),
    xlab = xlab %||% "Count", ylab = "sqrt(Frequency)", ...
  )
  graphics::rect(k - 0.4, bottom, k + 0.4, top, col = "grey80", border = "grey45")
  graphics::abline(h = 0, col = "grey20")
  graphics::lines(k, top, col = "#D55E00", lwd = 2)
  graphics::points(k, top, pch = 19, cex = 0.7, col = "#D55E00")
}

#' Pearson residuals of a panglm fit, in panel order
#'
#' Returns `list(residuals, type)`. With `fallback = TRUE` an undefined
#' variance gives response residuals with a message, otherwise an error.
#'
#' @keywords internal
#' @noRd
panglm_pearson <- function(object, fallback = FALSE) {
  mu <- object$fitted.values
  y <- object$y
  V <- panglm_variance(object, mu)
  if (is.null(V)) {
    fam <- object$family$family
    if (!fallback) {
      stop("Pearson residuals are not defined for family = '", fam, "'", call. = FALSE)
    }
    message("Pearson residuals are not defined for family = '", fam,
            "'; using response residuals")
    return(list(residuals = y - mu, type = "response"))
  }
  list(residuals = (y - mu) / sqrt(V), type = "pearson")
}

#' Variance function of a panglm fit at the means `mu`
#'
#' Marginal variance for random-intercept Poisson-gamma and
#' beta-negative-binomial fits; `NULL` when undefined.
#'
#' @keywords internal
#' @noRd
panglm_variance <- function(object, mu) {
  fam <- object$family$family
  random_intercept <- identical(object$model, "random") && is.null(object$random)
  switch(fam,
    gaussian = rep(1, length(mu)),
    poisson = {
      d <- object$dispersion_param
      if (random_intercept && length(d) == 1L && is.finite(d) && d > 0) mu + mu^2 / d else mu
    },
    negbin = {
      nb <- negbin_parameters(object)
      if (length(nb$theta) == 1L && is.finite(nb$theta)) {
        mu + mu^2 / nb$theta
      } else if (!is.null(nb$ab) && nb$ab[["a"]] > 2) {
        a <- nb$ab[["a"]]; b <- nb$ab[["b"]]
        lambda <- mu * (a - 1) / b
        lambda * b * (lambda + a - 1) * (a + b - 1) / ((a - 2) * (a - 1)^2)
      } else {
        mu
      }
    },
    binomial = mu * (1 - mu),
    gamma = mu^2,
    beta = {
      disp <- object$dispersion
      if (length(disp) == 1L && is.finite(disp)) mu * (1 - mu) * disp else NULL
    },
    NULL
  )
}

#' Binned calibration table
#'
#' @keywords internal
#' @noRd
panglm_calibration <- function(y, p, bins = 10, level = 0.95) {
  if (length(bins) != 1L || !is.finite(bins) || bins < 1) {
    stop("'bins' must be a positive integer", call. = FALSE)
  }
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]
  p <- p[keep]
  breaks <- unique(stats::quantile(p, seq(0, 1, length.out = round(bins) + 1L), names = FALSE))
  bin <- if (length(breaks) < 2L) rep(1L, length(p)) else {
    cut(p, breaks, include.lowest = TRUE, labels = FALSE)
  }
  n <- tabulate(bin, max(bin))
  ok <- n > 0
  events <- as.numeric(rowsum(y, bin, reorder = TRUE))
  predicted <- as.numeric(rowsum(p, bin, reorder = TRUE)) / n[ok]
  n <- n[ok]
  alpha <- 1 - level
  successes <- round(events)
  lower <- ifelse(successes > 0, stats::qbeta(alpha / 2, successes, n - successes + 1), 0)
  upper <- ifelse(successes < n, stats::qbeta(1 - alpha / 2, successes + 1, n - successes), 1)
  data.frame(bin = seq_along(n), n = n, predicted = predicted, observed = events / n,
             lower = lower, upper = upper)
}

#' Observed and expected count frequencies of a count fit
#'
#' @keywords internal
#' @noRd
panglm_rootogram <- function(object, max = NULL) {
  fam <- object$family$family
  if (!fam %in% c("poisson", "negbin")) {
    stop("the rootogram is available for Poisson and negative binomial fits", call. = FALSE)
  }
  random <- identical(object$model, "random")
  if (random && !is.null(object$random)) {
    stop("the rootogram is not available for random-slope fits: their marginal count ",
         "distribution has no closed form", call. = FALSE)
  }
  if (identical(object$model, "between")) {
    stop("the rootogram is not available for between fits", call. = FALSE)
  }
  y <- object$y
  mu <- object$fitted.values
  keep <- is.finite(y) & is.finite(mu)
  if (!is.null(object$weights)) keep <- keep & object$weights > 0
  if (!any(keep)) stop("no finite fitted values are available", call. = FALSE)
  if (any(abs(y[keep] - round(y[keep])) > 1e-8)) {
    stop("the rootogram requires an integer count outcome", call. = FALSE)
  }
  kmax <- max %||% base::max(1, ceiling(stats::quantile(y[keep], 0.99, names = FALSE)))
  if (length(kmax) != 1L || !is.finite(kmax) || kmax < 0) {
    stop("'max' must be a nonnegative number", call. = FALSE)
  }
  k <- 0:round(kmax)
  pmf <- if (fam == "poisson") {
    if (random) {
      d <- object$dispersion_param
      function(kk) stats::dnbinom(kk, size = d, mu = mu[keep])
    } else {
      function(kk) stats::dpois(kk, mu[keep])
    }
  } else if (random) {
    a <- object$negbin_a; b <- object$negbin_b
    if (is.null(a) || is.null(b)) stop("the fit stores no negative binomial parameters", call. = FALSE)
    lambda <- exp(object$linear.predictors[keep])
    function(kk) exp(lgamma(lambda + kk) - lgamma(lambda) - lgamma(kk + 1) +
                       lbeta(a + lambda, b + kk) - lbeta(a, b))
  } else {
    theta <- object$theta
    if (length(theta) != 1L || !is.finite(theta)) {
      stop("the fit stores no negative binomial theta", call. = FALSE)
    }
    function(kk) stats::dnbinom(kk, size = theta, mu = mu[keep])
  }
  expected <- vapply(k, function(kk) sum(pmf(kk)), numeric(1))
  observed <- tabulate(y[keep] + 1, length(k))
  data.frame(count = k, observed = observed, expected = expected)
}

#' Plot a fitted panglm hurdle model
#'
#' Displays coefficient estimates and Wald confidence intervals for the
#' binary hurdle component, the positive-count component, or both, in the
#' dot-and-whisker layout of [plot.panglm()].
#'
#' @param x a fitted `"panglm_hurdle"` object
#' @param which one of `"both"`, `"zero"`, or `"count"`
#' @param level confidence level for coefficient intervals
#' @param parm names or positions of the terms to plot; by default every
#'   term except `"(Intercept)"`
#' @param ... additional arguments passed to [graphics::plot.default()]
#' @return `x`, invisibly
#' @family panhurdle methods
#' @examples
#' data(copd)
#' fit <- panhurdle(exacerbations ~ crp, data = copd,
#'                      index = c("id", "visit"))
#' plot(fit)
#' @export
plot.panglm_hurdle <- function(x, which = c("both", "zero", "count"),
                               level = 0.95, parm = NULL, ...) {
  which <- match.arg(which)
  parts <- if (which == "both") c("zero", "count") else which
  titles <- c(zero = "Positive-outcome component",
              count = paste("Positive count:",
                            if (identical(x$count_family, "negbin")) "NB2" else "Poisson"))
  panels <- lapply(parts, function(part) {
    fit <- x[[part]]
    list(main = titles[[part]],
         parts = list(list(est = fit$coefficients,
                           se = .panglm_coef_se(fit$coefficients, fit$vcov))))
  })
  .panglm_coef_panels(panels, parm, level, ...)
  invisible(x)
}

#' Draw one or several coefficient panels side by side
#'
#' `panels` is a list of `list(main, parts)`; `parts` as in
#' `.panglm_coef_matrices()`.
#'
#' @keywords internal
#' @noRd
.panglm_coef_panels <- function(panels, parm, level, xlab = "Estimate", ...) {
  if (length(panels) > 1L) {
    old_par <- graphics::par(mfrow = c(1, length(panels)))
    on.exit(graphics::par(old_par), add = TRUE)
  }
  for (p in panels) {
    tab <- .panglm_coef_matrices(p$parts, parm, level, strict = length(panels) == 1L)
    .panglm_dotwhisker(tab, main = p$main, xlab = xlab, ...)
  }
  invisible(NULL)
}

#' Single-model coefficient plot with every term (internal interface)
#'
#' Kept for internal callers written against the former helper; draws with
#' the shared dot-and-whisker layout.
#'
#' @keywords internal
#' @noRd
.panglm_plot_coefficients <- function(coefficients, vcov, level, main, xlab, ...) {
  tab <- .panglm_coef_matrices(
    list(list(est = coefficients, se = .panglm_coef_se(coefficients, vcov))),
    seq_along(coefficients), level
  )
  .panglm_dotwhisker(tab, main = main, xlab = xlab, ...)
}
