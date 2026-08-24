#' Plot a fitted panglm model
#'
#' Provides base-graphics coefficient and diagnostic plots. The coefficient
#' view shows Wald confidence intervals. Diagnostic views are available when
#' the fitted object stores fitted values.
#'
#' @param x a fitted `"panglm"` object
#' @param which one of `"coefficients"`, `"residuals"`, or `"fitted"`
#' @param level confidence level for coefficient intervals
#' @param main optional plot title
#' @param xlab optional horizontal-axis label
#' @param ... additional arguments passed to [graphics::plot.default()]
#' @return `x`, invisibly
#' @examples
#' data(copd)
#' fit <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
#'               model = "within", family = "gaussian")
#' plot(fit)
#' @export
plot.panglm <- function(x, which = c("coefficients", "residuals", "fitted"),
                        level = 0.95, main = NULL, xlab = NULL, ...) {
  which <- match.arg(which)

  if (which == "coefficients") {
    .panglm_plot_coefficients(
      x$coefficients, x$vcov, level = level,
      main = main %||% "Coefficient estimates",
      xlab = xlab %||% "Estimate", ...
    )
  } else {
    if (is.null(x$fitted.values) || is.null(x$y)) {
      stop("diagnostic plots require stored fitted values", call. = FALSE)
    }
    keep <- is.finite(x$fitted.values) & is.finite(x$y)
    if (!any(keep)) stop("no finite fitted values are available to plot", call. = FALSE)
    fitted <- x$fitted.values[keep]
    observed <- x$y[keep]

    if (which == "residuals") {
      graphics::plot.default(
        fitted, observed - fitted,
        main = main %||% "Residuals versus fitted values",
        xlab = xlab %||% "Fitted values", ylab = "Residuals", ...
      )
      graphics::abline(h = 0, lty = 2, col = "grey40")
    } else {
      lim <- range(c(fitted, observed), finite = TRUE)
      graphics::plot.default(
        fitted, observed,
        main = main %||% "Observed versus fitted values",
        xlab = xlab %||% "Fitted values", ylab = "Observed", xlim = lim,
        ylim = lim, ...
      )
      graphics::abline(a = 0, b = 1, lty = 2, col = "grey40")
    }
  }

  invisible(x)
}

#' Plot a fitted panglm hurdle model
#'
#' Displays coefficient estimates and Wald confidence intervals for the
#' binary hurdle component, the positive-count component, or both.
#'
#' @param x a fitted `"panglm_hurdle"` object
#' @param which one of `"both"`, `"zero"`, or `"count"`
#' @param level confidence level for coefficient intervals
#' @param ... additional arguments passed to [graphics::plot.default()]
#' @return `x`, invisibly
#' @examples
#' data(copd)
#' fit <- panglm_hurdle(exacerbations ~ crp, data = copd,
#'                      index = c("id", "visit"))
#' plot(fit)
#' @export
plot.panglm_hurdle <- function(x, which = c("both", "zero", "count"),
                               level = 0.95, ...) {
  which <- match.arg(which)
  parts <- if (which == "both") c("zero", "count") else which
  old_par <- NULL
  if (length(parts) == 2L) {
    old_par <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(old_par), add = TRUE)
  }

  for (part in parts) {
    fit <- x[[part]]
    title <- if (part == "zero") {
      "Positive-outcome component"
    } else {
      paste("Positive count:", if (identical(x$count_family, "negbin")) "NB2" else "Poisson")
    }
    .panglm_plot_coefficients(
      fit$coefficients, fit$vcov, level = level,
      main = title, xlab = "Estimate", ...
    )
  }

  invisible(x)
}

.panglm_plot_coefficients <- function(coefficients, vcov, level, main, xlab, ...) {
  if (!length(coefficients)) stop("the model has no coefficients to plot", call. = FALSE)
  se <- sqrt(diag(vcov))
  critical <- stats::qnorm(1 - (1 - level) / 2)
  lower <- coefficients - critical * se
  upper <- coefficients + critical * se
  finite_limits <- range(c(lower, upper, coefficients), finite = TRUE)
  if (!all(is.finite(finite_limits))) finite_limits <- range(coefficients, finite = TRUE)
  if (diff(finite_limits) == 0) finite_limits <- finite_limits + c(-0.5, 0.5)

  positions <- rev(seq_along(coefficients))
  graphics::plot.default(
    coefficients, positions, xlim = finite_limits,
    ylim = c(0.5, length(coefficients) + 0.5), yaxt = "n",
    ylab = "", xlab = xlab, main = main, pch = 19, ...
  )
  graphics::segments(lower, positions, upper, positions)
  graphics::axis(2, at = positions, labels = names(coefficients), las = 1)
  graphics::abline(v = 0, lty = 2, col = "grey40")
}
