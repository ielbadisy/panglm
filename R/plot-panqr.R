#' Plot a panel quantile regression fit
#'
#' Base-graphics plot of the coefficients of a [panqr()] fit across the
#' estimated quantiles, one panel per term: the estimates joined by a line,
#' a shaded pointwise Wald confidence band (from the same variances as
#' `confint()`), and a dashed reference line at zero. With a single quantile
#' each panel shows the estimate and its interval. Without a covariance
#' (`vcov = "none"`) only the estimates are drawn.
#'
#' @param x a `"panqr"` fit
#' @param parm names or positions of the terms to plot; all by default
#' @param level confidence level of the bands
#' @param ncol number of panel columns; by default up to three
#' @param col colour of the estimates
#' @param band colour of the confidence bands
#' @param main optional titles, recycled over the terms; the term names by
#'   default
#' @param xlab horizontal-axis label
#' @param ylab vertical-axis label
#' @param ... additional arguments passed to [graphics::plot.default()]
#' @return `x`, invisibly
#' @examples
#' set.seed(2)
#' d <- data.frame(id = rep(1:100, each = 5), t = rep(1:5, 100))
#' d$x <- rnorm(500) + rep(rnorm(100), each = 5)
#' d$y <- d$x + rep(rnorm(100), each = 5) + (1 + 0.2 * abs(d$x)) * rnorm(500)
#' fit <- panqr(y ~ x, data = d, index = c("id", "t"),
#'              tau = c(0.1, 0.25, 0.5, 0.75, 0.9), model = "mmqr")
#' plot(fit)
#' @export
plot.panqr <- function(x, parm = NULL, level = 0.95, ncol = NULL, col = "black",
                       band = "grey85", main = NULL, xlab = expression(tau),
                       ylab = "Estimate", ...) {
  tidy_check_level(level)
  B <- x$coefficients
  terms <- rownames(B)
  if (is.null(parm)) {
    parm <- terms
  } else if (is.numeric(parm)) {
    parm <- terms[parm]
  }
  if (!length(parm) || anyNA(parm) || !all(parm %in% terms)) {
    stop("'parm' must name or index terms of the fit", call. = FALSE)
  }
  tau <- x$tau
  crit <- stats::qnorm(1 - (1 - level) / 2)
  se <- if (is.null(x$vcov)) {
    NULL
  } else {
    vapply(x$vcov, function(V) sqrt(diag(V))[terms], numeric(length(terms)))
  }
  if (!is.null(se)) se <- matrix(se, nrow = length(terms), dimnames = list(terms, NULL))
  main <- rep_len(if (is.null(main)) parm else main, length(parm))

  if (length(parm) > 1L) {
    ncol <- ncol %||% min(3L, length(parm))
    old_par <- graphics::par(mfrow = c(ceiling(length(parm) / ncol), ncol))
    on.exit(graphics::par(old_par), add = TRUE)
  }
  for (k in seq_along(parm)) {
    est <- B[parm[k], ]
    lower <- upper <- NULL
    if (!is.null(se)) {
      lower <- est - crit * se[parm[k], ]
      upper <- est + crit * se[parm[k], ]
    }
    ylim <- range(c(est, lower, upper, 0), finite = TRUE)
    xlim <- if (length(tau) > 1L) range(tau) else tau + c(-0.05, 0.05)
    graphics::plot.default(tau, est, type = "n", xlim = xlim, ylim = ylim, xlab = xlab,
                           ylab = ylab, main = main[k], ...)
    if (!is.null(lower)) {
      if (length(tau) > 1L) {
        ok <- is.finite(lower) & is.finite(upper)
        graphics::polygon(c(tau[ok], rev(tau[ok])), c(lower[ok], rev(upper[ok])),
                          col = band, border = NA)
      } else {
        graphics::segments(tau, lower, tau, upper, col = col)
      }
    }
    graphics::abline(h = 0, lty = 2, col = "grey40")
    graphics::lines(tau, est, col = col)
    graphics::points(tau, est, pch = 19, col = col)
  }
  invisible(x)
}
