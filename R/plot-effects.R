#' Plot average partial effects
#'
#' Dot-and-whisker plot of the average partial effects returned by
#' [partialeffects()]: one row per continuous or binary regressor and per
#' non-reference factor level, the estimate as a dot, its Wald confidence
#' interval as a horizontal segment and a dotted reference line at zero. The
#' left margin is widened to fit the row labels and restored on exit. For
#' two-part models the part of the mean (`"total"`, `"binary"` or
#' `"positive"`) is shown in the axis label; plot one object per part to
#' compare them.
#'
#' @param x a `"panglm_ape"` object from [partialeffects()]
#' @param labels row labels; by default the variable name for continuous
#'   regressors and `"variable: contrast"` for discrete ones
#' @param xlab horizontal-axis label
#' @param main plot title
#' @param col colour of the dots and whiskers
#' @param pch plotting symbol of the estimates
#' @param ... passed to [graphics::plot.default()]
#' @return `x`, invisibly
#' @family panglm methods
#' @examples
#' data(copd)
#' fit <- panglm(hospitalized ~ treatment + age + crp + factor(smoker),
#'               data = copd, index = c("id", "visit"), family = "probit")
#' plot(partialeffects(fit))
#' @export
plot.panglm_ape <- function(x, labels = NULL, xlab = NULL, main = NULL, col = "black",
                            pch = 19, ...) {
  n <- nrow(x)
  if (!n) stop("nothing to plot", call. = FALSE)
  labels <- labels %||% ifelse(x$contrast == "dY/dX", x$term,
                               paste0(x$term, ": ", x$contrast))
  if (length(labels) != n) stop("'labels' must have one entry per row of 'x'", call. = FALSE)
  part <- attr(x, "part")
  xlab <- xlab %||% if (is.null(part)) {
    "Average partial effect"
  } else {
    paste0("Average partial effect (", part, " mean)")
  }
  old_par <- graphics::par(no.readonly = TRUE)["mai"]
  on.exit(graphics::par(old_par), add = TRUE)
  mai <- old_par$mai
  mai[2L] <- max(graphics::strwidth(labels, units = "inches")) + 0.3
  graphics::par(mai = mai)
  y <- rev(seq_len(n))
  xlim <- range(x$estimate, x$conf.low, x$conf.high, 0, finite = TRUE)
  graphics::plot.default(x$estimate, y, type = "n", xlim = xlim, ylim = c(0.5, n + 0.5),
                         yaxt = "n", xlab = xlab, ylab = "", main = main, ...)
  graphics::abline(v = 0, lty = 3)
  graphics::segments(x$conf.low, y, x$conf.high, y, col = col, lwd = 1.5)
  graphics::points(x$estimate, y, pch = pch, col = col)
  graphics::axis(2, at = y, labels = labels, las = 1, tick = FALSE)
  invisible(x)
}
