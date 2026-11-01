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

#' Average prediction (and its slope) when a regressor is set for everyone
#'
#' The engine behind `effectcurve(scale = "response" | "slope")`. For each
#' value, the regressor (and the `by` variable) is set to that value for
#' every observation of the fitted data, and the mean used by
#' [partialeffects()] is averaged over the sample. The gradient with respect
#' to the coefficients is analytic in the design (`mean(g'(eta) x)`), with
#' `g'` taken by a central difference on the linear index.
#'
#' @keywords internal
#' @noRd
response_curve <- function(object, variable, values, vcov, level, scale, by) {
  data <- object$data
  if (is.null(data)) {
    stop("this fit does not store its data; re-fit it with the current panglm", call. = FALSE)
  }
  if (scale == "response" && object$model == "within") {
    stop("the level of the mean is not identified in a fixed-effects (within) fit; ",
         "use scale = \"link\" for the centred curve or scale = \"slope\" for the ",
         "average marginal effect", call. = FALSE)
  }
  x <- data[[variable]]
  discrete <- curve_discrete(x)
  if (discrete && scale == "slope") {
    stop("scale = \"slope\" needs a continuous regressor; use scale = \"response\" ",
         "or partialeffects() for '", variable, "'", call. = FALSE)
  }
  if (discrete) {
    values <- values %||% curve_levels(x, object$xlevels[[variable]])
  } else {
    if (!is.numeric(x)) stop("'", variable, "' is not a numeric regressor of the fit",
                             call. = FALSE)
    values <- values %||% seq(stats::quantile(x, 0.01, na.rm = TRUE),
                              stats::quantile(x, 0.99, na.rm = TRUE), length.out = 100L)
  }
  by_values <- list(NULL)
  if (!is.null(by)) {
    if (!is.character(by) || length(by) != 1L || is.null(data[[by]])) {
      stop("'by' must name one regressor of the model", call. = FALSE)
    }
    if (identical(by, variable)) stop("'by' must differ from 'variable'", call. = FALSE)
    z <- data[[by]]
    if (!curve_discrete(z) && length(unique(z)) > 10L) {
      stop("'by' must be a factor, character, logical or have at most 10 distinct values",
           call. = FALSE)
    }
    by_values <- as.list(curve_levels(z, object$xlevels[[by]]))
  }

  theta <- object$coefficients
  V <- if (is.null(vcov)) stats::vcov(object) else stats::vcov(object, type = vcov)
  V <- V[names(theta), names(theta), drop = FALSE]
  parts <- panglm_index_predictor(object)
  base <- parts$meanfun(data)(as.numeric(parts$design(data) %*% theta) + parts$shift)
  w <- object$weights %||% rep(1, nrow(data))
  use <- is.finite(base) & w > 0
  w <- w[use] / sum(w[use])
  average <- function(d) {
    X <- parts$design(d)[use, , drop = FALSE]
    g <- parts$meanfun(d)
    eta <- as.numeric(X %*% theta) + if (length(parts$shift) > 1L) parts$shift[use] else parts$shift
    # meanfun may need the full-length data (random slopes), so evaluate on all rows
    mean_at <- function(e) {
      full <- rep(0, nrow(d))
      full[use] <- e
      g(full)[use]
    }
    h <- 1e-6 * pmax(1, abs(eta))
    dmu <- (mean_at(eta + h) - mean_at(eta - h)) / (2 * h)
    list(estimate = sum(w * mean_at(eta)), gradient = colSums(w * dmu * X))
  }
  step <- if (scale == "slope") 1e-4 * stats::sd(x, na.rm = TRUE) else 0
  if (scale == "slope" && (!is.finite(step) || step == 0)) {
    stop("'", variable, "' does not vary in the fitted data", call. = FALSE)
  }
  point <- function(v, b) {
    d <- data
    if (!is.null(b)) d[[by]] <- ape_set(data[[by]], b)
    if (scale == "response") {
      d[[variable]] <- ape_set(x, v)
      return(average(d))
    }
    d[[variable]] <- rep(v + step, nrow(d))
    up <- average(d)
    d[[variable]] <- rep(v - step, nrow(d))
    dn <- average(d)
    list(estimate = (up$estimate - dn$estimate) / (2 * step),
         gradient = (up$gradient - dn$gradient) / (2 * step))
  }

  crit <- stats::qnorm(1 - (1 - level) / 2)
  rows <- lapply(by_values, function(b) {
    res <- lapply(values, point, b = b)
    estimate <- vapply(res, `[[`, numeric(1), "estimate")
    G <- do.call(rbind, lapply(res, `[[`, "gradient"))
    se <- sqrt(pmax(rowSums((G %*% V) * G), 0))
    out <- data.frame(value = if (discrete) factor(values, levels = values) else values,
                      estimate = estimate, se = se,
                      lower = estimate - crit * se, upper = estimate + crit * se)
    if (!is.null(b)) out <- cbind(out[1L], by = as.character(b), out[-1L])
    out
  })
  out <- do.call(rbind, rows)
  if (!is.null(by)) {
    out$by <- factor(out$by, levels = vapply(by_values, as.character, character(1)))
  }
  rownames(out) <- NULL
  attr(out, "variable") <- variable
  attr(out, "level") <- level
  attr(out, "scale") <- scale
  attr(out, "by") <- by
  attr(out, "discrete") <- discrete
  attr(out, "rug") <- if (!discrete) curve_rug(x)
  class(out) <- c("panglm_effectcurve", "data.frame")
  out
}

curve_discrete <- function(x) {
  is.factor(x) || is.character(x) || is.logical(x) ||
    (is.numeric(x) && all(x[!is.na(x)] %in% c(0, 1)))
}

curve_levels <- function(x, xlevels = NULL) {
  if (is.factor(x)) return(levels(x))
  xlevels %||% sort(unique(x[!is.na(x)]))
}

# compact summary of the observed regressor for a rug: its percentiles
curve_rug <- function(x) {
  unique(as.numeric(stats::quantile(x, seq(0, 1, by = 0.01), na.rm = TRUE, names = FALSE)))
}

# Okabe-Ito colours without yellow, for curves by group
curve_palette <- function(n) {
  rep_len(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"), n)
}

# translucent version of hex colours (band fills); other colours are kept
curve_translucent <- function(col, alpha = "33") {
  ifelse(grepl("^#[0-9A-Fa-f]{6}$", col), paste0(col, alpha), col)
}

# the plot corner holding the fewest drawn points, for the legend
curve_legend_corner <- function(x, y) {
  usr <- graphics::par("usr")
  u <- (x - usr[1L]) / (usr[2L] - usr[1L])
  v <- (y - usr[3L]) / (usr[4L] - usr[3L])
  corners <- c("topleft", "topright", "bottomleft", "bottomright")
  load <- c(sum(u < 0.35 & v > 0.7), sum(u > 0.65 & v > 0.7),
            sum(u < 0.35 & v < 0.3), sum(u > 0.65 & v < 0.3))
  corners[which.min(load)]
}
