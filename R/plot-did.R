#' Legend text widths: per entry, padded, where R supports it (R >= 4.1)
#'
#' @keywords internal
#' @noRd
did_textwidth <- function(lab, cex) {
  if (getRversion() < "4.1.0") return(NULL)
  graphics::strwidth(lab, cex = cex) + graphics::strwidth("mm", cex = cex)
}

#' Standard table for the event-study drawer
#'
#' Returns a data frame with `level`, `estimate`, `lower`, `upper`, `ref`
#' (the omitted reference period) and, when available, `band_lower` and
#' `band_upper`, restricted to `min_e <= level <= max_e`.
#'
#' @keywords internal
#' @noRd
did_plottab <- function(tab, ref = NULL, min_e = -Inf, max_e = Inf) {
  out <- data.frame(level = tab$level, estimate = tab$estimate,
                    lower = tab$lower, upper = tab$upper)
  if (!is.null(tab$band_lower)) {
    out$band_lower <- tab$band_lower
    out$band_upper <- tab$band_upper
  }
  out$ref <- if (is.null(ref)) rep(FALSE, nrow(out)) else out$level %in% ref
  out <- out[out$level >= min_e & out$level <= max_e, , drop = FALSE]
  if (!nrow(out)) stop("no estimate in the requested window", call. = FALSE)
  rownames(out) <- NULL
  out
}

#' Check a min_e / max_e pair
#'
#' @keywords internal
#' @noRd
did_check_window <- function(min_e, max_e) {
  ok <- function(v) is.numeric(v) && length(v) == 1L && !is.na(v)
  if (!ok(min_e) || !ok(max_e) || min_e > max_e) {
    stop("'min_e' and 'max_e' must be numbers with min_e <= max_e", call. = FALSE)
  }
  invisible(TRUE)
}

#' Shared drawer for event-study and aggregated DiD plots
#'
#' One style for [plot.panglm_attgt()], [plot.panglm_etwfe()] and
#' [plot.panglm_eventstudy()]: hollow points before treatment, filled
#' points after, an open diamond at the omitted reference period, black
#' pointwise intervals, uniform bands as light grey bars behind them, a
#' horizontal zero line and, for event time, a dashed line at treatment
#' onset.
#'
#' @param tab a table from did_plottab()
#' @param event logical: the levels are event times
#' @param xaxis logical: draw the x axis at the levels
#' @return `tab`, invisibly, with attribute `band` (were uniform bands drawn?)
#' @keywords internal
#' @noRd
did_drawes <- function(tab, event = TRUE, xlab = NULL, ylab = NULL, main = NULL,
                       ylim = NULL, legend = FALSE, level = 0.95, xaxis = TRUE, ...) {
  if (!is.logical(legend) || length(legend) != 1L || is.na(legend)) {
    stop("'legend' must be TRUE or FALSE", call. = FALSE)
  }
  ok <- is.finite(tab$estimate)
  if (!any(ok)) stop("no finite estimate to plot", call. = FALSE)
  x <- tab$level
  est <- tab$estimate
  lo <- ifelse(is.finite(tab$lower), tab$lower, est)
  hi <- ifelse(is.finite(tab$upper), tab$upper, est)
  band <- !is.null(tab$band_lower) && any(is.finite(tab$band_lower[ok]))
  if (band) {
    blo <- ifelse(is.finite(tab$band_lower), tab$band_lower, est)
    bhi <- ifelse(is.finite(tab$band_upper), tab$band_upper, est)
  }
  ref <- tab$ref & ok
  pre <- event & x < 0 & !ref
  if (is.null(ylim)) {
    ylim <- range(lo[ok], hi[ok], est[ok], 0, if (band) c(blo[ok], bhi[ok]))
    if (legend) ylim[2] <- ylim[2] + 0.18 * diff(ylim)
  }
  dx <- if (length(unique(x)) > 1L) min(diff(sort(unique(x)))) else 1
  dots <- list(...)
  if (is.null(dots$xlim)) dots$xlim <- range(x[ok]) + c(-0.5, 0.5) * dx
  do.call(graphics::plot, c(list(x[ok], est[ok], type = "n", xlab = xlab, ylab = ylab,
                                 main = main, ylim = ylim, xaxt = "n"), dots))
  if (xaxis) graphics::axis(1, at = sort(unique(x)))
  graphics::abline(h = 0, col = "grey40", lty = 3)
  if (event) graphics::abline(v = -0.5, lty = 2, col = "grey40")
  if (band) {
    b <- ok & !ref
    graphics::rect(x[b] - 0.2 * dx, blo[b], x[b] + 0.2 * dx, bhi[b], col = "grey85",
                   border = NA)
  }
  ci <- ok & !ref
  graphics::segments(x[ci], lo[ci], x[ci], hi[ci], lwd = 1.5)
  graphics::points(x[ci], est[ci], pch = ifelse(pre[ci], 21, 19), bg = "white", cex = 1.1)
  if (any(ref)) graphics::points(x[ref], est[ref], pch = 23, bg = "white", cex = 1.1)
  if (legend) {
    pct <- paste0(format(100 * level), "% ")
    lab <- c(if (event && any(pre)) "Pre-treatment", if (event) "Post-treatment" else "Estimate",
             if (any(ref)) "Reference", paste0(pct, "pointwise CI"),
             if (band) paste0(pct, "uniform band"))
    pch <- c(if (event && any(pre)) 21, 19, if (any(ref)) 23, NA, if (band) 15)
    lty <- c(if (event && any(pre)) 0, 0, if (any(ref)) 0, 1, if (band) 0)
    col <- c(rep("black", length(lab) - band), if (band) "grey85")
    graphics::legend("top", legend = lab, pch = pch, lty = lty, lwd = 1.5, col = col,
                     pt.bg = "white", pt.cex = ifelse(seq_along(lab) == length(lab) & band,
                                                      2, 1.1),
                     horiz = TRUE, bty = "n", cex = 0.75, seg.len = 1.2,
                     x.intersp = 0.6, text.width = did_textwidth(lab, 0.75))
  }
  attr(tab, "band") <- band
  invisible(tab)
}

#' Event-study table of an eventstudy() fit
#'
#' @keywords internal
#' @noRd
es_plottab <- function(x, min_e = -Inf, max_e = Inf) {
  tab <- x$estimates
  names(tab)[names(tab) == "event"] <- "level"
  did_plottab(tab, ref = x$ref, min_e = min_e, max_e = max_e)
}

#' Plot table of an etwfe() summary; the reference period has a zero
#' effect and a zero standard error
#'
#' @keywords internal
#' @noRd
etwfe_plottab <- function(x, min_e = -Inf, max_e = Inf) {
  tab <- x$estimates
  ref <- if (x$type == "event") tab$level[tab$estimate == 0 & tab$se == 0] else NULL
  did_plottab(tab, ref = ref, min_e = min_e, max_e = max_e)
}
