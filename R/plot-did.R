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

#' Plot the treatment timing of a staggered design
#'
#' Draws one row per individual and one column per period, sorted by
#' treatment cohort (earliest adopters at the top, never-treated units at
#' the bottom), with treated cells shaded. The right axis gives each cohort
#' and its size. With many individuals the unit labels are dropped and the
#' tiles are drawn as a raster, so the plot stays readable for thousands of
#' units. Periods in which a unit is not observed are left blank.
#'
#' @param data a data frame in long format
#' @param index `c(individual, time)`: column names, the time numeric
#' @param treat name of the column with the first treated period, as in
#'   [attgt()] and [etwfe()] (`0`, `NA` or `Inf` for never-treated units)
#' @param labels logical: show the individual labels on the left axis. The
#'   default shows them for at most 60 individuals.
#' @param col colours of the treated, not-yet-treated and never-treated
#'   cells
#' @param legend logical: add a legend below the plot
#' @param xlab,ylab,main axis labels and title
#' @param ... passed to [graphics::image()]
#' @return a data frame with one row per cohort (`cohort`, `n`; `NA` for
#'   the never treated), invisibly
#' @family policy evaluation
#' @examples
#' d <- expand.grid(id = 1:30, year = 2010:2019)
#' d$first <- rep(c(2013, 2015, 2017, 0, 0, 0), each = 5)[d$id]
#' treatmentplot(d, index = c("id", "year"), treat = "first")
#' @export
treatmentplot <- function(data, index, treat, labels = NULL,
                          col = c("#08519C", "#C6DBEF", "grey88"), legend = TRUE,
                          xlab = "Period", ylab = NULL, main = NULL, ...) {
  data <- as.data.frame(data)
  timing <- did_timing(data, index, treat)
  if (length(col) != 3L) stop("'col' must give three colours", call. = FALSE)
  id <- data[[index[1]]]
  ok <- !is.na(id) & is.finite(timing$time)
  id <- id[ok]
  time <- timing$time[ok]
  first <- timing$first[ok]
  uid <- unique(id)
  uid <- uid[order(uid)]
  ufirst <- first[match(uid, id)]
  ord <- order(is.na(ufirst), ufirst, seq_along(uid))
  uid <- uid[ord]
  ufirst <- ufirst[ord]
  n <- length(uid)
  if (is.null(labels)) labels <- n <= 60L
  periods <- sort(unique(time))
  # code 1 treated, 2 not yet treated, 3 never treated; NA not observed
  z <- matrix(NA_real_, length(periods), n)
  r <- match(id, uid)
  code <- ifelse(is.na(first), 3, ifelse(time >= first, 1, 2))
  z[cbind(match(time, periods), n + 1L - r)] <- code
  cohorts <- unique(ufirst)
  sizes <- vapply(cohorts, function(g) sum(if (is.na(g)) is.na(ufirst) else ufirst %in% g),
                  numeric(1))
  out <- data.frame(cohort = cohorts, n = sizes)
  coh_lab <- paste0(ifelse(is.na(cohorts), "Never", format(cohorts)), " (n = ", sizes, ")")
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  lab_w <- if (labels) max(graphics::strwidth(format(uid), units = "inches", cex = 0.6)) else 0
  coh_w <- max(graphics::strwidth(coh_lab, units = "inches", cex = 0.8))
  csi <- graphics::par("csi")
  graphics::par(mar = c(if (legend) 5.8 else 4, lab_w / csi + if (labels) 1.5 else 1.5,
                        if (is.null(main)) 1 else 3, coh_w / csi + 1))
  graphics::image(seq_along(periods), seq_len(n), z, col = col, breaks = c(0.5, 1.5, 2.5, 3.5),
                  xlab = "", ylab = "", main = main, axes = FALSE, useRaster = n > 500L,
                  ...)
  graphics::title(xlab = xlab, line = 2.2)
  graphics::box()
  at <- unique(round(pretty(seq_along(periods), n = min(length(periods), 12L))))
  at <- at[at >= 1 & at <= length(periods)]
  if (length(periods) <= 25L) at <- seq_along(periods)
  graphics::axis(1, at = at, labels = periods[at], cex.axis = 0.8)
  if (labels) {
    graphics::axis(2, at = n + 1L - seq_len(n), labels = format(uid), las = 1, cex.axis = 0.6,
                   tick = FALSE, line = -0.6)
  }
  if (!is.null(ylab)) graphics::mtext(ylab, side = 2, line = lab_w / csi + 0.5)
  ends <- cumsum(sizes)
  mid <- n + 1 - (ends - (sizes - 1) / 2)
  if (length(ends) > 1L) {
    graphics::abline(h = n + 0.5 - ends[-length(ends)], col = "white", lwd = 1.5)
  }
  # cohort labels on the right, pushed apart when small cohorts sit close
  usr <- graphics::par("usr")
  h <- 1.3 * graphics::strheight("M", cex = 0.8)
  k <- length(mid)
  if (k > 1L) for (i in 2:k) mid[i] <- min(mid[i], mid[i - 1L] - h)
  mid[k] <- max(mid[k], usr[3] + h / 2)
  if (k > 1L) for (i in (k - 1L):1) mid[i] <- max(mid[i], mid[i + 1L] + h)
  graphics::text(usr[2] + graphics::strwidth("m", cex = 0.8), mid, coh_lab, adj = 0,
                 cex = 0.8, xpd = NA)
  if (legend) {
    graphics::par(xpd = NA)
    graphics::legend(x = mean(graphics::par("usr")[1:2]),
                     y = graphics::grconvertY(graphics::grconvertY(0, "npc", "inches") -
                                                4.2 * csi, "inches", "user"),
                     legend = c("Treated", "Not yet treated", "Never treated"),
                     fill = col, border = NA, horiz = TRUE, bty = "n", cex = 0.8,
                     xjust = 0.5, yjust = 0.5,
                     text.width = did_textwidth(c("Treated", "Not yet treated", "Never treated"),
                                                0.8))
  }
  invisible(out)
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
