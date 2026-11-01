#' Legend text widths: per entry, padded, where R supports it (R >= 4.1)
#'
#' @keywords internal
#' @noRd
did_textwidth <- function(lab, cex) {
  if (getRversion() < "4.1.0") return(NULL)
  graphics::strwidth(lab, cex = cex) + graphics::strwidth("mm", cex = cex)
}

#' Legend layout that fits the plot region
#'
#' Chooses between one horizontal row and a few columns so that a legend at
#' the top of the plot region is never wider than the region (for example in
#' a narrow panel of `par(mfrow = c(1, 2))`). Widths are measured in inches
#' from the current device and margins, so this is called before the plot is
#' drawn. Returns `ncol`, `nrow` and `height` (inches).
#'
#' @keywords internal
#' @noRd
did_leglayout <- function(lab, cex, seg.len = 1.2, x.intersp = 1) {
  k <- length(lab)
  xc <- cex * graphics::par("cin")[1]
  yc <- cex * graphics::par("cin")[2]
  item <- graphics::strwidth(lab, units = "inches", cex = cex) +
    (seg.len + 2 * x.intersp + 0.5) * xc
  avail <- 0.97 * graphics::par("pin")[1]
  ncol <- if (sum(item) + xc <= avail) k else max(1L, min(k, floor(avail / max(item))))
  nrow <- ceiling(k / ncol)
  list(ncol = ncol, nrow = nrow, height = (nrow + 0.6) * yc)
}

#' Upper y limit with room for a top legend of a given layout
#'
#' @keywords internal
#' @noRd
did_legroom <- function(ylim, lay) {
  pin <- graphics::par("pin")[2]
  h <- min(lay$height, 0.6 * pin)
  ylim[2] + diff(ylim) * h / max(0.92 * pin - h, 0.2 * pin)
}

#' Draw a top legend with a layout from did_leglayout()
#'
#' @keywords internal
#' @noRd
did_toplegend <- function(lab, lay, cex, ...) {
  horiz <- lay$ncol == length(lab) && lay$nrow == 1L
  tw <- if (horiz) did_textwidth(lab, cex) else {
    if (getRversion() < "4.1.0") NULL else max(graphics::strwidth(lab, cex = cex))
  }
  if (!horiz && !is.null(tw)) tw <- tw + graphics::strwidth("m", cex = cex)
  graphics::legend("top", legend = lab, horiz = horiz, ncol = if (horiz) 1 else lay$ncol,
                   bty = "o", bg = "white", box.col = NA, cex = cex, text.width = tw, ...)
  graphics::box()
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

#' Resolve the style arguments of the event-study plots
#'
#' `onset = NULL` puts the treatment-onset line at `0` for
#' `style = "classic"` and at `-0.5` otherwise; `onset_label = NA` labels it
#' "Intervention" for `"classic"` only; `caps = NULL` draws whisker caps for
#' `"classic"` only.
#'
#' @keywords internal
#' @noRd
did_style <- function(style, onset = NULL, onset_label = NA, caps = NULL) {
  style <- match.arg(style, c("points", "line", "classic"))
  classic <- style == "classic"
  if (is.null(onset)) onset <- if (classic) 0 else -0.5
  if (!is.numeric(onset) || length(onset) != 1L || !is.finite(onset)) {
    stop("'onset' must be a number", call. = FALSE)
  }
  if (length(onset_label) == 1L && is.na(onset_label)) {
    onset_label <- if (classic) "Intervention" else NULL
  }
  if (!is.null(onset_label) && (!is.character(onset_label) || length(onset_label) != 1L)) {
    stop("'onset_label' must be a character string or NULL", call. = FALSE)
  }
  if (identical(onset_label, "")) onset_label <- NULL
  if (is.null(caps)) caps <- classic
  if (!is.logical(caps) || length(caps) != 1L || is.na(caps)) {
    stop("'caps' must be TRUE or FALSE", call. = FALSE)
  }
  list(style = style, onset = onset, onset_label = onset_label, caps = caps)
}

#' Mark the reference (normalized) period of an event-study table
#'
#' A row is the reference when its level is in `ref` and its estimate is
#' exactly zero: an omitted or normalized period. Without such a row (for
#' example `attgt()` with `base = "varying"`), nothing is marked.
#'
#' @keywords internal
#' @noRd
did_refmark <- function(tab, ref) {
  if (!is.null(ref) && (!is.numeric(ref) || anyNA(ref))) {
    stop("'ref' must be numeric event times or NULL", call. = FALSE)
  }
  tab$ref <- tab$level %in% ref & is.finite(tab$estimate) & tab$estimate == 0
  tab
}

#' Shaded polygons through consecutive finite points
#'
#' @keywords internal
#' @noRd
did_ribbon <- function(x, lo, hi, col) {
  ok <- is.finite(x) & is.finite(lo) & is.finite(hi)
  run <- cumsum(!ok)
  for (r in unique(run[ok])) {
    s <- ok & run == r
    did_polygon(x[s], lo[s], hi[s], col)
  }
}

#' Vertical intervals, with caps when asked
#'
#' @keywords internal
#' @noRd
did_whiskers <- function(x, lo, hi, caps = FALSE, col = "black", lwd = 1.5) {
  if (!length(x)) return(invisible(NULL))
  graphics::segments(x, lo, x, hi, col = col, lwd = lwd)
  if (caps) {
    w <- 0.045 * graphics::par("cex")
    w <- graphics::grconvertX(w, "inches", "user") - graphics::grconvertX(0, "inches", "user")
    graphics::segments(x - w, lo, x + w, lo, col = col, lwd = lwd)
    graphics::segments(x - w, hi, x + w, hi, col = col, lwd = lwd)
  }
  invisible(NULL)
}

#' Vertical treatment-onset line with an optional label above the plot
#'
#' @keywords internal
#' @noRd
did_onset <- function(st) {
  graphics::abline(v = st$onset, lty = 2, col = "grey40")
  if (!is.null(st$onset_label)) {
    graphics::mtext(st$onset_label, side = 3, at = st$onset, line = 0.15, cex = 0.8,
                    col = "grey30")
  }
}

#' Shared drawer for event-study and aggregated DiD plots
#'
#' One look for [plot.panglm_attgt()], [plot.panglm_etwfe()] and
#' [plot.panglm_eventstudy()], in three styles. `"points"`: hollow points
#' before treatment, filled points after, an open diamond at the omitted
#' reference period, black pointwise intervals, uniform bands as light
#' grey bars behind them, a dotted zero line and a dashed line at treatment
#' onset (`-0.5`). `"line"`: the estimates joined by a line, the pointwise
#' intervals as a shaded polygon and uniform bands as a lighter one.
#' `"classic"`: filled points with capped whiskers, an open circle at the
#' reference period, a solid zero line and a dashed line at the onset (`0`)
#' labelled "Intervention".
#'
#' @param tab a table from did_plottab()
#' @param event logical: the levels are event times
#' @param xaxis logical: draw the x axis at the levels
#' @param style,onset,onset_label,caps see did_style()
#' @return `tab`, invisibly, with attribute `band` (were uniform bands drawn?)
#' @keywords internal
#' @noRd
did_drawes <- function(tab, event = TRUE, xlab = NULL, ylab = NULL, main = NULL,
                       ylim = NULL, legend = FALSE, level = 0.95, xaxis = TRUE,
                       style = "points", onset = NULL, onset_label = NA, caps = NULL, ...) {
  if (!is.logical(legend) || length(legend) != 1L || is.na(legend)) {
    stop("'legend' must be TRUE or FALSE", call. = FALSE)
  }
  st <- did_style(style, onset, onset_label, caps)
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
  if (legend) {
    pct <- paste0(format(100 * level), "% ")
    leg <- switch(st$style,
      points = list(
        lab = c(if (event && any(pre)) "Pre-treatment",
                if (event) "Post-treatment" else "Estimate",
                if (any(ref)) "Reference", paste0(pct, "pointwise CI"),
                if (band) paste0(pct, "uniform band")),
        pch = c(if (event && any(pre)) 21, 19, if (any(ref)) 23, NA, if (band) 15),
        lty = c(if (event && any(pre)) 0, 0, if (any(ref)) 0, 1, if (band) 0)),
      classic = list(
        lab = c(paste0("Estimate (", pct, "CI)"), if (any(ref)) "Reference period",
                if (band) paste0(pct, "uniform band")),
        pch = c(19, if (any(ref)) 1, if (band) 15),
        lty = c(1, if (any(ref)) 0, if (band) 0)),
      line = list(
        lab = c("Estimate", paste0(pct, "pointwise CI"), if (band) paste0(pct, "uniform band")),
        pch = c(NA, 15, if (band) 15),
        lty = c(1, 0, if (band) 0)))
    leg$col <- switch(st$style,
      line = c("black", grDevices::adjustcolor("grey30", 0.35),
               if (band) grDevices::adjustcolor("grey30", 0.15)),
      c(rep("black", length(leg$lab) - band), if (band) "grey85"))
    leg$pt.cex <- ifelse(leg$pch %in% 15, 2, 1.1)
    lay <- did_leglayout(leg$lab, 0.75, x.intersp = 0.6)
  }
  if (is.null(ylim)) {
    ylim <- range(lo[ok], hi[ok], est[ok], 0, if (band) c(blo[ok], bhi[ok]))
    if (legend) ylim[2] <- did_legroom(ylim, lay)
  }
  dx <- if (length(unique(x)) > 1L) min(diff(sort(unique(x)))) else 1
  dots <- list(...)
  if (is.null(dots$xlim)) dots$xlim <- range(x[ok]) + c(-0.5, 0.5) * dx
  do.call(graphics::plot, c(list(x[ok], est[ok], type = "n", xlab = xlab, ylab = ylab,
                                 main = main, ylim = ylim, xaxt = "n"), dots))
  if (xaxis) graphics::axis(1, at = sort(unique(x)))
  if (st$style == "classic") {
    graphics::abline(h = 0, col = "grey30")
  } else {
    graphics::abline(h = 0, col = "grey40", lty = 3)
  }
  if (event) did_onset(st)
  ci <- ok & !ref
  if (st$style == "line") {
    if (band) did_ribbon(x[ok], blo[ok], bhi[ok], grDevices::adjustcolor("grey30", 0.15))
    did_ribbon(x[ok], lo[ok], hi[ok], grDevices::adjustcolor("grey30", 0.35))
    graphics::lines(x[ok], est[ok], lwd = 1.8)
    graphics::points(x[ci], est[ci], pch = ifelse(pre[ci], 21, 19), bg = "white", cex = 0.8)
    if (any(ref)) graphics::points(x[ref], est[ref], pch = 23, bg = "white", cex = 0.9)
  } else {
    if (band) {
      b <- ok & !ref
      graphics::rect(x[b] - 0.2 * dx, blo[b], x[b] + 0.2 * dx, bhi[b], col = "grey85",
                     border = NA)
    }
    if (st$style == "classic") {
      did_whiskers(x[ci], lo[ci], hi[ci], caps = st$caps, lwd = 1.2)
      graphics::points(x[ci], est[ci], pch = 19, cex = 1.1)
      if (any(ref)) graphics::points(x[ref], est[ref], pch = 21, bg = "white", cex = 1.2)
    } else {
      did_whiskers(x[ci], lo[ci], hi[ci], caps = st$caps)
      graphics::points(x[ci], est[ci], pch = ifelse(pre[ci], 21, 19), bg = "white", cex = 1.1)
      if (any(ref)) graphics::points(x[ref], est[ref], pch = 23, bg = "white", cex = 1.1)
    }
  }
  if (legend) {
    did_toplegend(leg$lab, lay, 0.75, pch = leg$pch, lty = leg$lty, lwd = 1.5, col = leg$col,
                  pt.bg = "white", pt.cex = leg$pt.cex, seg.len = 1.2, x.intersp = 0.6)
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
  op <- savepar()
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

#' Compare event studies on one plot
#'
#' Draws the event-time effects of several fits on one axis, offset
#' horizontally, with pointwise confidence intervals and a legend: for
#' example the two-way fixed-effects event study of [eventstudy()], the
#' Callaway and Sant'Anna estimator of [attgt()] and the extended two-way
#' fixed-effects estimator of [etwfe()]. Uniform bands are not drawn
#' here; use the `plot()` method of a single fit for them.
#'
#' The styles are those of [plot.panglm_attgt()]. With `style = "classic"`
#' (the default) the points are filled, the whiskers capped, each
#' reference period is a hollow symbol at zero and the dashed onset line is
#' at event time `0`, labelled "Intervention". With `style = "points"` the
#' points are hollow before treatment and filled after, and the onset line
#' is at `-0.5`. With `style = "line"` each fit is a coloured line with its
#' pointwise intervals as a semi-transparent band (the fits are not
#' offset).
#'
#' @param ... fits from [eventstudy()], [attgt()] or [etwfe()] (or the
#'   event-time summaries of the last two). Names of the arguments are used
#'   as labels.
#' @param labels optional labels, one per fit (default: the argument names,
#'   else the estimator)
#' @param min_e,max_e event window to plot
#' @param col,pch colours and (filled) plotting symbols, recycled over the
#'   fits; the hollow version of each symbol is used before treatment
#' @param offset horizontal distance between the fits, in event-time units
#' @param legend logical: add a legend
#' @param xlab,ylab,main,ylim axis labels, title and y-axis limits
#' @param ref event time of the reference period, marked when the estimate
#'   of a fit is exactly zero there; the default (`NULL`) uses the omitted
#'   period of each fit (`-1`, or the `ref` of an [eventstudy()] fit)
#' @inheritParams plot.panglm_attgt
#' @return a data frame of the plotted estimates (`model`, `event`,
#'   `estimate`, `lower`, `upper`), invisibly
#' @family policy evaluation
#' @examples
#' set.seed(2)
#' d <- expand.grid(id = 1:150, t = 1:8)
#' d$first <- rep(c(4, 6, NA), each = 50)[d$id]
#' d$y <- rnorm(150)[d$id] + 0.2 * d$t +
#'   ifelse(!is.na(d$first) & d$t >= d$first, 1, 0) + rnorm(nrow(d))
#' es <- eventstudy(y ~ 1, d, index = c("id", "t"), treat = "first",
#'                  leads = 3, lags = 3)
#' cs <- attgt(y ~ 1, d, index = c("id", "t"), treat = "first")
#' ew <- etwfe(y ~ 1, d, index = c("id", "t"), treat = "first", control = "never")
#' eventplot(TWFE = es, CS = cs, ETWFE = ew, min_e = -3, max_e = 3)
#' eventplot(TWFE = es, ETWFE = ew, min_e = -3, max_e = 3, style = "points")
#' eventplot(TWFE = es, CS = cs, min_e = -3, max_e = 3, style = "line")
#' @export
eventplot <- function(..., labels = NULL, min_e = -Inf, max_e = Inf,
                      col = c("black", "#D55E00", "#0072B2", "#009E73", "#CC79A7"),
                      pch = c(19, 17, 15, 18, 25), offset = 0.2, legend = TRUE,
                      xlab = "Event time", ylab = "Estimate", main = NULL, ylim = NULL,
                      style = c("classic", "points", "line"), ref = NULL, onset = NULL,
                      onset_label = NA, caps = NULL) {
  fits <- list(...)
  if (!length(fits)) stop("give at least one fit", call. = FALSE)
  did_check_window(min_e, max_e)
  st <- did_style(match.arg(style), onset, onset_label, caps)
  tabs <- lapply(fits, did_eventtab, min_e = min_e, max_e = max_e, ref = ref)
  k <- length(fits)
  if (is.null(labels)) {
    labels <- names(fits)
    if (is.null(labels)) labels <- rep("", k)
    auto <- vapply(fits, did_estimator, character(1))
    labels[labels == ""] <- auto[labels == ""]
    labels <- make.unique(labels, sep = " ")
  }
  if (length(labels) != k) stop("give one label per fit", call. = FALSE)
  col <- rep_len(col, k)
  pch <- rep_len(pch, k)
  hollow <- c(`19` = 21, `16` = 21, `17` = 24, `15` = 22, `18` = 23, `25` = 25)
  pch_pre <- ifelse(as.character(pch) %in% names(hollow), hollow[as.character(pch)], pch)
  if (st$style == "line") offset <- 0
  shift <- (seq_len(k) - (k + 1) / 2) * offset
  all <- do.call(rbind, lapply(tabs, `[`, c("level", "estimate", "lower", "upper")))
  ok <- is.finite(all$estimate)
  if (legend) lay <- did_leglayout(labels, 0.8)
  if (is.null(ylim)) {
    ylim <- range(all$lower[ok], all$upper[ok], all$estimate[ok], 0, na.rm = TRUE)
    if (legend) ylim[2] <- did_legroom(ylim, lay)
  }
  xr <- range(all$level) + c(-1, 1) * max(abs(shift), 0.3)
  graphics::plot(NA, xlim = xr, ylim = ylim, xlab = xlab, ylab = ylab, main = main,
                 xaxt = "n")
  graphics::axis(1, at = sort(unique(all$level)))
  if (st$style == "classic") {
    graphics::abline(h = 0, col = "grey30")
  } else {
    graphics::abline(h = 0, col = "grey40", lty = 3)
  }
  did_onset(st)
  for (j in seq_len(k)) {
    t <- tabs[[j]]
    f <- is.finite(t$estimate)
    x <- t$level + shift[j]
    ci <- f & !t$ref & is.finite(t$lower) & is.finite(t$upper)
    if (st$style == "line") {
      lo <- ifelse(t$ref, t$estimate, t$lower)
      hi <- ifelse(t$ref, t$estimate, t$upper)
      did_ribbon(x, lo, hi, grDevices::adjustcolor(col[j], 0.2))
      graphics::lines(x[f], t$estimate[f], col = col[j], lwd = 1.8)
      p <- ifelse(t$level < 0 | t$ref, pch_pre[j], pch[j])
      graphics::points(x[f], t$estimate[f], pch = p[f], col = col[j], bg = "white", cex = 0.8)
    } else {
      did_whiskers(x[ci], t$lower[ci], t$upper[ci], caps = st$caps, col = col[j],
                   lwd = if (st$style == "classic") 1.2 else 1.5)
      p <- if (st$style == "classic") ifelse(t$ref, pch_pre[j], pch[j]) else
        ifelse(t$level < 0 | t$ref, pch_pre[j], pch[j])
      graphics::points(x[f], t$estimate[f], pch = p[f], col = col[j], bg = "white", cex = 1.1)
    }
  }
  if (legend) {
    did_toplegend(labels, lay, 0.8, col = col, pch = pch, lty = 1,
                  lwd = if (st$style == "line") 1.8 else 1.5, seg.len = 1.2)
  }
  out <- do.call(rbind, lapply(seq_len(k), function(j) {
    data.frame(model = labels[j], event = tabs[[j]]$level, estimate = tabs[[j]]$estimate,
               lower = tabs[[j]]$lower, upper = tabs[[j]]$upper)
  }))
  rownames(out) <- NULL
  invisible(out)
}

#' Event-time table of a DiD fit (for eventplot)
#'
#' @keywords internal
#' @noRd
did_eventtab <- function(x, min_e = -Inf, max_e = Inf, ref = NULL) {
  if (inherits(x, "panglm_eventstudy")) return(es_plottab(x, min_e, max_e, ref))
  if (inherits(x, "panglm_attgt")) {
    x <- summary.panglm_attgt(x, type = "dynamic", min_e = min_e, max_e = max_e)
  }
  if (inherits(x, "panglm_etwfe")) {
    x <- summary.panglm_etwfe(x, type = "event", byxvar = FALSE)
  }
  if (inherits(x, "summary.panglm_attgt") && identical(x$type, "dynamic")) {
    tab <- did_plottab(x$estimates, min_e = min_e, max_e = max_e)
    return(did_refmark(tab, if (is.null(ref)) -1 else ref))
  }
  if (inherits(x, "summary.panglm_etwfe") && identical(x$type, "event") &&
      is.null(x$xvar)) {
    return(etwfe_plottab(x, min_e, max_e, ref))
  }
  stop("eventplot() takes fits from eventstudy(), attgt() or etwfe(), or their ",
       "event-time summaries", call. = FALSE)
}

#' Default label of a DiD fit
#'
#' @keywords internal
#' @noRd
did_estimator <- function(x) {
  if (inherits(x, "panglm_eventstudy")) return("TWFE event study")
  if (inherits(x, c("panglm_attgt", "summary.panglm_attgt"))) return("Callaway-Sant'Anna")
  if (inherits(x, c("panglm_etwfe", "summary.panglm_etwfe"))) return("ETWFE")
  "model"
}

#' Event-study table of an eventstudy() fit
#'
#' @keywords internal
#' @noRd
es_plottab <- function(x, min_e = -Inf, max_e = Inf, ref = NULL) {
  tab <- x$estimates
  names(tab)[names(tab) == "event"] <- "level"
  tab <- did_plottab(tab, min_e = min_e, max_e = max_e)
  did_refmark(tab, if (is.null(ref)) x$ref else ref)
}

#' Plot table of an etwfe() summary; the reference period (event time -1
#' with control = "never") has a zero effect
#'
#' @keywords internal
#' @noRd
etwfe_plottab <- function(x, min_e = -Inf, max_e = Inf, ref = NULL) {
  tab <- did_plottab(x$estimates, min_e = min_e, max_e = max_e)
  if (x$type != "event") return(tab)
  did_refmark(tab, if (is.null(ref)) -1 else ref)
}

#' Diverging blue-white-red colours
#'
#' @keywords internal
#' @noRd
did_diverging <- function(n = 101L) {
  grDevices::colorRampPalette(c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"))(n)
}

#' Vertical colour bar in the right margin of the current plot
#'
#' @keywords internal
#' @noRd
did_colorbar <- function(pal, zlim, title = NULL, frac = c(0.25, 1)) {
  usr <- graphics::par("usr")
  xw <- graphics::grconvertX(c(0, 0.18), "inches", "user")
  xw <- diff(xw)
  x0 <- usr[2] + diff(graphics::grconvertX(c(0, 0.15), "inches", "user"))
  yb <- usr[3] + frac[1] * (usr[4] - usr[3])
  yt <- usr[3] + frac[2] * (usr[4] - usr[3])
  n <- length(pal)
  yy <- seq(yb, yt, length.out = n + 1L)
  graphics::rect(x0, yy[-(n + 1L)], x0 + xw, yy[-1L], col = pal, border = NA, xpd = NA)
  graphics::rect(x0, yb, x0 + xw, yt, border = "grey40", xpd = NA)
  at <- pretty(zlim, n = 5)
  at <- at[at >= zlim[1] & at <= zlim[2]]
  ya <- yb + (at - zlim[1]) / diff(zlim) * (yt - yb)
  graphics::segments(x0 + xw, ya, x0 + 1.3 * xw, ya, xpd = NA)
  graphics::text(x0 + 1.5 * xw, ya, format(at), adj = 0, cex = 0.7, xpd = NA)
  if (!is.null(title)) {
    graphics::text(x0, yt + 0.6 * graphics::strheight("M", cex = 0.75), title, adj = c(0, 0),
                   cex = 0.75, xpd = NA)
  }
  invisible(c(x0 = x0, xw = xw))
}

#' Cohort by period heat map of ATT(g, t)
#'
#' @param est data frame with `group`, `time`, `estimate`, `lower`, `upper`
#' @param sizes named cohort sizes (names are the cohorts), or NULL
#' @return `est` with `pre`, invisibly
#' @keywords internal
#' @noRd
did_heatmap <- function(est, sizes = NULL, cells = FALSE, digits = 2, xlab = "Period",
                        ylab = "Treatment cohort", main = NULL, zlim = NULL, legend = TRUE,
                        title = "ATT", ...) {
  for (a in c("cells", "legend")) {
    v <- get(a)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop("'", a, "' must be TRUE or FALSE", call. = FALSE)
    }
  }
  est <- est[is.finite(est$estimate), , drop = FALSE]
  if (!nrow(est)) stop("no ATT(g, t) to plot", call. = FALSE)
  groups <- sort(unique(est$group))
  periods <- sort(unique(est$time))
  ng <- length(groups)
  np <- length(periods)
  est$pre <- est$time < est$group
  if (is.null(zlim)) zlim <- c(-1, 1) * max(abs(est$estimate))
  if (zlim[1] == zlim[2]) zlim <- zlim + c(-1, 1)
  pal <- did_diverging(101L)
  z <- matrix(NA_real_, np, ng)
  ci <- match(est$time, periods)
  ri <- ng + 1L - match(est$group, groups)
  z[cbind(ci, ri)] <- pmin(pmax(est$estimate, zlim[1]), zlim[2])
  rlab <- format(groups)
  if (!is.null(sizes)) rlab <- paste0(rlab, " (n = ", sizes[as.character(groups)], ")")
  op <- savepar()
  on.exit(graphics::par(op))
  csi <- graphics::par("csi")
  lab_w <- max(graphics::strwidth(rlab, units = "inches", cex = 0.8))
  key_w <- if (legend) {
    max(graphics::strwidth(c(format(pretty(zlim)), "Pre-treatment", title), units = "inches",
                           cex = 0.75)) + 0.5
  } else 0
  graphics::par(mar = c(4, lab_w / csi + 2.6, if (is.null(main)) 1.5 else 3,
                        key_w / csi + 1))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.5, np + 0.5), ylim = c(0.5, ng + 0.5), xaxs = "i",
                        yaxs = "i")
  graphics::rect(0.5, 0.5, np + 0.5, ng + 0.5, col = "grey88", border = NA)
  graphics::image(seq_len(np), seq_len(ng), z, col = pal,
                  breaks = seq(zlim[1], zlim[2], length.out = length(pal) + 1L), add = TRUE,
                  ...)
  if (any(est$pre)) {
    p <- est$pre
    graphics::rect(ci[p] - 0.5, ri[p] - 0.5, ci[p] + 0.5, ri[p] + 0.5, density = 10,
                   angle = 45, col = "grey45", border = NA, lwd = 0.6)
  }
  graphics::abline(h = seq_len(ng - 1L) + 0.5, v = seq_len(np - 1L) + 0.5, col = "white",
                   lwd = 1)
  # adoption: a black edge before the first treated period of each cohort
  for (j in seq_along(groups)) {
    x0 <- sum(periods < groups[j]) + 0.5
    r <- ng + 1L - j
    if (x0 > 0.5 && x0 < np + 0.5) {
      graphics::segments(x0, r - 0.5, x0, r + 0.5, lwd = 2.5)
    }
  }
  if (cells) {
    sig <- is.finite(est$lower) & is.finite(est$upper) & (est$lower > 0 | est$upper < 0)
    lab <- paste0(formatC(est$estimate, format = "f", digits = digits), ifelse(sig, "*", ""))
    dark <- abs(est$estimate) > 0.6 * max(abs(zlim))
    cex <- min(0.75, 0.9 * graphics::par("pin")[1] / np /
                 max(graphics::strwidth(lab, units = "inches", cex = 1)))
    graphics::text(ci, ri, lab, cex = cex, font = ifelse(sig, 2, 1),
                   col = ifelse(dark, "white", "black"))
  }
  graphics::box()
  at <- seq_len(np)
  if (np > 20L) at <- at[periods %in% pretty(periods)]
  graphics::axis(1, at = at, labels = periods[at], cex.axis = 0.8)
  graphics::axis(2, at = ng + 1L - seq_len(ng), labels = rlab, las = 1, cex.axis = 0.8,
                 tick = FALSE, line = -0.5)
  graphics::title(xlab = xlab, line = 2.4)
  graphics::mtext(ylab, side = 2, line = lab_w / csi + 1.2)
  if (!is.null(main)) graphics::title(main = main)
  if (legend) {
    bar <- did_colorbar(pal, zlim, title = title, frac = c(if (any(est$pre)) 0.3 else 0, 0.9))
    if (any(est$pre)) {
      usr <- graphics::par("usr")
      h <- 0.08 * (usr[4] - usr[3])
      yb <- usr[3] + 0.05 * (usr[4] - usr[3])
      graphics::rect(bar["x0"], yb, bar["x0"] + bar["xw"], yb + h, col = "white",
                     border = "grey40", xpd = NA)
      graphics::rect(bar["x0"], yb, bar["x0"] + bar["xw"], yb + h, density = 10, angle = 45,
                     col = "grey45", border = NA, xpd = NA, lwd = 0.6)
      graphics::text(bar["x0"] + 1.5 * bar["xw"], yb + h / 2, "Pre-treatment", adj = 0,
                     cex = 0.7, xpd = NA)
    }
  }
  invisible(est)
}
