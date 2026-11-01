#' Plot a panel as a unit by period heat map
#'
#' Draws one row per individual and one column per period. With `value`
#' the cells are coloured by that variable: a sequential palette with a
#' colour bar for a numeric variable, one colour per level for a factor,
#' character or logical variable. Without `value` the plot shows which
#' periods each individual is observed in, so the gaps of an unbalanced
#' panel stand out. Periods in which an individual is not observed, and
#' missing values of `value`, are drawn in neutral grey.
#'
#' With `treat` the individuals are sorted by treatment cohort as in
#' [treatmentplot()] (earliest adopters at the top, never-treated units at
#' the bottom), the cohorts are separated by white lines and labelled on
#' the right with their sizes, and a black edge marks the first treated
#' period of each cohort. Without `treat` the individuals are in the order
#' of their identifiers. The unit labels are shown for at most 60
#' individuals by default; with many individuals the tiles are drawn as a
#' raster.
#'
#' @param data a data frame in long format
#' @param index `c(individual, time)`: column names; the time must be
#'   numeric when `treat` is given
#' @param value optional name of the column that colours the cells; `NULL`
#'   (the default) shows observed versus missing periods
#' @param treat optional name of the column with the first treated period,
#'   as in [attgt()] (`0`, `NA` or `Inf` for never-treated units)
#' @param labels logical: show the individual labels on the left axis. The
#'   default shows them for at most 60 individuals.
#' @param col colours: for a numeric `value`, the palette (default: a
#'   sequential blue palette); for a discrete `value`, one colour per level;
#'   without `value`, the colour of the observed cells
#' @param zlim limits of the colour scale for a numeric `value` (default:
#'   its range)
#' @param legend logical: add a legend (or colour bar) below the plot
#' @param xlab,ylab,main axis labels and title (default `ylab`: none)
#' @param ... passed to [graphics::image()]
#' @return a data frame with one row per individual in plotting order (`id`,
#'   `cohort` when `treat` is given, `n_periods` observed and `n_missing`
#'   missing values of `value`), invisibly
#' @family policy evaluation
#' @seealso [treatmentplot()] for the treatment timing alone, [trendplot()]
#'   for the mean outcome paths
#' @examples
#' set.seed(1)
#' d <- expand.grid(id = 1:30, year = 2010:2019)
#' d$first <- rep(c(2013, 2015, 0), each = 10)[d$id]
#' d$y <- rnorm(30)[d$id] + 0.1 * (d$year - 2010) + rnorm(nrow(d), sd = 0.3)
#' d$y[sample(nrow(d), 15)] <- NA
#' # observed and missing periods of an unbalanced panel
#' u <- d[-sample(nrow(d), 40), ]
#' panelplot(u, index = c("id", "year"))
#' # outcome values, sorted by treatment cohort
#' panelplot(d, index = c("id", "year"), value = "y", treat = "first")
#' # a discrete variable
#' d$high <- d$y > 0.5
#' panelplot(d, index = c("id", "year"), value = "high")
#' @export
panelplot <- function(data, index, value = NULL, treat = NULL, labels = NULL, col = NULL,
                      zlim = NULL, legend = TRUE, xlab = "Period", ylab = NULL, main = NULL,
                      ...) {
  data <- as.data.frame(data)
  if (!is.character(index) || length(index) != 2L) {
    stop("'index' must give the individual and time columns", call. = FALSE)
  }
  for (v in c(index, value, treat)) {
    if (!is.character(v) || length(v) != 1L || !v %in% names(data)) {
      stop("column '", v, "' not found in data", call. = FALSE)
    }
  }
  for (a in c("legend", if (!is.null(labels)) "labels")) {
    v <- get(a)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop("'", a, "' must be TRUE or FALSE", call. = FALSE)
    }
  }
  id <- data[[index[1]]]
  if (is.factor(id)) id <- droplevels(id)
  time <- data[[index[2]]]
  first <- if (is.null(treat)) NULL else did_timing(data, index, treat)$first
  ok <- !is.na(id) & !is.na(time)
  if (!any(ok)) stop("no observed individual and period", call. = FALSE)
  id <- id[ok]
  time <- time[ok]
  if (!is.null(first)) first <- first[ok]
  if (anyDuplicated(data.frame(id, time))) {
    stop("more than one row for some individual and period", call. = FALSE)
  }
  uid <- unique(id)
  uid <- uid[order(uid)]
  if (!is.null(first)) {
    ufirst <- first[match(uid, id)]
    ord <- order(is.na(ufirst), ufirst, seq_along(uid))
    uid <- uid[ord]
    ufirst <- ufirst[ord]
  }
  n <- length(uid)
  if (is.null(labels)) labels <- n <= 60L
  periods <- sort(unique(time))
  np <- length(periods)
  r <- match(id, uid)
  cell <- cbind(match(time, periods), n + 1L - r)

  # cell codes and colours
  numeric_value <- FALSE
  if (is.null(value)) {
    if (is.null(col)) col <- "#2C7FB8"
    if (length(col) != 1L) stop("'col' must give one colour without 'value'", call. = FALSE)
    code <- rep(1, length(id))
    pal <- col
    breaks <- c(0.5, 1.5)
    key <- c("Observed", "Not observed")
    keycol <- c(col, "grey88")
    miss <- rep(FALSE, length(id))
  } else {
    x <- data[[value]][ok]
    miss <- is.na(x)
    if (is.numeric(x)) {
      numeric_value <- TRUE
      if (!any(is.finite(x))) stop("'", value, "' has no finite value", call. = FALSE)
      if (is.null(zlim)) zlim <- range(x[is.finite(x)])
      if (!is.numeric(zlim) || length(zlim) != 2L || any(!is.finite(zlim)) ||
          zlim[1] > zlim[2]) {
        stop("'zlim' must be two increasing finite numbers", call. = FALSE)
      }
      if (zlim[1] == zlim[2]) zlim <- zlim + c(-0.5, 0.5)
      pal <- if (is.null(col)) grDevices::hcl.colors(100L, "Blues 3", rev = TRUE) else col
      if (length(pal) < 2L) stop("'col' must give at least two colours", call. = FALSE)
      breaks <- seq(zlim[1], zlim[2], length.out = length(pal) + 1L)
      code <- pmin(pmax(x, zlim[1]), zlim[2])
      code[!is.finite(x)] <- NA
      miss <- !is.finite(x)
    } else {
      f <- if (is.factor(x)) droplevels(x) else factor(x)
      lev <- levels(f)
      k <- length(lev)
      if (!k) stop("'", value, "' has no observed value", call. = FALSE)
      pal <- if (is.null(col)) {
        if (k == 1L) "#2C7FB8" else grDevices::hcl.colors(k + 1L, "YlGnBu", rev = TRUE)[-1L]
      } else col
      if (length(pal) != k) {
        stop("'col' must give one colour per level of '", value, "' (", k, ")", call. = FALSE)
      }
      code <- as.integer(f)
      breaks <- seq(0.5, k + 0.5)
      key <- c(lev, "Missing")
      keycol <- c(pal, "grey88")
    }
  }
  z <- matrix(NA_real_, np, n)
  z[cell] <- code

  out <- data.frame(id = uid)
  if (!is.null(first)) out$cohort <- ufirst
  out$n_periods <- tabulate(r, n)
  out$n_missing <- tabulate(r[miss], n)

  op <- savepar()
  on.exit(graphics::par(op))
  csi <- graphics::par("csi")
  lab_w <- if (labels) max(graphics::strwidth(format(uid), units = "inches", cex = 0.6)) else 0
  if (!is.null(first)) {
    cohorts <- unique(ufirst)
    sizes <- vapply(cohorts, function(g) sum(if (is.na(g)) is.na(ufirst) else ufirst %in% g),
                    numeric(1))
    coh_lab <- paste0(ifelse(is.na(cohorts), "Never", trimws(format(cohorts))), " (n = ",
                      sizes, ")")
    coh_w <- max(graphics::strwidth(coh_lab, units = "inches", cex = 0.8))
  } else {
    coh_w <- 0
  }
  graphics::par(mar = c(if (legend) 6 else 4, lab_w / csi + 1.5 + if (is.null(ylab)) 0 else 1,
                        if (is.null(main)) 1 else 3, coh_w / csi + 1))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.5, np + 0.5), ylim = c(0.5, n + 0.5), xaxs = "i",
                        yaxs = "i")
  graphics::rect(0.5, 0.5, np + 0.5, n + 0.5, col = "grey88", border = NA)
  graphics::image(seq_len(np), seq_len(n), z, col = pal, breaks = breaks, add = TRUE,
                  useRaster = n > 500L, ...)
  if (n <= 200L && np <= 100L) {
    graphics::abline(h = seq_len(n - 1L) + 0.5, v = seq_len(np - 1L) + 0.5, col = "white",
                     lwd = if (n > 60L) 0.3 else 0.8)
  }
  if (!is.null(main)) graphics::title(main = main)
  graphics::title(xlab = xlab, line = 2.2)
  at <- seq_len(np)
  if (np > 25L) {
    at <- unique(round(pretty(at, n = 12L)))
    at <- at[at >= 1 & at <= np]
  }
  graphics::axis(1, at = at, labels = periods[at], cex.axis = 0.8)
  if (labels) {
    graphics::axis(2, at = n + 1L - seq_len(n), labels = format(uid), las = 1, cex.axis = 0.6,
                   tick = FALSE, line = -0.6)
  }
  if (!is.null(ylab)) graphics::mtext(ylab, side = 2, line = lab_w / csi + 0.8)
  if (!is.null(first)) {
    ends <- cumsum(sizes)
    starts <- c(0, ends[-length(ends)])
    if (length(ends) > 1L) {
      graphics::abline(h = n + 0.5 - ends[-length(ends)], col = "white", lwd = 2.5)
    }
    for (j in seq_along(cohorts)) {
      if (is.na(cohorts[j])) next
      x0 <- sum(periods < cohorts[j]) + 0.5
      if (x0 > 0.5 && x0 < np + 0.5) {
        graphics::segments(x0, n + 0.5 - starts[j], x0, n + 0.5 - ends[j], lwd = 2.5)
      }
    }
    mid <- n + 1 - (ends - (sizes - 1) / 2)
    usr <- graphics::par("usr")
    h <- 1.3 * graphics::strheight("M", cex = 0.8)
    k <- length(mid)
    if (k > 1L) for (i in 2:k) mid[i] <- min(mid[i], mid[i - 1L] - h)
    mid[k] <- max(mid[k], usr[3] + h / 2)
    if (k > 1L) for (i in (k - 1L):1) mid[i] <- max(mid[i], mid[i + 1L] + h)
    graphics::text(usr[2] + graphics::strwidth("m", cex = 0.8), mid, coh_lab, adj = 0,
                   cex = 0.8, xpd = NA)
  }
  graphics::box()
  if (legend) {
    y_in <- graphics::grconvertY(0, "npc", "inches") - 4.3 * csi
    if (numeric_value) {
      panel_colorbar(pal, zlim, value, y_in, any(miss) || any(is.na(z)))
    } else {
      graphics::par(xpd = NA)
      graphics::legend(x = mean(graphics::par("usr")[1:2]),
                       y = graphics::grconvertY(y_in, "inches", "user"), legend = key,
                       fill = keycol, border = NA, horiz = TRUE, bty = "n", cex = 0.8,
                       xjust = 0.5, yjust = 0.5, text.width = did_textwidth(key, 0.8))
    }
  }
  invisible(out)
}

#' Horizontal colour bar below a plot
#'
#' @keywords internal
#' @noRd
panel_colorbar <- function(pal, zlim, title, y_in, missing = FALSE) {
  usr <- graphics::par("usr")
  w <- min(3, 0.5 * graphics::par("pin")[1])
  x_mid <- graphics::grconvertX(mean(usr[1:2]), "user", "inches")
  x0 <- graphics::grconvertX(x_mid - w / 2, "inches", "user")
  x1 <- graphics::grconvertX(x_mid + w / 2, "inches", "user")
  y0 <- graphics::grconvertY(y_in - 0.06, "inches", "user")
  y1 <- graphics::grconvertY(y_in + 0.06, "inches", "user")
  n <- length(pal)
  xx <- seq(x0, x1, length.out = n + 1L)
  graphics::rect(xx[-(n + 1L)], y0, xx[-1L], y1, col = pal, border = NA, xpd = NA)
  graphics::rect(x0, y0, x1, y1, border = "grey40", xpd = NA)
  at <- pretty(zlim, n = 4)
  at <- at[at >= zlim[1] & at <= zlim[2]]
  xa <- x0 + (at - zlim[1]) / diff(zlim) * (x1 - x0)
  dy <- y0 - graphics::grconvertY(y_in - 0.1, "inches", "user")
  graphics::segments(xa, y0, xa, y0 - dy, xpd = NA)
  graphics::text(xa, y0 - 1.5 * dy, format(at), adj = c(0.5, 1), cex = 0.7, xpd = NA)
  graphics::text(x0 - graphics::strwidth("m", cex = 0.8), (y0 + y1) / 2, title, adj = c(1, 0.5),
                 cex = 0.8, xpd = NA)
  if (missing) {
    sx <- x1 + graphics::strwidth("mm", cex = 0.8)
    sw <- diff(graphics::grconvertX(c(0, 0.15), "inches", "user"))
    graphics::rect(sx, y0, sx + sw, y1, col = "grey88", border = "grey40", xpd = NA)
    graphics::text(sx + sw + graphics::strwidth("m", cex = 0.8), (y0 + y1) / 2, "Missing",
                   adj = c(0, 0.5), cex = 0.8, xpd = NA)
  }
  invisible(NULL)
}
