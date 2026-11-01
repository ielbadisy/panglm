#' Plot raw outcome trends of a difference-in-differences design
#'
#' Draws the mean outcome over time, the classical first look at a
#' difference-in-differences design before any model is fitted.
#'
#' * `type = "cohort"`: one line per adoption cohort and one for the
#'   never-treated units, against calendar time. Points are hollow before
#'   the cohort adopts and filled from adoption on.
#' * `type = "treated"`: the ever-treated and never-treated means against
#'   calendar time. A dashed vertical line marks a single adoption date
#'   (half a period before it); with staggered adoption each adoption date
#'   gets a dotted line (the ever-treated mean then mixes treated and not
#'   yet treated observations).
#' * `type = "event"`: the treated units aligned on event time (period
#'   minus first treated period), with a dashed line at treatment onset
#'   (`-0.5`). Never-treated units have no event time, so their mean cannot
#'   be aligned directly; the never-treated line is instead their mean in
#'   the same calendar periods as the treated observations: at event time
#'   `e`, each treated observation contributes the never-treated mean of
#'   its own period (a calendar-matched comparison, as in stacked event
#'   studies). The treated mean at each event time pools cohorts observed
#'   in different calendar periods, so its composition changes with `e`.
#'
#' Means are weighted by `weights` when given. The standard error of each
#' mean is cluster-robust by individual,
#' \deqn{\widehat{se}(\bar y_c) = \sqrt{\frac{G_c}{G_c - 1} \sum_{g}
#'   \Big(\sum_{i \in g} w_i (y_i - \bar y_c)\Big)^2 \Big/ W_c^2},}{se =
#'   sqrt(Gc / (Gc - 1) * sum_g (sum_i w_i (y_i - mean))^2) / Wc,}
#' with `G_c` the number of individuals contributing to the mean and `W_c`
#' the sum of their weights: the usual standard error of a mean (with one
#' observation per individual it is `s / sqrt(n)`), and the one of
#' `fixest::feols(y ~ 1, cluster = ~id)` fitted on the observations of the
#' mean. The calendar-matched never-treated line of `type = "event"` is a
#' fixed linear combination of period means; its standard error uses the
#' same formula with the influence of each never-treated observation.
#'
#' @param data a data frame in long format
#' @param index `c(individual, time)`: column names, the time numeric
#' @param treat name of the column with the first treated period, as in
#'   [attgt()] and [treatmentplot()] (`0`, `NA` or `Inf` for never-treated
#'   units)
#' @param outcome name of the outcome column (numeric or logical)
#' @param type `"cohort"`, `"treated"` or `"event"` (see Details)
#' @param weights optional sampling weights: the name of a column of `data`
#'   or a non-negative numeric vector with one value per row
#' @param band logical: add pointwise confidence bands (shaded) around each
#'   line
#' @param level confidence level of the bands
#' @param min_e,max_e event window for `type = "event"`
#' @param col line colours, one per line (cohorts first, the never treated
#'   last). The default gives ordered colours to the cohorts and black to
#'   the never treated.
#' @param legend logical: add a legend
#' @param xlab,ylab,main axis labels and title
#' @param ylim y-axis limits
#' @param ... passed to [graphics::plot()]
#' @return a data frame of the plotted means, invisibly: `group` (the line
#'   label), `time` (calendar period, or event time for `type = "event"`),
#'   `mean`, `se` and `n` (number of observations)
#' @family policy evaluation
#' @seealso [treatmentplot()] for the treatment timing, [panelplot()] for the
#'   unit by period values, [attgt()] and
#'   [etwfe()] for the estimators.
#' @examples
#' set.seed(2)
#' d <- expand.grid(id = 1:90, t = 1:8)
#' d$first <- rep(c(4, 6, 0), each = 30)[d$id]
#' d$y <- rnorm(90)[d$id] + 0.2 * d$t +
#'   ifelse(d$first > 0 & d$t >= d$first, 1, 0) + rnorm(nrow(d))
#' trendplot(d, index = c("id", "t"), treat = "first", outcome = "y")
#' trendplot(d, c("id", "t"), "first", "y", type = "treated", band = TRUE)
#' m <- trendplot(d, c("id", "t"), "first", "y", type = "event")
#' head(m)
#' @export
trendplot <- function(data, index, treat, outcome, type = c("cohort", "treated", "event"),
                      weights = NULL, band = FALSE, level = 0.95, min_e = -Inf, max_e = Inf,
                      col = NULL, legend = TRUE, xlab = NULL, ylab = NULL, main = NULL,
                      ylim = NULL, ...) {
  type <- match.arg(type)
  data <- as.data.frame(data)
  timing <- did_timing(data, index, treat)
  if (!is.character(outcome) || length(outcome) != 1L || !outcome %in% names(data)) {
    stop("'outcome' must be the name of one column of data", call. = FALSE)
  }
  y <- data[[outcome]]
  if (is.logical(y)) y <- as.numeric(y)
  if (!is.numeric(y)) stop("the outcome must be numeric or logical", call. = FALSE)
  for (a in c("band", "legend")) {
    v <- get(a)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop("'", a, "' must be TRUE or FALSE", call. = FALSE)
    }
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("'level' must be a number between 0 and 1", call. = FALSE)
  }
  did_check_window(min_e, max_e)
  w <- did_weights(weights, data)
  if (is.null(w)) w <- rep(1, nrow(data))
  id <- data[[index[1]]]
  ok <- !is.na(y) & !is.na(w) & !is.na(id) & is.finite(timing$time)
  id <- id[ok]
  y <- y[ok]
  w <- w[ok]
  time <- timing$time[ok]
  first <- timing$first[ok]
  if (!length(y)) stop("no observation with a non-missing outcome", call. = FALSE)
  never <- is.na(first)
  cohorts <- sort(unique(first[!never]))
  if (type == "cohort") {
    key <- ifelse(never, "Never treated", as.character(first))
    levs <- c(as.character(cohorts), if (any(never)) "Never treated")
    tab <- did_cellmeans(y, w, id, key, time)
  } else if (type == "treated") {
    key <- ifelse(never, "Never treated", "Treated")
    levs <- c(if (!all(never)) "Treated", if (any(never)) "Never treated")
    tab <- did_cellmeans(y, w, id, key, time)
  } else {
    if (all(never)) stop("type = 'event' needs treated units", call. = FALSE)
    e <- time - first
    use <- !never & e >= min_e & e <= max_e
    if (!any(use)) stop("no treated observation in the requested window", call. = FALSE)
    tab <- did_cellmeans(y[use], w[use], id[use], rep("Treated", sum(use)), e[use])
    levs <- "Treated"
    if (any(never)) {
      nv <- did_matchednever(y, w, id, time, never, e, use)
      if (!is.null(nv)) {
        tab <- rbind(tab, nv)
        levs <- c(levs, "Never treated")
      }
    }
  }
  tab <- tab[order(match(tab$group, levs), tab$time), , drop = FALSE]
  rownames(tab) <- NULL
  k <- length(levs)
  if (is.null(col)) {
    ncoh <- k - any(levs == "Never treated")
    ccol <- if (type == "cohort") {
      grDevices::hcl.colors(ncoh + 1L, "Plasma")[seq_len(ncoh)]
    } else {
      rep("#D55E00", ncoh)
    }
    col <- c(ccol, if (any(levs == "Never treated")) "black")
  }
  col <- rep_len(col, k)
  z <- stats::qnorm(1 - (1 - level) / 2)
  lo <- tab$mean - z * tab$se
  hi <- tab$mean + z * tab$se
  ukey <- if (type == "cohort") ifelse(never, "Never treated", as.character(first)) else
    ifelse(never, "Never treated", "Treated")
  sizes <- vapply(levs, function(g) length(unique(id[ukey == g])), numeric(1))
  lab <- paste0(if (type == "cohort") ifelse(levs == "Never treated", "", "Cohort "), levs,
                " (n = ", sizes, ")")
  if (is.null(xlab)) xlab <- if (type == "event") "Event time" else "Period"
  if (is.null(ylab)) ylab <- paste0("Mean ", outcome)
  if (legend) lay <- did_leglayout(lab, 0.8)
  if (is.null(ylim)) {
    ylim <- range(tab$mean, if (band) c(lo, hi), na.rm = TRUE, finite = TRUE)
    if (legend) ylim[2] <- did_legroom(ylim, lay)
  }
  xs <- sort(unique(tab$time))
  graphics::plot(tab$time, tab$mean, type = "n", xlab = xlab, ylab = ylab, main = main,
                 ylim = ylim, xaxt = "n", ...)
  at <- if (length(xs) <= 20L) xs else pretty(xs)
  graphics::axis(1, at = at)
  if (type == "event") graphics::abline(v = -0.5, lty = 2, col = "grey40")
  if (type == "treated" && length(cohorts)) {
    if (length(cohorts) == 1L) {
      graphics::abline(v = cohorts - 0.5, lty = 2, col = "grey40")
    } else {
      graphics::abline(v = cohorts - 0.5, lty = 3, col = "grey40")
    }
  }
  if (band) {
    for (j in seq_len(k)) {
      s <- tab$group == levs[j] & is.finite(lo) & is.finite(hi)
      did_polygon(tab$time[s], lo[s], hi[s], grDevices::adjustcolor(col[j], 0.18))
    }
  }
  for (j in seq_len(k)) {
    s <- tab$group == levs[j]
    tj <- tab$time[s]
    mj <- tab$mean[s]
    nv <- levs[j] == "Never treated"
    graphics::lines(tj, mj, col = col[j], lwd = 1.8)
    pre <- if (type == "cohort" && !nv) tj < as.numeric(levs[j]) else
      if (type == "event" && !nv) tj < 0 else rep(FALSE, length(tj))
    graphics::points(tj, mj, col = col[j], pch = if (nv) 15 else ifelse(pre, 21, 19),
                     bg = "white", cex = if (nv) 0.8 else 0.9)
  }
  if (legend) {
    did_toplegend(lab, lay, 0.8, col = col, lwd = 1.8, pch = ifelse(levs == "Never treated",
                                                                    15, 19),
                  seg.len = 1.6)
  }
  invisible(tab)
}

#' Weighted means by group and time with cluster-robust standard errors
#'
#' One row per `group` by `time` cell: `mean` (weighted), `se` (clustered
#' by `id`, small-sample factor `G / (G - 1)` with `G` the number of
#' individuals in the cell) and `n` (observations).
#'
#' @keywords internal
#' @noRd
did_cellmeans <- function(y, w, id, group, time) {
  cell <- paste(group, time, sep = "\r")
  u <- !duplicated(cell)
  key <- data.frame(group = group[u], time = time[u], cell = cell[u], stringsAsFactors = FALSE)
  f <- match(cell, key$cell)
  W <- rowsum(w, f)[, 1]
  m <- rowsum(w * y, f)[, 1] / W
  infl <- w * (y - m[f]) / W[f]
  out <- data.frame(group = key$group, time = key$time, mean = unname(m), se = NA_real_,
                    n = tabulate(f, nrow(key)), stringsAsFactors = FALSE)
  out$se <- did_clusterse(infl, f, id, nrow(key))
  out
}

#' Cluster-robust standard errors of linear statistics from influence values
#'
#' `infl` holds the influence of each observation on statistic `f`; the
#' variance of statistic `j` is `G_j / (G_j - 1)` times the sum over
#' clusters of the squared cluster sums (`NA` with fewer than two
#' clusters).
#'
#' @keywords internal
#' @noRd
did_clusterse <- function(infl, f, id, k) {
  cl <- paste(f, as.character(id), sep = "\r")
  u <- !duplicated(cl)
  ps <- as.numeric(rowsum(infl, cl, reorder = FALSE))
  fc <- f[u]
  G <- tabulate(fc, k)
  v <- numeric(k)
  s2 <- rowsum(ps^2, fc, reorder = TRUE)
  v[as.integer(rownames(s2))] <- s2[, 1]
  se <- sqrt(G / (G - 1) * v)
  se[G < 2] <- NA_real_
  se
}

#' Calendar-matched never-treated means by event time
#'
#' At event time `e`, the average over the treated observations at `e` of
#' the never-treated mean of their own period, with its cluster-robust
#' standard error (the period weights are held fixed).
#'
#' @keywords internal
#' @noRd
did_matchednever <- function(y, w, id, time, never, e, use) {
  nt <- never
  periods <- sort(unique(time[nt]))
  W <- as.numeric(tapply(w[nt], factor(time[nt], periods), sum))
  M <- as.numeric(tapply(w[nt] * y[nt], factor(time[nt], periods), sum)) / W
  es <- sort(unique(e[use]))
  rows <- lapply(es, function(ev) {
    s <- use & e == ev
    if (!all(time[s] %in% periods)) return(NULL)
    a <- tapply(w[s], factor(time[s], periods), sum)
    a[is.na(a)] <- 0
    a <- as.numeric(a) / sum(w[s])
    if (any(a > 0 & !is.finite(M))) return(NULL)
    p <- match(time[nt], periods)
    infl <- a[p] * w[nt] * (y[nt] - M[p]) / W[p]
    keep <- a[p] > 0
    se <- did_clusterse(infl[keep], rep(1L, sum(keep)), id[nt][keep], 1L)
    data.frame(group = "Never treated", time = ev, mean = sum(a[a > 0] * M[a > 0]), se = se,
               n = sum(keep), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Shaded band between two curves, broken at gaps
#'
#' @keywords internal
#' @noRd
did_polygon <- function(x, lo, hi, col, border = NA, ...) {
  if (length(x) < 2L) {
    if (length(x) == 1L) graphics::segments(x, lo, x, hi, col = col, lwd = 6, lend = 1)
    return(invisible(NULL))
  }
  o <- order(x)
  x <- x[o]
  lo <- lo[o]
  hi <- hi[o]
  graphics::polygon(c(x, rev(x)), c(lo, rev(hi)), col = col, border = border, ...)
  invisible(NULL)
}
