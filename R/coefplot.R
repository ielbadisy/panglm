#' Compare coefficient estimates of several fitted models
#'
#' Dot-and-whisker plot of the coefficients of one or more fitted models,
#' terms on the vertical axis and Wald confidence intervals as horizontal
#' segments. Several models are drawn side by side for each term (dodged
#' vertically, with distinct plotting symbols and colours and a legend), so
#' that, for example, pooled, fixed-effects, random-effects and correlated
#' random effects fits of the same specification can be compared in one
#' figure. Every term present in at least one model is shown; a model that
#' does not contain a term has no mark on that row. The left margin is
#' widened to the longest term label so that no label is clipped.
#'
#' Any fitted object with `coef()` and `vcov()` methods can be used, e.g.
#' [panglm()], [paniv()], [pangmm()], [panfmm()] or [twopart()] fits, and
#' also [stats::glm()] fits. When `vcov()` fails or has no entry for a term,
#' only the estimate is drawn.
#'
#' @param ... fitted models, or a single unclassed list of fitted models.
#'   Argument names are used as legend labels.
#' @param parm names or positions (in the order the terms first appear
#'   across the models) of the terms to plot; by default every term except
#'   `"(Intercept)"`
#' @param level confidence level of the intervals
#' @param labels legend labels of the models; by default the argument names,
#'   or the deparsed arguments
#' @param main plot title
#' @param xlab horizontal-axis label
#' @param col colours of the models, recycled
#' @param pch plotting symbols of the models, recycled
#' @param ref position of the dashed reference line; `NULL` for none
#' @return the models, invisibly, as a named list
#' @family panglm methods
#' @examples
#' data(copd)
#' f <- exacerbations ~ crp + treatment + age
#' pool <- panglm(f, data = copd, index = c("id", "visit"),
#'                model = "pooling", family = "poisson")
#' re <- panglm(f, data = copd, index = c("id", "visit"),
#'              model = "random", family = "poisson")
#' cre <- panglm(f, data = copd, index = c("id", "visit"),
#'               model = "random", family = "poisson", cre = TRUE)
#' coefplot(Pooled = pool, RE = re, CRE = cre)
#' @export
coefplot <- function(..., parm = NULL, level = 0.95, labels = NULL,
                     main = "Coefficient estimates", xlab = "Estimate",
                     col = NULL, pch = NULL, ref = 0) {
  models <- list(...)
  if (length(models) == 1L && is.list(models[[1L]]) && is.null(oldClass(models[[1L]]))) {
    models <- models[[1L]]
    default_labels <- names(models)
  } else {
    default_labels <- names(models)
    exprs <- vapply(as.list(substitute(list(...)))[-1L],
                    function(e) paste(deparse(e, width.cutoff = 50L), collapse = ""),
                    character(1))
    if (is.null(default_labels)) default_labels <- exprs
    default_labels[default_labels == ""] <- exprs[default_labels == ""]
  }
  if (!length(models)) stop("supply at least one fitted model", call. = FALSE)
  if (is.null(default_labels)) default_labels <- paste("Model", seq_along(models))
  labels <- labels %||% default_labels
  if (length(labels) != length(models)) {
    stop("'labels' must have one entry per model", call. = FALSE)
  }
  names(models) <- labels
  tidy_check_level(level)
  parts <- lapply(models, function(m) {
    est <- tryCatch(stats::coef(m), error = function(e) NULL)
    if (!is.numeric(est) || is.null(names(est))) {
      stop("every model must have a named coef() vector", call. = FALSE)
    }
    V <- tryCatch(stats::vcov(m), error = function(e) NULL)
    list(est = est, se = .panglm_coef_se(est, V))
  })
  tab <- .panglm_coef_matrices(parts, parm, level)
  .panglm_dotwhisker(tab, main = main, xlab = xlab, col = col, pch = pch, ref = ref)
  invisible(models)
}

#' Standard errors matched to the coefficient names
#'
#' @keywords internal
#' @noRd
.panglm_coef_se <- function(est, V) {
  se <- rep(NA_real_, length(est))
  names(se) <- names(est)
  if (is.null(V) || !is.matrix(V)) return(se)
  d <- suppressWarnings(sqrt(diag(V)))
  if (!is.null(rownames(V)) && all(names(est) %in% rownames(V))) {
    se[] <- d[match(names(est), rownames(V))]
  } else if (length(d) == length(est)) {
    se[] <- d
  }
  se
}

#' Term-by-model estimate and interval matrices
#'
#' `parts` is a named list of `list(est, se)`. With `strict = FALSE`, names
#' in `parm` absent from all models are dropped instead of an error (used
#' when one `parm` is applied to several panels).
#'
#' @keywords internal
#' @noRd
.panglm_coef_matrices <- function(parts, parm, level, strict = TRUE) {
  terms <- unique(unlist(lapply(parts, function(p) names(p$est)), use.names = FALSE))
  if (!length(terms)) stop("the model has no coefficients to plot", call. = FALSE)
  if (is.null(parm)) {
    keep <- terms[!grepl("(Intercept)", terms, fixed = TRUE)]
    terms <- if (length(keep)) keep else terms
  } else if (is.numeric(parm)) {
    if (anyNA(parm) || any(parm < 1 | parm > length(terms))) {
      stop("'parm' must name or index terms of the fit", call. = FALSE)
    }
    terms <- terms[parm]
  } else {
    if (strict && !all(parm %in% terms)) {
      stop("'parm' must name or index terms of the fit", call. = FALSE)
    }
    terms <- parm[parm %in% terms]
    if (!length(terms)) stop("none of 'parm' is a term of the fit", call. = FALSE)
  }
  M <- length(parts)
  est <- se <- matrix(NA_real_, length(terms), M, dimnames = list(terms, names(parts)))
  for (m in seq_len(M)) {
    hit <- match(terms, names(parts[[m]]$est))
    est[, m] <- parts[[m]]$est[hit]
    se[, m] <- parts[[m]]$se[hit]
  }
  crit <- stats::qnorm(1 - (1 - level) / 2)
  list(est = est, lower = est - crit * se, upper = est + crit * se)
}

#' Shared dot-and-whisker drawing routine of all panglm coefficient plots
#'
#' Terms run from top to bottom; models are dodged within each term row.
#' The left margin is set from the longest label (in inches) and the top
#' margin makes room for the legend; both are restored on exit.
#'
#' @keywords internal
#' @noRd
.panglm_dotwhisker <- function(tab, main = NULL, xlab = "Estimate", col = NULL,
                               pch = NULL, ref = 0, legend = NULL, ...) {
  est <- tab$est
  lower <- tab$lower
  upper <- tab$upper
  terms <- rownames(est)
  M <- ncol(est)
  K <- length(terms)
  col <- rep_len(col %||% if (M == 1L) "black" else .panglm_palette, M)
  pch <- rep_len(pch %||% if (M == 1L) 19 else c(19, 17, 15, 18, 1, 2, 0, 5), M)
  legend <- legend %||% (M > 1L)

  xlim <- range(c(est, lower, upper, ref), finite = TRUE)
  if (!all(is.finite(xlim))) stop("no finite estimates to plot", call. = FALSE)
  if (diff(xlim) == 0) xlim <- xlim + c(-0.5, 0.5)

  op <- graphics::par(c("mai", "xpd"))
  on.exit(graphics::par(op), add = TRUE)
  mai <- op$mai
  cex_lab <- graphics::par("cex") * graphics::par("cex.axis")
  label_width <- max(graphics::strwidth(terms, units = "inches", cex = cex_lab))
  fin <- graphics::par("fin")[1]
  mai[2] <- min(max(mai[2], label_width + 0.25), 0.6 * fin)
  legend_rows <- 0L
  line_height <- graphics::par("csi")
  if (legend) {
    entry <- graphics::strwidth(colnames(est), units = "inches", cex = graphics::par("cex")) +
      0.45 * graphics::par("cex")
    avail <- fin - mai[2] - mai[4]
    per_row <- max(1L, min(M, floor(avail / max(entry))))
    legend_rows <- ceiling(M / per_row)
    mai[3] <- mai[3] + legend_rows * line_height * 1.1
  }
  graphics::par(mai = mai)

  positions <- rev(seq_len(K))
  offsets <- if (M == 1L) 0 else seq(0.28, -0.28, length.out = M)
  graphics::plot.default(NA, NA, xlim = xlim, ylim = c(0.5, K + 0.5), yaxt = "n",
                         ylab = "", xlab = xlab, main = "", ...)
  graphics::abline(h = positions, col = "grey92", lwd = 1)
  if (!is.null(ref)) graphics::abline(v = ref, lty = 2, col = "grey40")
  for (m in seq_len(M)) {
    y <- positions + offsets[m]
    ok <- is.finite(lower[, m]) & is.finite(upper[, m])
    graphics::segments(lower[ok, m], y[ok], upper[ok, m], y[ok], col = col[m], lwd = 1.5)
    graphics::points(est[, m], y, pch = pch[m], col = col[m], bg = col[m])
  }
  graphics::axis(2, at = positions, labels = terms, las = 1, tick = FALSE)
  if (!is.null(main)) {
    # the legend occupies legend_rows * 1.1 margin lines above the plot
    graphics::title(main = main, line = 1.5 + legend_rows * 1.1)
  }
  if (legend) {
    usr <- graphics::par("usr")
    graphics::par(xpd = NA)
    graphics::legend(mean(usr[1:2]), usr[4], legend = colnames(est), col = col, pch = pch,
                     lwd = 1.5, ncol = ceiling(M / legend_rows), bty = "n",
                     xjust = 0.5, yjust = 0, seg.len = 1.2, x.intersp = 0.6)
  }
  invisible(NULL)
}

# Okabe-Ito colours (colour-blind safe), black first
.panglm_palette <- c("#000000", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
                     "#56B4E9", "#F0E442")
