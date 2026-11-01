#' Heatmap of the correlation matrix of the estimates
#'
#' Draws `cov2cor(V)` as a heatmap with a diverging palette centred at
#' zero (blue negative, white zero, red positive) on the fixed range
#' `[-1, 1]`, the term labels on both axes with margins sized to them, an
#' optional printed value per cell and a colour key in the right margin.
#' Terms with a zero or undefined variance are drawn grey.
#'
#' `parm` selects rows and columns of `V` by name or position; by default
#' every term whose name does not contain `"(Intercept)"`.
#'
#' @keywords internal
#' @noRd
.panglm_corrplot <- function(V, parm = NULL, cells = FALSE, main = NULL, ...) {
  V <- as.matrix(V)
  if (is.null(V) || !length(V) || nrow(V) != ncol(V)) {
    stop("the fit has no covariance matrix of the estimates", call. = FALSE)
  }
  terms <- rownames(V) %||% paste0("b", seq_len(nrow(V)))
  dimnames(V) <- list(terms, terms)
  if (is.null(parm)) {
    keep <- !grepl("(Intercept)", terms, fixed = TRUE)
    sel <- if (sum(keep) >= 2L) which(keep) else seq_along(terms)
  } else if (is.numeric(parm)) {
    if (anyNA(parm) || any(parm < 1 | parm > length(terms))) {
      stop("'parm' must name or index terms of the fit", call. = FALSE)
    }
    sel <- parm
  } else {
    if (!all(parm %in% terms)) stop("'parm' must name or index terms of the fit", call. = FALSE)
    sel <- match(parm, terms)
  }
  if (length(sel) < 2L) stop("the correlation plot needs at least two terms", call. = FALSE)
  V <- V[sel, sel, drop = FALSE]
  d <- diag(V)
  ok <- is.finite(d) & d > 0
  R <- matrix(NA_real_, nrow(V), ncol(V), dimnames = dimnames(V))
  R[ok, ok] <- stats::cov2cor(V[ok, ok, drop = FALSE])
  R[R > 1] <- 1
  R[R < -1] <- -1
  k <- nrow(R)
  labs <- rownames(R)

  pal <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  lab_cex <- if (k > 25) 0.6 else if (k > 12) 0.75 else 0.85
  lab_w <- max(graphics::strwidth(labs, units = "inches", cex = lab_cex))
  mai <- graphics::par("mai")
  graphics::par(mai = c(lab_w + 0.35, lab_w + 0.3, mai[3], 1.05), pty = "s")
  # row i of R at the top: image() draws z[x, y] with y upwards
  z <- t(R[rev(seq_len(k)), , drop = FALSE])
  graphics::image(seq_len(k), seq_len(k), z, zlim = c(-1, 1), col = pal,
                  axes = FALSE, xlab = "", ylab = "",
                  main = main %||% "Correlation of the estimates", ...)
  if (anyNA(z)) {
    na <- which(is.na(z), arr.ind = TRUE)
    graphics::rect(na[, 1] - 0.5, na[, 2] - 0.5, na[, 1] + 0.5, na[, 2] + 0.5,
                   col = "grey70", border = NA)
  }
  graphics::abline(h = seq_len(k + 1L) - 0.5, v = seq_len(k + 1L) - 0.5, col = "white", lwd = 0.8)
  graphics::box(col = "grey45")
  graphics::axis(1, at = seq_len(k), labels = labs, las = 2, tick = FALSE, cex.axis = lab_cex,
                 line = -0.6)
  graphics::axis(2, at = seq_len(k), labels = rev(labs), las = 1, tick = FALSE,
                 cex.axis = lab_cex, line = -0.6)
  if (isTRUE(cells)) {
    xy <- expand.grid(x = seq_len(k), y = seq_len(k))
    val <- z[cbind(xy$x, xy$y)]
    show <- is.finite(val)
    cell_cex <- max(0.35, min(0.9, 7 / k))
    txt <- sub("^-(0\\.00)$", "\\1", sprintf("%.2f", val[show]))
    graphics::text(xy$x[show], xy$y[show], txt, cex = cell_cex,
                   col = ifelse(abs(val[show]) > 0.6, "white", "grey15"))
  }
  # colour key in the right margin
  usr <- graphics::par("usr")
  span <- usr[2] - usr[1]
  x0 <- usr[2] + 0.05 * span
  x1 <- x0 + 0.045 * span
  yb <- seq(usr[3], usr[4], length.out = length(pal) + 1L)
  graphics::rect(x0, yb[-length(yb)], x1, yb[-1L], col = pal, border = NA, xpd = NA)
  graphics::rect(x0, usr[3], x1, usr[4], border = "grey45", xpd = NA)
  at <- c(-1, -0.5, 0, 0.5, 1)
  ya <- usr[3] + (at + 1) / 2 * (usr[4] - usr[3])
  graphics::segments(x1, ya, x1 + 0.012 * span, ya, xpd = NA, col = "grey30")
  graphics::text(x1 + 0.02 * span, ya, format(at), adj = 0, cex = 0.75, xpd = NA)
  invisible(R)
}
