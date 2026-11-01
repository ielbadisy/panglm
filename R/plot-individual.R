#' Plot the distribution of estimated individual effects
#'
#' Fixed effects: histogram and kernel density of [fixef()]. Random
#' effects: histogram and kernel density of the posterior means [ranef()]
#' with the assumed mixing density overlaid; with random slopes one panel
#' per effect followed by a scatter of the intercept against each slope
#' with the 95% contour of the estimated normal distribution.
#'
#' @keywords internal
#' @noRd
.panglm_plot_effects <- function(x, bw = "nrd0", main = NULL, xlab = NULL, ...) {
  if (identical(x$model, "within")) {
    fe <- fixef(x)
    if (is.list(fe)) fe <- fe$individual
    .panglm_effect_hist(fe, bw, density = NULL,
                        main = main %||% "Estimated fixed effects",
                        xlab = xlab %||% "Individual effect",
                        sub = sprintf("%d individuals, SD %.3g", length(fe), stats::sd(fe)), ...)
    return(invisible(NULL))
  }
  if (!identical(x$model, "random")) {
    stop("the effects plot needs a fixed-effects (model = 'within') or random-effects ",
         "(model = 'random') fit", call. = FALSE)
  }
  re <- ranef(x)
  fam <- x$family$family
  if (ncol(re) == 1L) {
    prior <- panglm_effect_prior(x)
    b <- re[, 1]
    .panglm_effect_hist(
      b, bw, density = prior$density,
      main = main %||% "Predicted random effects",
      xlab = xlab %||% prior$label,
      sub = sprintf("%d individuals; SD of predictions %.3g, model SD %.3g (shrinkage)",
                    length(b), stats::sd(b), prior$sd), ...)
    return(invisible(NULL))
  }
  Sigma <- x$varcomp
  q <- ncol(re)
  labs <- colnames(re) %||% paste0("b", seq_len(q))
  npanel <- 2L * q - 1L
  op <- graphics::par(mfrow = c(ceiling(npanel / 3), min(npanel, 3L)))
  on.exit(graphics::par(op), add = TRUE)
  for (j in seq_len(q)) {
    sj <- sqrt(Sigma[j, j])
    .panglm_effect_hist(
      re[, j], bw, density = function(v) stats::dnorm(v, 0, sj),
      main = if (j == 1L) main %||% "Predicted random effects" else "",
      xlab = paste("Random effect:", labs[j]),
      sub = sprintf("SD of predictions %.3g, model SD %.3g", stats::sd(re[, j]), sj), ...)
  }
  ang <- seq(0, 2 * pi, length.out = 200)
  radius <- sqrt(stats::qchisq(0.95, 2))
  for (j in seq_len(q)[-1L]) {
    S <- Sigma[c(1L, j), c(1L, j)]
    e <- eigen(S, symmetric = TRUE)
    ell <- radius * cbind(cos(ang), sin(ang)) %*% t(e$vectors %*% diag(sqrt(pmax(e$values, 0))))
    graphics::plot.default(
      re[, 1], re[, j], xlim = range(re[, 1], ell[, 1]), ylim = range(re[, j], ell[, 2]),
      xlab = labs[1], ylab = labs[j], main = "",
      pch = 16, cex = 0.6, col = .panglm_scatter_col(nrow(re)), ...)
    graphics::lines(ell[, 1], ell[, 2], col = "#D55E00", lwd = 2)
    graphics::abline(h = 0, v = 0, lty = 3, col = "grey45")
    graphics::mtext(sprintf("correlation of predictions %.2f, model %.2f",
                            stats::cor(re[, 1], re[, j]), S[1, 2] / sqrt(S[1, 1] * S[2, 2])),
                    side = 3, line = 0.3, cex = 0.75, col = "grey30")
  }
  invisible(NULL)
}

#' Assumed mixing density of a random intercept, on the scale of ranef()
#'
#' @keywords internal
#' @noRd
panglm_effect_prior <- function(x) {
  fam <- x$family$family
  if (fam == "poisson") {
    d <- x$dispersion_param
    return(list(density = function(v) exp(d * log(d) - lgamma(d) + d * v - d * exp(v)),
                sd = sqrt(trigamma(d)),
                label = "Random intercept, log(nu) (nu ~ Gamma(d, d))"))
  }
  if (fam == "negbin") {
    a <- x$negbin_a; b <- x$negbin_b
    return(list(density = function(v) {
      exp(b * stats::plogis(v, log.p = TRUE) + a * stats::plogis(-v, log.p = TRUE) - lbeta(b, a))
    }, sd = sqrt(trigamma(a) + trigamma(b)),
    label = "Random intercept, log(delta)"))
  }
  s <- if (fam == "gaussian") sqrt(x$sigma_mu2) else x$sigma
  list(density = function(v) stats::dnorm(v, 0, s), sd = s,
       label = "Random intercept (linear-predictor scale)")
}

#' Histogram with kernel density and an optional reference density
#'
#' @keywords internal
#' @noRd
.panglm_effect_hist <- function(b, bw, density, main, xlab, sub = NULL, ...) {
  b <- b[is.finite(b)]
  if (length(b) < 2L) stop("too few finite effects to plot", call. = FALSE)
  # Freedman-Diaconis bins for many effects, Sturges otherwise
  fd <- 2 * stats::IQR(b) * length(b)^(-1 / 3)
  nb <- if (length(b) > 200L && fd > 0) {
    min(60L, ceiling(diff(range(b)) / fd))
  } else {
    ceiling(log2(length(b)) + 1)
  }
  h <- graphics::hist(b, breaks = pretty(range(b), max(nb, 5L)), plot = FALSE)
  kd <- stats::density(b, bw = bw)
  xr <- range(h$breaks, kd$x)
  ref <- NULL
  if (!is.null(density)) {
    grid <- seq(xr[1], xr[2], length.out = 400)
    ref <- list(x = grid, y = density(grid))
    ref$y[!is.finite(ref$y)] <- NA
  }
  top <- max(h$density, kd$y, ref$y, na.rm = TRUE) * 1.1
  graphics::plot.default(NA, NA, xlim = xr, ylim = c(0, top), main = main, xlab = xlab,
                         ylab = "Density", yaxs = "i", ...)
  graphics::rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1L], h$density,
                 col = "grey80", border = "grey45")
  graphics::lines(kd$x, kd$y, col = "#0072B2", lwd = 2)
  if (!is.null(ref)) graphics::lines(ref$x, ref$y, col = "#D55E00", lwd = 2, lty = 2)
  if (!is.null(sub)) graphics::mtext(sub, side = 3, line = 0.3, cex = 0.75, col = "grey30")
  graphics::legend("topright", bty = "n", cex = 0.8,
                   legend = c("Kernel density", if (!is.null(ref)) "Assumed distribution"),
                   col = c("#0072B2", if (!is.null(ref)) "#D55E00"),
                   lty = c(1, if (!is.null(ref)) 2), lwd = 2)
}
