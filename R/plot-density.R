#' Model-implied marginal density of the outcome
#'
#' Returns `list(f, masses, y, weights, continuous, support, skewed, title)`:
#' `f(t)` is the average over observations of the model density of a
#' single outcome at the points `t` (for random-effects fits, integrated
#' over the random effects); `masses` is a data frame of point masses
#' (`at`, `predicted`, `observed`, `label`), empty when there are none;
#' `continuous` flags the observations of the continuous part.
#'
#' @keywords internal
#' @noRd
panglm_outcome_density <- function(object) {
  fam <- object$family$family
  model <- object$model
  if (fam %in% c("poisson", "negbin")) {
    stop("the density plot is for continuous outcomes; for counts use which = \"rootogram\"",
         call. = FALSE)
  }
  if (fam == "binomial" || inherits(object, "panordered")) {
    stop("the density plot is for continuous outcomes; for a ",
         if (fam == "binomial") "binary outcome use which = \"calibration\"" else
           "discrete outcome compare observed and predicted category shares",
         call. = FALSE)
  }
  if (!fam %in% c("gaussian", "gamma", "beta", "tobit")) {
    stop("the density plot is not available for family = '", fam, "'", call. = FALSE)
  }
  if (identical(model, "between")) {
    stop("the density plot is not available for between fits, which model ",
         "individual means rather than single outcomes", call. = FALSE)
  }
  y <- object$y
  eta <- object$linear.predictors
  n <- length(y)
  w <- object$weights %||% rep(1, n)
  random <- identical(model, "random")
  # standard deviation of the normal random effect on the link scale, per observation
  s_re <- if (!random) {
    rep(0, n)
  } else if (!is.null(object$random)) {
    sqrt(pmax(rowSums((object$Z %*% object$varcomp) * object$Z), 0))
  } else if (fam == "gaussian") {
    rep(sqrt(object$sigma_mu2), n)
  } else {
    rep(object$sigma, n)
  }
  avg <- function(M) as.numeric(crossprod(w, M)) / sum(w)
  masses <- data.frame(at = numeric(0), predicted = numeric(0), observed = numeric(0),
                       label = character(0))
  continuous <- rep(TRUE, n)
  support <- c(-Inf, Inf)
  skewed <- FALSE

  if (fam %in% c("gaussian", "tobit")) {
    base <- if (fam == "tobit") object$sigma_e else if (random) {
      sqrt(object$sigma_v2)
    } else {
      sqrt(object$dispersion)
    }
    sd <- sqrt(base^2 + s_re^2)
    f <- function(t) avg(stats::dnorm(outer(eta, t, function(e, v) (v - e)) / sd) / sd)
    if (fam == "tobit") {
      left <- object$family$left; right <- object$family$right
      support <- c(left, right)
      status <- tobit_status(y, left, right)
      continuous <- status == 0L
      if (is.finite(left)) {
        masses <- rbind(masses, data.frame(
          at = left, predicted = sum(w * stats::pnorm((left - eta) / sd)) / sum(w),
          observed = sum(w[status == -1L]) / sum(w), label = "left limit"))
      }
      if (is.finite(right)) {
        masses <- rbind(masses, data.frame(
          at = right, predicted = sum(w * stats::pnorm((eta - right) / sd)) / sum(w),
          observed = sum(w[status == 1L]) / sum(w), label = "right limit"))
      }
      title <- "Outcome density (Tobit)"
    } else {
      title <- "Outcome density (Gaussian)"
    }
  } else {
    gh <- marginal_quadrature(object$quadrature_nodes %||% 21L)
    z <- if (random) sqrt(2) * gh$nodes else 0
    qw <- if (random) gh$weights else 1
    if (fam == "gamma") {
      shape <- if (random) object$shape else 1 / object$dispersion
      cond <- function(t, e) stats::dgamma(t, shape = shape, rate = shape * exp(-e))
      support <- c(0, Inf)
      skewed <- TRUE
      title <- "Outcome density (gamma)"
    } else {
      phi <- object$precision
      mu_fun <- beta_link(object$family$link)$mu
      cond <- function(t, e) {
        mu <- mu_fun(e)
        stats::dbeta(t, mu * phi, (1 - mu) * phi)
      }
      support <- c(0, 1)
      title <- "Outcome density (beta)"
    }
    f <- function(t) {
      vapply(t, function(v) {
        dens <- 0
        for (q in seq_along(z)) dens <- dens + qw[q] * cond(v, eta + s_re * z[q])
        sum(w * dens) / sum(w)
      }, numeric(1))
    }
  }
  list(f = f, masses = masses, y = y, weights = w, continuous = continuous,
       support = support, skewed = skewed, title = title)
}

#' Model-implied marginal density of a two-part outcome
#'
#' Same value as `panglm_outcome_density()`, with the point mass at zero.
#'
#' @keywords internal
#' @noRd
twopart_outcome_density <- function(object) {
  data <- object$data
  y <- object$y
  n <- length(y)
  w <- object$weights %||% rep(1, n)
  a <- as.numeric(panglm_design(object$binary, data) %*% object$binary$coefficients)
  b <- as.numeric(panglm_design(object$positive, data) %*% object$positive$coefficients)
  Fb <- if (identical(object$link, "probit")) stats::pnorm else stats::plogis
  random <- identical(object$model, "random")
  gh <- marginal_quadrature(object$binary$quadrature_nodes %||% 21L)
  z <- sqrt(2) * gh$nodes
  qw <- gh$weights
  pos <- object$positive
  lognormal <- object$positive_type == "lognormal"

  if (isTRUE(object$correlated)) {
    j <- object$joint
    s1 <- j$sigma_binary
    shift <- j$rho * j$sigma_positive
    sdl <- sqrt(j$sigma_positive^2 * (1 - j$rho^2) + j$sigma_e^2)
    pq <- vapply(z, function(zz) Fb(a + s1 * zz), numeric(n))
    p <- as.numeric(pq %*% qw)
    f <- function(t) {
      vapply(t, function(v) {
        dens <- 0
        for (q in seq_along(z)) dens <- dens + qw[q] * pq[, q] * stats::dlnorm(v, b + shift * z[q], sdl)
        sum(w * dens) / sum(w)
      }, numeric(1))
    }
  } else {
    p <- if (random) {
      s1 <- object$binary$sigma
      as.numeric(vapply(z, function(zz) Fb(a + s1 * zz), numeric(n)) %*% qw)
    } else {
      Fb(a)
    }
    if (lognormal) {
      sdl <- if (random) sqrt(pos$sigma_v2 + pos$sigma_mu2) else sqrt(pos$dispersion)
      f <- function(t) {
        vapply(t, function(v) sum(w * p * stats::dlnorm(v, b, sdl)) / sum(w), numeric(1))
      }
    } else {
      shape <- if (random) pos$shape else 1 / pos$dispersion
      zz <- if (random) z * pos$sigma else 0
      ww <- if (random) qw else 1
      f <- function(t) {
        vapply(t, function(v) {
          dens <- 0
          for (q in seq_along(zz)) {
            dens <- dens + ww[q] * stats::dgamma(v, shape = shape, rate = shape * exp(-(b + zz[q])))
          }
          sum(w * p * dens) / sum(w)
        }, numeric(1))
      }
    }
  }
  masses <- data.frame(at = 0, predicted = 1 - sum(w * p) / sum(w),
                       observed = sum(w[y == 0]) / sum(w), label = "zero")
  list(f = f, masses = masses, y = y, weights = w, continuous = y > 0,
       support = c(0, Inf), skewed = TRUE,
       title = paste0("Outcome density (two-part, ",
                      if (lognormal) "lognormal" else "gamma", ")"))
}

#' Draw observed against model-implied outcome density
#'
#' The observed kernel density of the continuous part is scaled by its
#' share of the observations so that, with the point masses, it is on the
#' scale of the model-implied density. Positive, right-skewed outcomes are
#' drawn on a logarithmic axis by default (`log = "x"` forces it, `log = ""`
#' turns it off): the curves are then densities of `log(y)`, the implied one
#' `exp(s) f(exp(s))`, with the axis labelled in outcome units.
#'
#' @keywords internal
#' @noRd
.panglm_plot_density <- function(dens, bw = "nrd0", main = NULL, xlab = NULL, n = 512L, ...) {
  dots <- list(...)
  y <- dens$y
  w <- dens$weights
  cont <- dens$continuous & is.finite(y)
  yc <- y[cont]
  if (length(yc) < 2L) stop("too few continuous outcomes to estimate a density", call. = FALSE)
  share <- sum(w[cont]) / sum(w)
  masses <- dens$masses
  has_mass <- nrow(masses) > 0L

  logx <- if (is.null(dots$log)) {
    all(yc > 0) && dens$support[1] >= 0 &&
      stats::quantile(yc, 0.99, names = FALSE) > 20 * stats::median(yc)
  } else {
    grepl("x", dots$log, fixed = TRUE)
  }
  dots$log <- NULL
  if (logx && (any(yc <= 0) || dens$support[1] < 0)) {
    stop("a logarithmic axis needs a positive outcome", call. = FALSE)
  }
  tr <- if (logx) log else identity
  s <- tr(yc)
  lower <- if (logx) {
    if (dens$support[1] > 0) log(dens$support[1]) else -Inf
  } else {
    dens$support[1]
  }
  upper <- tr(dens$support[2])

  range_c <- if (!is.null(dots$xlim)) {
    tr(dots$xlim)
  } else {
    c(min(s), if (dens$skewed && !logx) stats::quantile(s, 0.99, names = FALSE) else max(s))
  }
  from <- max(range_c[1], lower)
  to <- min(range_c[2], upper)
  mass_at <- tr(pmax(masses$at, if (logx) 0 else -Inf))
  mass_at[!is.finite(mass_at)] <- -Inf
  xlim <- c(from, to)
  pad <- 0.14 * diff(xlim)
  if (has_mass && is.null(dots$xlim)) {
    if (any(mass_at <= from)) xlim[1] <- xlim[1] - pad
    if (any(mass_at >= to)) xlim[2] <- xlim[2] + pad
  }
  dots$xlim <- NULL

  wc <- w[cont] / sum(w[cont])
  if (is.character(bw)) bw <- stats::density(s, bw = bw)$bw
  kd <- if (all(wc == wc[1])) {
    stats::density(s, bw = bw, from = from, to = to, n = n)
  } else {
    stats::density(s, bw = bw, weights = wc, from = from, to = to, n = n)
  }
  observed <- kd$y * share
  grid <- kd$x
  # open interval at the limits so that densities with a pole stay finite
  inner <- grid
  eps <- 1e-6 * diff(range(grid))
  inner[1] <- inner[1] + eps
  inner[length(inner)] <- inner[length(inner)] - eps
  implied <- if (logx) exp(inner) * dens$f(exp(inner)) else dens$f(inner)
  implied[!is.finite(implied)] <- NA
  top <- max(c(observed, implied), na.rm = TRUE) * 1.12

  op <- savepar()
  on.exit(graphics::par(op), add = TRUE)
  if (has_mass) {
    mar <- graphics::par("mar")
    graphics::par(mar = c(mar[1:3], max(mar[4], 4.1)))
  }
  do.call(graphics::plot.default, c(list(
    NA, NA, xlim = xlim, ylim = c(0, top),
    main = main %||% dens$title,
    xlab = xlab %||% if (logx) "Outcome (log scale)" else "Outcome",
    ylab = if (logx) "Density of log(outcome)" else "Density",
    yaxs = "i", xaxt = if (logx) "n" else "s"), dots))
  if (logx) {
    k <- seq(floor(from / log(10)), ceiling(to / log(10)))
    at <- k * log(10)
    keep <- at >= from - 1e-8 & at <= to + 1e-8
    graphics::axis(1, at = at[keep],
                   labels = formatC(10^k[keep], format = "fg", big.mark = ","))
  }
  graphics::polygon(c(grid[1], grid, grid[length(grid)]), c(0, observed, 0),
                    col = "grey80", border = "grey45")
  graphics::lines(grid, implied, col = "#D55E00", lwd = 2)

  if (has_mass) {
    mass_top <- max(c(masses$predicted, masses$observed)) * 1.25
    ticks <- pretty(c(0, mass_top))
    ticks <- ticks[ticks <= mass_top]
    scale <- top / mass_top
    width <- 0.3 * pad
    for (i in seq_len(nrow(masses))) {
      left_side <- mass_at[i] <= from
      edge <- if (left_side) from else to
      centre <- edge + (if (left_side) -1 else 1) * pad / 2
      if (is.finite(mass_at[i]) && mass_at[i] >= xlim[1] && mass_at[i] <= xlim[2]) {
        graphics::abline(v = mass_at[i], lty = 3, col = "grey45")
      }
      xo <- centre - width - 0.02 * pad
      xp <- centre + 0.02 * pad
      graphics::rect(xo, 0, xo + width, masses$observed[i] * scale, col = "grey80",
                     border = "grey45")
      graphics::rect(xp, 0, xp + width, masses$predicted[i] * scale, col = "#D55E00",
                     border = "#A04000")
      graphics::text(c(xo, xp) + width / 2,
                     c(masses$observed[i], masses$predicted[i]) * scale,
                     sprintf("%.0f%%", 100 * c(masses$observed[i], masses$predicted[i])),
                     pos = 3, cex = 0.62, offset = 0.2, xpd = NA)
      if (logx || !is.finite(mass_at[i])) {
        graphics::axis(1, at = centre, labels = format(masses$at[i]), tick = FALSE)
      }
    }
    graphics::axis(4, at = ticks * scale, labels = format(ticks))
    graphics::mtext("Probability (point mass)", side = 4, line = 2.6)
  }
  graphics::legend(
    "topright", bty = "n", cex = 0.85,
    legend = c("Observed", "Model-implied"),
    fill = c("grey80", if (has_mass) "#D55E00" else NA),
    border = c("grey45", if (has_mass) "#A04000" else NA),
    lty = c(NA, 1), lwd = c(NA, 2), col = c(NA, "#D55E00")
  )
  invisible(list(x = grid, observed = observed, implied = implied, masses = masses,
                 log = logx))
}
