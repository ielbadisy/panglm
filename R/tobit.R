#' @details `tobitfamily()` is the censored normal (Tobit) model for
#'   outcomes piled up at a floor or a ceiling, such as health scores at
#'   the ceiling of a scale (motor function scores, EQ-5D utility at 1) or
#'   spending at zero: the latent outcome is `y* = x'beta + u_i + e_it`
#'   with `e_it ~ N(0, sigma_e^2)`, and `y = left` when `y* <= left`,
#'   `y = right` when `y* >= right`, `y = y*` otherwise. Use `left = -Inf`
#'   or `right = Inf` for one-sided censoring. `model = "pooling"` and
#'   `model = "random"` (Gaussian random intercept `u_i ~ N(0, sigma_u^2)`,
#'   adaptive Gauss-Hermite quadrature) are available. There is no
#'   fixed-effects Tobit (incidental parameters); use `model = "random"`
#'   with `cre = TRUE`. `predict(type = "link")` returns the latent mean
#'   `x'beta` and `predict(type = "response")` the mean of the observed,
#'   censored outcome (averaged over the random intercept), which is the
#'   scale of [partialeffects()]. `fit$sigma_e` and `fit$sigma` (random
#'   intercept SD) are stored on the fit.
#' @param left,right censoring limits of `tobitfamily()`
#' @references Tobin, J. (1958). Estimation of relationships for limited
#'   dependent variables. *Econometrica*, 26, 24-36. \doi{10.2307/1907382}
#' @rdname panglm-families
#' @export
tobitfamily <- function(left = 0, right = Inf) {
  if (!is.numeric(left) || !is.numeric(right) || length(left) != 1L ||
      length(right) != 1L || is.na(left) || is.na(right) || left >= right) {
    stop("'left' and 'right' must be single numbers with left < right", call. = FALSE)
  }
  if (!is.finite(left) && !is.finite(right)) {
    stop("at least one censoring limit must be finite; use family = 'gaussian' otherwise",
         call. = FALSE)
  }
  list(family = "tobit", link = "identity", family_id = NA_integer_, link_id = 0L,
       left = left, right = right)
}

#' Censoring codes: -1 left-censored, 1 right-censored, 0 observed
#'
#' @keywords internal
#' @noRd
tobit_status <- function(y, left, right) {
  if (any(y < left - 1e-10 * max(1, abs(left)) | y > right + 1e-10 * max(1, abs(right)))) {
    stop("the response must lie within the censoring limits [left, right]", call. = FALSE)
  }
  ifelse(y <= left, -1L, ifelse(y >= right, 1L, 0L))
}

#' Inverse Mills ratio `lambda(z) = dnorm(z) / pnorm(z)`, stable in the tails
#'
#' Below `z = -25` both `lambda(z)` and `z + lambda(z)` come from the
#' asymptotic series `z + lambda(z) = 1/x - 2/x^3 + 10/x^5 - ...`
#' (`x = -z`): the ratio of the logged densities loses digits there
#' (relative error 5e-3 at `z = -1e7`), and `z + lambda(z)` cancels.
#'
#' @keywords internal
#' @noRd
tobit_mills <- function(z) {
  out <- exp(stats::dnorm(z, log = TRUE) - stats::pnorm(z, log.p = TRUE))
  tail <- z < -25
  if (any(tail)) out[tail] <- tobit_mills_tail(-z[tail]) - z[tail]
  out
}

tobit_mills_tail <- function(x) {
  x2 <- x^2
  (1 - (2 - (10 - (74 - (706 - (8162 - 110410 / x2) / x2) / x2) / x2) / x2) / x2) / x
}

#' `z + lambda(z)`, stable in the tail
#'
#' `-lambda(z) (z + lambda(z))` is the second derivative of
#' `log pnorm(z)`, always in `(-1, 0)`. For `z -> -Inf`, `lambda(z) ~ -z`
#' and the sum cancels: at `z = -1e5` it came out negative, the curvature
#' positive and the quadrature scale `NaN`.
#'
#' @keywords internal
#' @noRd
tobit_mills_shift <- function(z, lam) {
  out <- z + lam
  tail <- z < -25
  if (any(tail)) out[tail] <- tobit_mills_tail(-z[tail])
  out
}

#' Observation log-likelihood and derivative with respect to the latent mean
#'
#' @keywords internal
#' @noRd
tobit_obs <- function(e, y, status, se, left, right) {
  ll <- numeric(length(y))
  d1 <- numeric(length(y))
  d2 <- numeric(length(y))
  o <- status == 0L
  if (any(o)) {
    r <- (y[o] - e[o]) / se
    ll[o] <- stats::dnorm(r, log = TRUE) - log(se)
    d1[o] <- r / se
    d2[o] <- -1 / se^2
  }
  l <- status == -1L
  if (any(l)) {
    z <- (left - e[l]) / se
    lam <- tobit_mills(z)
    ll[l] <- stats::pnorm(z, log.p = TRUE)
    d1[l] <- -lam / se
    d2[l] <- -lam * tobit_mills_shift(z, lam) / se^2
  }
  r_ <- status == 1L
  if (any(r_)) {
    z <- (e[r_] - right) / se
    lam <- tobit_mills(z)
    ll[r_] <- stats::pnorm(z, log.p = TRUE)
    d1[r_] <- lam / se
    d2[r_] <- -lam * tobit_mills_shift(z, lam) / se^2
  }
  list(ll = ll, d1 = d1, d2 = d2)
}

#' Pooled Tobit contributions, theta = (beta, log sigma_e)
#'
#' @keywords internal
#' @noRd
tobit_contrib_pooled <- function(X, y, status, left, right) {
  k <- ncol(X)
  function(theta, score) {
    if (!is.finite(theta[k + 1L]) || abs(theta[k + 1L]) > 25) {
      return(list(ll = rep(-Inf, length(y)), score = matrix(0, length(y), k + 1L)))
    }
    e <- as.numeric(X %*% theta[seq_len(k)])
    se <- exp(theta[k + 1L])
    ob <- tobit_obs(e, y, status, se, left, right)
    out <- list(ll = ob$ll)
    if (score) {
      # d ll / d log se, by the chain rule through the standardized residual
      dls <- numeric(length(y))
      o <- status == 0L
      dls[o] <- ((y[o] - e[o]) / se)^2 - 1
      cens <- status != 0L
      dls[cens] <- -ob$d1[cens] * se * ifelse(status[cens] == -1L,
                                              (left - e[cens]) / se * -1,
                                              (e[cens] - right) / se)
      out$score <- cbind(ob$d1 * X, dls)
    }
    out
  }
}

#' Random-intercept Tobit contributions by adaptive Gauss-Hermite quadrature
#'
#' theta = (beta, log sigma_e, log sigma_u). The integrand is log-concave in
#' u_i, so Newton steps find the mode; its curvature scales the rule. The
#' contributions and their exact scores are computed in C++
#' (`random_agq_cpp()`, kind 4, with log sigma_e in the place of the cut).
#'
#' @keywords internal
#' @noRd
tobit_contrib_random <- function(X, y, left, right, group_size, R) {
  k <- ncol(X); G <- length(group_size)
  gh <- gauss_hermite_quadrature(R)
  X <- unname(as.matrix(X)); storage.mode(X) <- "double"
  y <- as.numeric(y)
  group_size <- as.integer(group_size)
  group_start <- as.integer(c(0L, cumsum(group_size)[-G]))
  function(theta, score) {
    if (any(!is.finite(theta[k + 1:2])) || any(abs(theta[k + 1:2]) > 25)) {
      return(list(ll = rep(-Inf, G), score = matrix(0, G, k + 2L)))
    }
    res <- random_agq_cpp(theta[seq_len(k)], theta[k + 1L], exp(theta[k + 2L]), X, y,
                          group_start, group_size, gh$nodes, gh$weights, 4L, left, right)
    out <- list(ll = res$ll)
    if (score) out$score <- res$score
    out
  }
}

#' Mean of the censored outcome for a normal latent variable N(eta, s^2)
#'
#' @keywords internal
#' @noRd
tobit_censored_mean <- function(eta, s, left, right) {
  a <- (left - eta) / s
  b <- (right - eta) / s
  Pa <- if (is.finite(left)) stats::pnorm(a) else 0
  Pb <- if (is.finite(right)) stats::pnorm(b) else 1
  fa <- if (is.finite(left)) stats::dnorm(a) else 0
  fb <- if (is.finite(right)) stats::dnorm(b) else 0
  out <- eta * (Pb - Pa) + s * (fa - fb)
  if (is.finite(left)) out <- out + left * Pa
  if (is.finite(right)) out <- out + right * (1 - Pb)
  out
}

#' Total latent SD used for the (population-averaged) censored mean
#'
#' @keywords internal
#' @noRd
tobit_total_sd <- function(object) {
  if (object$model == "random") sqrt(object$sigma_e^2 + object$sigma^2) else object$sigma_e
}

#' Fit pooled or random-intercept Tobit
#'
#' The likelihood is maximized for the outcome divided by a power of two
#' close to its standard deviation (`y / c`, limits `left / c` and
#' `right / c`), with standardized regressors. Costs in dollars otherwise
#' put the slopes and the log standard deviations on scales that differ by
#' orders of magnitude, and the quasi-Newton search stalls. The estimates
#' are mapped back: `beta = c beta_s`, `log sigma = log c + log sigma_s`,
#' and the log-likelihood loses `log c` for every uncensored observation
#' (the density of `y` is that of `y / c` divided by `c`; the censoring
#' probabilities do not change). A power of two makes `y / c` exact.
#'
#' @keywords internal
#' @noRd
fit_tobit <- function(X, y, family, model, group_start, group_size, weights, R, maxit, tol) {
  k <- ncol(X)
  left <- family$left; right <- family$right
  status <- tobit_status(y, left, right)
  if (!any(status == 0L)) stop("the Tobit model needs uncensored observations", call. = FALSE)
  w_obs <- weights %||% rep(1, length(y))
  sc <- stats::sd(y)
  sc <- if (is.finite(sc) && sc > 0) 2^round(log2(sc)) else 1
  ys <- y / sc; ls <- left / sc; rs <- right / sc
  # starting values: least squares on the uncensored rows
  o <- status == 0L
  b0 <- as.numeric(stats::lm.wfit(X[o, , drop = FALSE], ys[o], w_obs[o])$coefficients)
  b0[!is.finite(b0)] <- 0
  s0 <- stats::sd(ys[o] - as.numeric(X[o, , drop = FALSE] %*% b0))
  if (!is.finite(s0) || s0 <= 0) s0 <- 1
  pooled <- ml_optimize(c(b0, log(s0)), tobit_contrib_pooled(X, ys, status, ls, rs),
                        w_obs, maxit, tol, scaling = ml_design_scaling(X, k + 1L),
                        warn = model == "pooling", what = "the Tobit fit")
  if (model == "pooling") {
    res <- pooled
    w_units <- w_obs
    w_uncens <- sum(w_obs[o])
  } else {
    w_units <- if (is.null(weights)) rep(1, length(group_size)) else individual_weights(weights, group_size)
    keep_ind <- w_units > 0
    keep_obs <- rep(keep_ind, group_size)
    contrib <- tobit_contrib_random(X[keep_obs, , drop = FALSE], ys[keep_obs], ls, rs,
                                    group_size[keep_ind], R)
    # split the pooled variance between the two components
    start <- c(pooled$theta[seq_len(k)], pooled$theta[k + 1L] + log(sqrt(0.6)),
               pooled$theta[k + 1L] + log(sqrt(0.4)))
    res <- ml_optimize(start, contrib, w_units[keep_ind], maxit, tol,
                       scaling = ml_design_scaling(X[keep_obs, , drop = FALSE], k + 2L),
                       what = "the random-effects Tobit fit")
    score <- matrix(0, length(w_units), length(res$theta))
    score[keep_ind, ] <- res$score
    res$score <- score
    w_uncens <- sum(rep(w_units, group_size)[o])
  }
  # back to the scale of y
  p <- length(res$theta)
  a <- c(rep(sc, k), rep(1, p - k))
  res$theta <- a * res$theta + c(rep(0, k), rep(log(sc), p - k))
  res$vcov <- res$vcov * outer(a, a)
  res$score <- sweep(res$score, 2, a, "/")
  res$loglik <- res$loglik - w_uncens * log(sc)
  coefs <- res$theta[seq_len(k)]
  names(coefs) <- colnames(X)
  vcov <- res$vcov[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(vcov) <- list(colnames(X), colnames(X))
  sigma_e <- exp(res$theta[k + 1L])
  sigma <- if (model == "random") exp(res$theta[k + 2L]) else NULL
  eta <- as.numeric(X %*% coefs)
  s_total <- if (model == "random") sqrt(sigma_e^2 + sigma^2) else sigma_e
  list(coefficients = coefs, vcov = vcov, bread = vcov,
       fitted.values = tobit_censored_mean(eta, s_total, left, right),
       linear.predictors = eta, loglik = res$loglik, sigma_e = sigma_e, sigma = sigma,
       dispersion = sigma_e^2, theta = res$theta, theta_vcov = res$vcov,
       theta_score = res$score, unit_weights = w_units,
       censoring = c(left = sum(status == -1L), uncensored = sum(status == 0L),
                     right = sum(status == 1L)),
       quadrature_nodes = if (model == "random") R else NULL,
       npar = length(res$theta), df.residual = sum(w_obs > 0) - length(res$theta),
       iterations = res$iterations, convergence = res$convergence)
}
