#' @details `tobit_family()` is the censored normal (Tobit) model for
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
#' @param left,right censoring limits of `tobit_family()`
#' @references Tobin, J. (1958). Estimation of relationships for limited
#'   dependent variables. *Econometrica*, 26, 24-36. \doi{10.2307/1907382}
#' @rdname panglm-families
#' @export
tobit_family <- function(left = 0, right = Inf) {
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

# log of the inverse Mills ratio dnorm(z) / pnorm(z), stable in the tails
tobit_mills <- function(z) exp(stats::dnorm(z, log = TRUE) - stats::pnorm(z, log.p = TRUE))

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
    d2[l] <- -lam * (z + lam) / se^2
  }
  r_ <- status == 1L
  if (any(r_)) {
    z <- (e[r_] - right) / se
    lam <- tobit_mills(z)
    ll[r_] <- stats::pnorm(z, log.p = TRUE)
    d1[r_] <- lam / se
    d2[r_] <- -lam * (z + lam) / se^2
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
#' u_i, so Newton steps find the mode; its curvature scales the rule.
#'
#' @keywords internal
#' @noRd
tobit_contrib_random <- function(X, y, status, left, right, group, R) {
  k <- ncol(X); G <- max(group)
  gh <- gauss_hermite_quadrature(R)
  log_w <- log(gh$weights * sqrt(pi)) + gh$nodes^2
  function(theta, score) {
    if (any(!is.finite(theta[k + 1:2])) || any(abs(theta[k + 1:2]) > 25)) {
      return(list(ll = rep(-Inf, G)))
    }
    eta <- as.numeric(X %*% theta[seq_len(k)])
    se <- exp(theta[k + 1L]); su <- exp(theta[k + 2L])
    h <- function(u) {
      as.numeric(rowsum(tobit_obs(eta + u[group], y, status, se, left, right)$ll, group)) +
        stats::dnorm(u, 0, su, log = TRUE)
    }
    u <- rep(0, G)
    for (step in seq_len(100L)) {
      ob <- tobit_obs(eta + u[group], y, status, se, left, right)
      g <- as.numeric(rowsum(ob$d1, group)) - u / su^2
      H <- as.numeric(rowsum(ob$d2, group)) - 1 / su^2
      delta <- pmax(-2, pmin(2, g / H))
      if (anyNA(delta)) return(list(ll = rep(-Inf, G)))
      u <- u - delta
      if (max(abs(delta)) < 1e-10) break
    }
    ob <- tobit_obs(eta + u[group], y, status, se, left, right)
    tau <- 1 / sqrt(-(as.numeric(rowsum(ob$d2, group)) - 1 / su^2))
    terms <- vapply(seq_along(gh$nodes), function(q) {
      log_w[q] + h(u + sqrt(2) * tau * gh$nodes[q])
    }, numeric(G))
    if (!is.matrix(terms)) terms <- matrix(terms, nrow = G)
    list(ll = log(sqrt(2) * tau) + log_row_sum_exp(terms))
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
#' @keywords internal
#' @noRd
fit_tobit <- function(X, y, family, model, group_start, group_size, weights, R, maxit, tol) {
  k <- ncol(X)
  left <- family$left; right <- family$right
  status <- tobit_status(y, left, right)
  if (!any(status == 0L)) stop("the Tobit model needs uncensored observations", call. = FALSE)
  w_obs <- weights %||% rep(1, length(y))
  # starting values: least squares on the uncensored rows
  o <- status == 0L
  b0 <- as.numeric(stats::lm.wfit(X[o, , drop = FALSE], y[o], w_obs[o])$coefficients)
  b0[!is.finite(b0)] <- 0
  s0 <- stats::sd(y[o] - as.numeric(X[o, , drop = FALSE] %*% b0))
  if (!is.finite(s0) || s0 <= 0) s0 <- 1
  pooled <- ml_optimize(c(b0, log(s0)), tobit_contrib_pooled(X, y, status, left, right),
                        w_obs, maxit, tol)
  if (model == "pooling") {
    res <- pooled
    w_units <- w_obs
  } else {
    w_units <- if (is.null(weights)) rep(1, length(group_size)) else individual_weights(weights, group_size)
    keep_ind <- w_units > 0
    keep_obs <- rep(keep_ind, group_size)
    group <- rep(seq_len(sum(keep_ind)), group_size[keep_ind])
    contrib <- tobit_contrib_random(X[keep_obs, , drop = FALSE], y[keep_obs], status[keep_obs],
                                    left, right, group, R)
    # split the pooled variance between the two components
    start <- c(pooled$theta[seq_len(k)], pooled$theta[k + 1L] + log(sqrt(0.6)),
               pooled$theta[k + 1L] + log(sqrt(0.4)))
    res <- ml_optimize(start, contrib, w_units[keep_ind], maxit, tol)
    score <- matrix(0, length(w_units), length(res$theta))
    score[keep_ind, ] <- res$score
    res$score <- score
  }
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
