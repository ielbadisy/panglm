#' Beta regression contributions
#'
#' `y ~ Beta(mu * phi, (1 - mu) * phi)` with `g(mu) = eta`, so that
#' `E(y) = mu` and `Var(y) = mu (1 - mu) / (1 + phi)`. The parameter vector
#' is `theta = (beta, log phi)` for the pooled model and
#' `(beta, log phi, log sigma)` with a Gaussian random intercept on the
#' link scale.
#'
#' @keywords internal
#' @noRd
beta_link <- function(link) {
  if (link == "logit") {
    list(mu = stats::plogis,
         d1 = function(eta) { m <- stats::plogis(eta); m * (1 - m) },
         d2 = function(eta) { m <- stats::plogis(eta); m * (1 - m) * (1 - 2 * m) })
  } else {
    list(mu = stats::pnorm,
         d1 = stats::dnorm,
         d2 = function(eta) -eta * stats::dnorm(eta))
  }
}

beta_obs_ll <- function(mu, phi, ly, l1y) {
  mu <- pmin(pmax(mu, 1e-12), 1 - 1e-12)
  lgamma(phi) - lgamma(mu * phi) - lgamma((1 - mu) * phi) +
    (mu * phi - 1) * ly + ((1 - mu) * phi - 1) * l1y
}

beta_contrib_pooled <- function(X, y, link) {
  k <- ncol(X)
  lf <- beta_link(link)
  ly <- log(y); l1y <- log1p(-y); ystar <- ly - l1y
  function(theta, score) {
    eta <- as.numeric(X %*% theta[seq_len(k)])
    if (!is.finite(theta[k + 1L]) || abs(theta[k + 1L]) > 25) {
      return(list(ll = rep(-Inf, length(y)), score = matrix(0, length(y), k + 1L)))
    }
    phi <- exp(theta[k + 1L])
    mu <- pmin(pmax(lf$mu(eta), 1e-12), 1 - 1e-12)
    out <- list(ll = beta_obs_ll(mu, phi, ly, l1y))
    if (score) {
      p1 <- digamma(mu * phi); p2 <- digamma((1 - mu) * phi)
      dmu <- phi * (ystar - (p1 - p2))
      dphi <- phi * (digamma(phi) - mu * p1 - (1 - mu) * p2 + mu * ly + (1 - mu) * l1y)
      out$score <- cbind((dmu * lf$d1(eta)) * X, dphi)
    }
    out
  }
}

#' Random-intercept beta contributions by adaptive Gauss-Hermite quadrature
#'
#' The mode of each individual's integrand is found by Newton steps on the
#' observed curvature, with a Fisher scoring step (always a negative
#' curvature) for the individuals where the observed curvature is not
#' negative, starting from the modes of the previous evaluation; the
#' expected curvature at the mode sets the scale of the quadrature rule.
#' Returns one log-likelihood per individual.
#'
#' The second derivative of the beta log-density with respect to `mu` does
#' not depend on `y`, so the observed curvature on the link scale is the
#' expected one plus `score_mu * g''`, at no extra cost.
#'
#' @keywords internal
#' @noRd
beta_contrib_random <- function(X, y, link, group, R) {
  k <- ncol(X); G <- max(group)
  lf <- beta_link(link)
  ly <- log(y); l1y <- log1p(-y); ystar <- ly - l1y
  gh <- gauss_hermite_quadrature(R)
  log_w <- log(gh$weights * sqrt(pi)) + gh$nodes^2
  state <- new.env(parent = emptyenv())
  function(theta, score) {
    eta <- as.numeric(X %*% theta[seq_len(k)])
    if (any(!is.finite(theta[k + 1:2])) || any(abs(theta[k + 1:2]) > 25)) {
      return(list(ll = rep(-Inf, G)))
    }
    phi <- exp(theta[k + 1L])
    sigma <- exp(theta[k + 2L])
    h <- function(u) {
      mu <- lf$mu(eta + u[group])
      as.numeric(rowsum(beta_obs_ll(mu, phi, ly, l1y), group, reorder = FALSE)) +
        stats::dnorm(u, 0, sigma, log = TRUE)
    }
    curvature <- function(u) {
      e <- eta + u[group]
      mu <- pmin(pmax(lf$mu(e), 1e-12), 1 - 1e-12)
      d1 <- lf$d1(e)
      a <- mu * phi; b <- (1 - mu) * phi
      score_mu <- phi * (ystar - (digamma(a) - digamma(b)))
      info_mu <- phi^2 * (trigamma(a) + trigamma(b))
      sums <- rowsum(cbind(score_mu * d1, info_mu * d1^2, score_mu * lf$d2(e)), group,
                     reorder = FALSE)
      H <- -sums[, 2L] - 1 / sigma^2
      list(g = sums[, 1L] - u / sigma^2, H = H, Hobs = H + sums[, 3L])
    }
    find_mode <- function(u) {
      for (step in seq_len(100L)) {
        cv <- curvature(u)
        Hs <- ifelse(cv$Hobs < 0, cv$Hobs, cv$H)
        delta <- pmax(-2, pmin(2, cv$g / Hs))
        if (anyNA(delta)) return(NULL)
        u <- u - delta
        if (max(abs(delta)) < 1e-9) return(u)
      }
      u
    }
    u <- NULL
    if (length(state$u) == G && all(is.finite(state$u))) u <- find_mode(state$u)
    if (is.null(u)) u <- find_mode(rep(0, G))
    if (is.null(u)) return(list(ll = rep(-Inf, G)))
    state$u <- u
    tau <- 1 / sqrt(-curvature(u)$H)
    terms <- vapply(seq_along(gh$nodes), function(q) {
      log_w[q] + h(u + sqrt(2) * tau * gh$nodes[q])
    }, numeric(G))
    if (!is.matrix(terms)) terms <- matrix(terms, nrow = G)
    list(ll = log(sqrt(2) * tau) + log_row_sum_exp(terms))
  }
}

#' Fit pooled or random-intercept beta regression
#'
#' @keywords internal
#' @noRd
fit_beta <- function(X, y, family, model, group_start, group_size, weights, R, maxit, tol) {
  k <- ncol(X)
  lf <- beta_link(family$link)
  w_obs <- weights %||% rep(1, length(y))
  # starting values: least squares on the link scale, moment estimate of phi
  ylink <- if (family$link == "logit") stats::qlogis(y) else stats::qnorm(y)
  b0 <- as.numeric(stats::lm.wfit(X, ylink, w_obs)$coefficients)
  b0[!is.finite(b0)] <- 0
  mu0 <- lf$mu(as.numeric(X %*% b0))
  phi0 <- max(1, stats::weighted.mean(mu0 * (1 - mu0) / pmax((y - mu0)^2, 1e-8), w_obs) - 1)
  phi0 <- min(phi0, 1e3)
  # the optimizer works on standardized columns: with a badly scaled
  # regressor the first quasi-Newton steps otherwise push the linear
  # predictor to where mu is 0 or 1 and the likelihood is flat
  pooled <- ml_optimize(c(b0, log(phi0)), beta_contrib_pooled(X, y, family$link), w_obs,
                        maxit, tol, scaling = ml_design_scaling(X, k + 1L),
                        warn = model == "pooling", what = "the beta regression fit")
  if (model == "pooling") {
    res <- pooled
    w_units <- w_obs
  } else {
    w_units <- if (is.null(weights)) rep(1, length(group_size)) else individual_weights(weights, group_size)
    keep_ind <- w_units > 0
    keep_obs <- rep(keep_ind, group_size)
    group <- rep(seq_len(sum(keep_ind)), group_size[keep_ind])
    contrib <- beta_contrib_random(X[keep_obs, , drop = FALSE], y[keep_obs], family$link, group, R)
    start <- c(pooled$theta, log(0.5))
    start[k + 1L] <- start[k + 1L] + log(1.5)
    res <- ml_optimize(start, contrib, w_units[keep_ind], maxit, tol,
                       scaling = ml_design_scaling(X, k + 2L),
                       what = "the random-effects beta regression fit")
    score <- matrix(0, length(w_units), length(res$theta))
    score[keep_ind, ] <- res$score
    res$score <- score
  }
  coefs <- res$theta[seq_len(k)]
  names(coefs) <- colnames(X)
  vcov <- res$vcov[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(vcov) <- list(colnames(X), colnames(X))
  phi <- exp(res$theta[k + 1L])
  sigma <- if (model == "random") exp(res$theta[k + 2L]) else NULL
  eta <- as.numeric(X %*% coefs)
  fitted <- if (model == "random") {
    random_binomial_marginal_mean(eta, sigma, family, R)
  } else {
    lf$mu(eta)
  }
  list(coefficients = coefs, vcov = vcov, bread = vcov, fitted.values = fitted,
       linear.predictors = eta, loglik = res$loglik, precision = phi, sigma = sigma,
       dispersion = 1 / (1 + phi), theta = res$theta, theta_vcov = res$vcov,
       theta_score = res$score, unit_weights = w_units,
       quadrature_nodes = if (model == "random") R else NULL,
       npar = length(res$theta), df.residual = sum(w_obs > 0) - length(res$theta),
       iterations = res$iterations, convergence = res$convergence)
}

#' Sandwich covariance for beta regression fits
#'
#' @keywords internal
#' @noRd
robust_vcov_beta <- function(object, type, cluster) {
  k <- length(object$coefficients)
  if (object$model == "random") return(robust_vcov_random(object, cluster))
  cl <- if (type == "cluster") panglm_cluster_vector(object, cluster) else NULL
  v <- ml_sandwich(object$theta_vcov, object$theta_score, object$unit_weights, cl)
  v <- v[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(v) <- dimnames(object$vcov)
  v
}
