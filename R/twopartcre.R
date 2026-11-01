#' Correlated random-effects two-part model: individual log-likelihoods
#'
#' Part 1: `P(y_it > 0 | u1) = F(x'gamma + u1)` (probit or logit).
#' Part 2: `log y_it | y_it > 0, u2 ~ N(z'beta + u2, sigma_e^2)`.
#' `(u1, u2)` bivariate normal with standard deviations `sigma1`, `sigma2`
#' and correlation `rho`. Given `u1`, `u2 ~ N(c u1, tau^2)` with
#' `c = rho sigma2 / sigma1` and `tau^2 = sigma2^2 (1 - rho^2)`, so the
#' positive part integrates over `u2` in closed form (a Gaussian random
#' intercept likelihood with shifted residuals); the remaining integral over
#' `u1` uses adaptive Gauss-Hermite quadrature around the exact mode (the
#' integrand is log-concave in `u1`). The likelihood is on the cost scale
#' (it includes the Jacobian `-sum(log y)`). The contributions and their
#' exact gradient (including the dependence of the mode and the curvature
#' on the parameters) are computed in C++ (`src/twopart_cre.cpp`).
#'
#' `theta = (gamma, beta, log sigma1, log sigma2, log sigma_e, atanh rho)`.
#'
#' @keywords internal
#' @noRd
twopart_cre_contrib <- function(Xb, d, group_b, Xp, ly, group_p, link, R) {
  kb <- ncol(Xb); kp <- ncol(Xp); G <- max(group_b)
  gh <- gauss_hermite_quadrature(R)
  log_w <- log(gh$weights * sqrt(pi)) + gh$nodes^2
  storage.mode(Xb) <- "double"; storage.mode(Xp) <- "double"
  sgn <- as.numeric(2 * d - 1)
  ly <- as.numeric(ly)
  group_b <- as.integer(group_b); group_p <- as.integer(group_p)
  probit <- link == "probit"
  function(theta, score) {
    lim <- theta[kb + kp + 1:4]
    if (any(!is.finite(lim)) || any(abs(lim[1:3]) > 15) || abs(lim[4]) > 7) {
      return(list(ll = rep(-Inf, G)))
    }
    out <- twopart_cre_cpp(Xb, sgn, group_b, Xp, ly, group_p, G, as.numeric(theta),
                           probit, gh$nodes, log_w, isTRUE(score))
    if (score && is.null(out$score)) {
      out$score <- matrix(0, G, length(theta))
    }
    out
  }
}

#' Joint fit of the correlated random-effects two-part model
#'
#' Starts from the separate random-effects fits (rho = 0) and replaces their
#' coefficients and covariances by the joint estimates.
#'
#' @keywords internal
#' @noRd
fit_twopart_cre <- function(fit_binary, fit_positive, index, link, vcov_type, R,
                            maxit, tol) {
  id_b <- fit_binary$data[[index[1]]]
  id_p <- fit_positive$data[[index[1]]]
  ids <- unique(id_b)
  group_b <- match(id_b, ids)
  group_p <- match(id_p, ids)
  Xb <- fit_binary$X; Xp <- fit_positive$X
  kb <- ncol(Xb); kp <- ncol(Xp)
  w_units <- if (is.null(fit_binary$weights)) rep(1, length(ids)) else {
    individual_weights(fit_binary$weights, fit_binary$group_size)
  }
  contrib <- twopart_cre_contrib(Xb, fit_binary$y, group_b, Xp, fit_positive$y, group_p,
                                 link, R)
  sb <- fit_binary$sigma
  sv <- fit_positive$sigma_v2; su <- max(fit_positive$sigma_mu2, 1e-4)
  start <- c(unname(fit_binary$coefficients), unname(fit_positive$coefficients),
             log(max(sb, 0.05)), 0.5 * log(su), 0.5 * log(sv), 0)
  res <- ml_optimize(start, contrib, w_units, maxit, tol)
  p <- length(res$theta)
  V <- res$vcov
  if (vcov_type == "cluster") {
    cl <- individual_clusters(fit_binary$cluster_id, fit_binary$group_size)
    V <- ml_sandwich(res$vcov, res$score, w_units, cl, clusters_only = TRUE)
  }
  nb <- c(paste0("binary:", names(fit_binary$coefficients)),
          paste0("positive:", names(fit_positive$coefficients)))
  extra <- res$theta[kb + kp + 1:4]
  # delta method for (sigma1, sigma2, sigma_e, rho)
  J <- diag(c(exp(extra[1:3]), 1 - tanh(extra[4])^2))
  Vx <- V[kb + kp + 1:4, kb + kp + 1:4]
  list(theta = res$theta, vcov = V[seq_len(kb + kp), seq_len(kb + kp), drop = FALSE],
       names = nb, kb = kb, kp = kp,
       sigma_binary = exp(extra[1]), sigma_positive = exp(extra[2]),
       sigma_e = exp(extra[3]), rho = tanh(extra[4]),
       variance_params = stats::setNames(c(exp(extra[1:3]), tanh(extra[4])),
                                         c("sigma_u1", "sigma_u2", "sigma_e", "rho")),
       variance_vcov = J %*% Vx %*% J,
       loglik = res$loglik, npar = p, score = res$score, theta_vcov = res$vcov,
       unit_weights = w_units, convergence = res$convergence, quadrature_nodes = R)
}

#' Population-averaged means of the correlated two-part model
#'
#' `E[y] = exp(z'beta + sigma_e^2 / 2 + tau^2 / 2) E_u1[F(x'gamma + u1) exp(c u1)]`
#' and `P(y > 0) = E_u1[F(x'gamma + u1)]`, by Gauss-Hermite quadrature.
#'
#' @keywords internal
#' @noRd
twopart_cre_means <- function(a, b, joint, link, R = 41L) {
  gh <- gauss_hermite_quadrature(R)
  s1 <- joint$sigma_binary; s2 <- joint$sigma_positive; rho <- joint$rho
  cc <- rho * s2 / s1; tau2 <- s2^2 * (1 - rho^2)
  Fn <- if (link == "probit") stats::pnorm else stats::plogis
  u <- sqrt(2) * s1 * gh$nodes
  wq <- gh$weights  # normalised: sum to 1
  P <- vapply(u, function(uq) Fn(a + uq), numeric(length(a)))
  if (!is.matrix(P)) P <- matrix(P, nrow = length(a))
  prob <- as.numeric(P %*% wq)
  total <- exp(b + joint$sigma_e^2 / 2 + tau2 / 2) * as.numeric(P %*% (wq * exp(cc * u)))
  list(prob = prob, total = total)
}
