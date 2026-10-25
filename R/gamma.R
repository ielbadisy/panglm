#' Fixed-effects gamma regression with a log link
#'
#' Maximum likelihood with one intercept per individual, computed by IRLS.
#' Under the log link the gamma working weights are all one, so each IRLS
#' step is a within (demeaned) least-squares regression of the working
#' response, and the individual intercepts follow in closed form. With a
#' multiplicative individual effect, the slope estimator is consistent for
#' fixed T under strict exogeneity (the concentrated score has mean zero),
#' like the Poisson fixed-effects estimator.
#'
#' @keywords internal
#' @noRd
fit_within_gamma <- function(X, y, group_start, group_size, maxit, tol) {
  n <- nrow(X); k <- ncol(X); G <- length(group_size)
  group <- rep(seq_len(G), group_size)
  beta <- rep(0, k)
  alpha <- log(as.numeric(rowsum(y, group)) / group_size)
  eta <- alpha[group]
  deviance <- function(mu) 2 * sum(-log(y / mu) + (y - mu) / mu)
  dev_old <- deviance(exp(eta))
  iter <- 0L
  for (iter in seq_len(maxit)) {
    mu <- exp(eta)
    z <- eta + (y - mu) / mu
    dm <- within_demean_cpp(X, z, group_start, group_size)
    beta_new <- as.numeric(solve(crossprod(dm$X), crossprod(dm$X, dm$y)))
    alpha_new <- as.numeric(rowsum(z - X %*% beta_new, group)) / group_size
    eta_new <- as.numeric(X %*% beta_new) + alpha_new[group]
    dev_new <- deviance(exp(eta_new))
    halving <- 0L
    while (!is.finite(dev_new) || dev_new > dev_old + 1e-12 * abs(dev_old)) {
      halving <- halving + 1L
      if (halving > 30L) break
      beta_new <- (beta + beta_new) / 2
      alpha_new <- (alpha + alpha_new) / 2
      eta_new <- as.numeric(X %*% beta_new) + alpha_new[group]
      dev_new <- deviance(exp(eta_new))
    }
    converged <- abs(dev_new - dev_old) / (abs(dev_old) + 0.1) < tol
    beta <- beta_new; alpha <- alpha_new; eta <- eta_new; dev_old <- dev_new
    if (converged) break
  }
  mu <- exp(eta)
  Xd <- within_demean_cpp(X, y, group_start, group_size)$X
  bread <- solve(crossprod(Xd))
  df_resid <- n - k - G
  dispersion <- sum(((y - mu) / mu)^2) / max(1, df_resid)
  vcov <- bread * dispersion
  names(beta) <- colnames(X)
  dimnames(vcov) <- dimnames(bread) <- list(colnames(X), colnames(X))
  list(coefficients = beta, vcov = vcov, bread = bread, fitted.values = mu,
       individual_effects = alpha, loglik = gamma_loglik(y, mu),
       dispersion = dispersion, df.residual = df_resid, iterations = iter)
}

#' Random-intercept gamma regression with a log link
#'
#' `y_it | u_i ~ Gamma(shape = nu, mean = exp(x_it'beta + u_i))` with
#' `u_i ~ N(0, sigma^2)`. The individual integral is computed by adaptive
#' Gauss-Hermite quadrature centred at the mode of each individual's
#' integrand (`R = 1` is the Laplace approximation). The integrand only
#' depends on the data through four individual sums, so each evaluation
#' costs O(N + G R).
#'
#' @keywords internal
#' @noRd
fit_random_gamma <- function(X, y, group_start, group_size, R, maxit, tol) {
  k <- ncol(X); G <- length(group_size)
  group <- rep(seq_len(G), group_size)
  gh <- gauss_hermite_quadrature(R)
  log_w <- log(gh$weights * sqrt(pi)) + gh$nodes^2
  Ti <- group_size
  slogy <- as.numeric(rowsum(log(y), group))

  loglik_at <- function(theta) {
    beta <- theta[seq_len(k)]
    sigma <- exp(theta[k + 1L])
    nu <- exp(theta[k + 2L])
    eta <- as.numeric(X %*% beta)
    seta <- as.numeric(rowsum(eta, group))
    sy <- as.numeric(rowsum(y * exp(-eta), group))
    h <- function(u) {
      Ti * (nu * log(nu) - lgamma(nu)) + (nu - 1) * slogy -
        nu * (seta + Ti * u) - nu * sy * exp(-u) +
        stats::dnorm(u, 0, sigma, log = TRUE)
    }
    # Newton for the mode of the (strictly concave) integrand in u_i
    u <- rep(0, G)
    for (step in seq_len(100L)) {
      g1 <- -nu * Ti + nu * sy * exp(-u) - u / sigma^2
      g2 <- -nu * sy * exp(-u) - 1 / sigma^2
      delta <- pmax(-5, pmin(5, g1 / g2))
      if (anyNA(delta)) return(-Inf)
      u <- u - delta
      if (max(abs(delta)) < 1e-10) break
    }
    tau <- 1 / sqrt(nu * sy * exp(-u) + 1 / sigma^2)
    terms <- vapply(seq_along(gh$nodes), function(q) {
      log_w[q] + h(u + sqrt(2) * tau * gh$nodes[q])
    }, numeric(G))
    if (!is.matrix(terms)) terms <- matrix(terms, nrow = G)
    top <- apply(terms, 1, max)
    sum(log(sqrt(2) * tau) + top + log(rowSums(exp(terms - top))))
  }
  objective <- function(theta) {
    value <- -loglik_at(theta)
    if (!is.finite(value)) 1e100 else value
  }
  gradient <- function(theta) {
    vapply(seq_along(theta), function(j) {
      d <- 1e-5 * max(1, abs(theta[j]))
      up <- theta; up[j] <- up[j] + d
      dn <- theta; dn[j] <- dn[j] - d
      (objective(up) - objective(dn)) / (2 * d)
    }, numeric(1))
  }

  pooled <- irls_fit_cpp(X, y, 3L, 1L, maxit, tol)
  beta0 <- as.numeric(pooled$coefficients)
  resid_mean <- as.numeric(rowsum(log(y) - X %*% beta0, group)) / Ti
  sigma0 <- max(0.05, stats::sd(resid_mean))
  nu0 <- 1 / max(1e-3, pooled$dispersion)
  theta0 <- c(beta0, log(sigma0), log(nu0))

  opt <- stats::optim(theta0, objective, gradient, method = "BFGS",
                      control = list(maxit = max(200L, maxit), reltol = 1e-14))
  hessian <- stats::optimHess(opt$par, objective, gradient)
  vcov_full <- tryCatch(solve(hessian), error = function(e) matrix(NA_real_, k + 2, k + 2))

  coefs <- opt$par[seq_len(k)]
  names(coefs) <- colnames(X)
  sigma <- exp(opt$par[k + 1L])
  shape <- exp(opt$par[k + 2L])
  vcov <- vcov_full[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(vcov) <- list(colnames(X), colnames(X))
  eta <- as.numeric(X %*% coefs)
  list(coefficients = coefs, vcov = vcov,
       fitted.values = exp(eta + sigma^2 / 2),
       loglik = -opt$value, sigma = sigma, shape = shape, dispersion = 1 / shape,
       quadrature_nodes = R, df.residual = nrow(X) - k - 2,
       iterations = opt$counts[[1]], convergence = opt$convergence)
}
