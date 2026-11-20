#' Pooled GLM by IRLS in R, with an offset and step-halving
#'
#' Used for the pooled fits the compiled IRLS does not cover: a nonzero
#' offset (any family) and the log-binomial model. Each step is halved until
#' the log-likelihood does not decrease. NB2 alternates IRLS for beta with a
#' one-dimensional maximization in theta, as MASS::glm.nb. The log-binomial
#' model is fitted by logbin_newton().
#'
#' @keywords internal
#' @noRd
fit_pooled_r <- function(X, y, family, off, w, maxit, tol) {
  fam <- family$family; lid <- family$link_id
  w <- w %||% rep(1, length(y)); pos <- w > 0; k <- ncol(X)
  nb <- fam == "negbin"; logbin <- fam == "binomial" && lid == 1L
  theta <- if (nb) 1 else Inf
  vfun <- switch(fam, gaussian = function(mu) rep(1, length(mu)), poisson = function(mu) mu,
                 binomial = function(mu) mu * (1 - mu), gamma = function(mu) mu^2,
                 negbin = function(mu) mu + mu^2 / theta)
  llfun <- function(mu) sum((w * switch(fam,
    gaussian = -0.5 * (y - mu)^2, poisson = y * log(mu) - mu - lgamma(y + 1),
    binomial = y * log(mu) + (1 - y) * log1p(-mu), gamma = -(y / mu + log(mu)),
    negbin = stats::dnbinom(y, size = theta, mu = mu, log = TRUE)))[pos])
  ok <- function(eta) all(is.finite(eta[pos])) && (!logbin || all(eta[pos] < 0))
  linkinv <- function(eta) linkinv_r(eta, lid)
  # start: least squares of the linearized mean, or the intercept-only fit
  mu0 <- switch(fam, binomial = (y + 0.5) / 2, gaussian = y, pmax(y, 0.1))
  beta <- stats::lm.wfit(X, linkfun_r(mu0, lid) - off, w)$coefficients
  beta[is.na(beta)] <- 0
  if (!ok(X %*% beta + off)) {
    beta <- numeric(k); int <- match("(Intercept)", colnames(X))
    b0 <- log(sum(w * y) / sum(w * exp(off)))
    if (!is.na(int)) beta[int] <- min(b0, -max(off[pos]) - 1e-4)
    if (!ok(X %*% beta + off)) {
      stop("no admissible starting value for the log-binomial model: every fitted risk ",
           "must be below 1", call. = FALSE)
    }
  }
  eta <- as.numeric(X %*% beta + off); mu <- linkinv(eta); ll <- llfun(mu)
  conv <- FALSE; nbound <- 0L
  if (logbin) {
    nf <- logbin_newton(X, y, w, off, beta, maxit)
    beta <- nf$beta; eta <- nf$eta; mu <- linkinv(eta); it <- nf$it; conv <- nf$conv
    nbound <- nf$nbound; ll <- nf$ll
  } else for (it in seq_len(maxit)) {
    if (nb) {
      theta <- exp(stats::optimize(function(lt) { theta <<- exp(lt); llfun(mu) }, c(-10, 15),
                                   maximum = TRUE, tol = 1e-12)$maximum)
      ll <- llfun(mu)
    }
    d <- mu_eta_r(eta, lid)
    ww <- w * d^2 / vfun(mu)
    step <- stats::lm.wfit(X, eta - off + (y - mu) / d, ww)$coefficients - beta
    step[is.na(step)] <- 0
    for (h in 0:40) {
      b_new <- beta + step / 2^h
      e_new <- as.numeric(X %*% b_new + off)
      if (ok(e_new) && (ll_new <- llfun(linkinv(e_new))) >= ll - 1e-12 * abs(ll)) break
    }
    if (h == 40L) break
    beta <- b_new; eta <- e_new; mu <- linkinv(eta)
    # a full step that is small, not a halved one, signals convergence
    done <- h == 0L && abs(ll_new - ll) / (abs(ll) + 0.1) < tol &&
      max(abs(step)) < 1e-8 * (1 + max(abs(beta)))
    ll <- ll_new
    if (done && it > 1L) { conv <- TRUE; break }
  }
  if (!conv) {
    warning("the pooled ", fam, " (", family$link, " link) fit did not converge in ", it,
            if (logbin) " Newton" else " IRLS", " iterations",
            if (logbin) paste0("; the largest fitted risk is ", format(max(mu[pos]), digits = 6)),
            call. = FALSE)
  }
  if (conv && nbound > 0L) {
    warning("the log-binomial fit reached its maximum on the parameter boundary (fitted risk at 1) ",
            "for ", nbound, " observation", if (nbound > 1L) "s", ": the estimates are the constrained ",
            "maximum, but the standard errors assume an interior maximum and are not reliable. ",
            "The modified Poisson model (family = \"poisson\", vcov = \"HC1\") estimates the ",
            "risk ratios without that constraint", call. = FALSE)
  }
  # at the boundary the expected information is infinite (risk 1): use the observed one there
  if (logbin) mu <- pmin(mu, 1 - 1e-10)
  ww <- if (nbound > 0L) w * (1 - y) * mu / (1 - mu)^2 else w * mu_eta_r(eta, lid)^2 / vfun(mu)
  bread <- chol2inv(chol(crossprod(X * sqrt(ww))))
  n_pos <- sum(pos); r <- (y - mu)[pos]
  phi <- switch(fam, gaussian = sum((w * (y - mu)^2)[pos]) / max(1, n_pos - k),
                gamma = sum(w[pos] * (r / mu[pos])^2) / max(1, n_pos - k), 1)
  if (fam == "gaussian") {
    rss <- sum((w * (y - mu)^2)[pos])
    ll <- -0.5 * n_pos * (log(2 * pi * rss / n_pos) + 1) + 0.5 * sum(log(w[pos]))
  }
  if (fam == "gamma") ll <- gamma_loglik(y, mu, w)
  names(beta) <- colnames(X)
  dimnames(bread) <- list(colnames(X), colnames(X))
  list(coefficients = beta, vcov = bread * phi, bread = bread, fitted.values = mu,
       linear.predictors = eta, loglik = ll, dispersion = phi,
       theta = if (nb) theta, df.residual = n_pos - k - nb, iterations = it, converged = conv,
       boundary = if (logbin) nbound > 0L)
}

#' Log-binomial maximum likelihood by active-set Newton
#'
#' Maximizes the log-binomial log-likelihood, which is concave in beta, under
#' the constraints eta_i <= 0 (fitted risk at most 1); only observations with
#' y = 1 can reach the bound, since a risk of 1 for y = 0 has likelihood 0.
#' Newton steps use the exact Hessian and stay in the null space of the
#' active constraints, are cut at the first blocking constraint, which then
#' becomes active, and are halved until the log-likelihood does not
#' decrease. At a stationary point an active constraint with a negative
#' Lagrange multiplier is released. Without active constraints this is
#' plain Newton-Raphson, so an interior maximum is the glm() estimate.
#' Returns nbound, the number of observations with a fitted risk of 1.
#'
#' @keywords internal
#' @noRd
logbin_newton <- function(X, y, w, off, beta, maxit) {
  pos <- w > 0; one <- pos & y == 1; zer <- pos & y < 1; k <- ncol(X); n <- length(y)
  ll_of <- function(e) if (any(e[zer] >= 0) || any(e[one] > 1e-10)) -Inf else
    sum((w * y * pmin(e, 0))[pos]) + sum((w * (1 - y) * log1p(-exp(e)))[zer])
  eta <- as.numeric(X %*% beta + off); ll <- ll_of(eta)
  act <- released <- integer(); conv <- FALSE
  for (it in seq_len(maxit)) {
    mu <- exp(eta); g <- hw <- numeric(n)
    g[pos] <- (w * y)[pos]; g[zer] <- g[zer] - (w * (1 - y) * mu / (1 - mu))[zer]
    hw[zer] <- (w * (1 - y) * mu / (1 - mu)^2)[zer]
    g <- crossprod(X, g)[, 1]
    act <- union(act, setdiff(which(one & eta > -1e-9), released)); released <- integer()
    Z <- if (length(act)) {
      q <- qr(t(X[act, , drop = FALSE]))
      qr.Q(q, complete = TRUE)[, -seq_len(q$rank), drop = FALSE]
    } else diag(k)
    d <- if (ncol(Z)) {
      XZ <- (X %*% Z) * sqrt(hw)
      as.numeric(Z %*% qr.solve(crossprod(XZ), crossprod(Z, g), tol = 1e-12))
    } else numeric(k)
    if (max(abs(d)) < 1e-8 * (1 + max(abs(beta)))) {
      if (!length(act)) { conv <- TRUE; break }
      lam <- qr.coef(qr(t(X[act, , drop = FALSE])), g); lam[is.na(lam)] <- 0
      if (all(lam >= -1e-6 * (1 + max(abs(g))))) { conv <- TRUE; break }
      released <- act[which.min(lam)]; act <- setdiff(act, released)
      next
    }
    xd <- as.numeric(X %*% d); blk <- one & xd > 0 & !(seq_len(n) %in% act)
    t <- min(1, if (any(blk)) min(-eta[blk] / xd[blk]))
    for (h in 0:50) {
      e_new <- pmin(as.numeric(X %*% (beta + t * d) + off), ifelse(one, 0, Inf))
      if ((ll_new <- ll_of(e_new)) >= ll - 1e-12 * abs(ll)) break
      t <- t / 2
    }
    if (h == 50L) break
    beta <- beta + t * d; eta <- e_new; ll <- ll_new
  }
  list(beta = beta, eta = eta, ll = ll, it = it, conv = conv, nbound = sum(one & eta > -1e-9))
}
