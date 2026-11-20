#' Pooled GLM by IRLS in R, with an offset and step-halving
#'
#' Used for the pooled fits the compiled IRLS does not cover: a nonzero
#' offset (any family) and the log-binomial model. Each step is halved until
#' the linear index is admissible (fitted risks below 1 for the log-binomial
#' model) and the log-likelihood does not decrease, so the iterates stay in
#' the parameter space. NB2 alternates IRLS for beta with a one-dimensional
#' maximization in theta, as MASS::glm.nb.
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
  conv <- FALSE
  for (it in seq_len(maxit)) {
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
            " IRLS iterations",
            if (logbin) paste0("; the largest fitted risk is ", format(max(mu[pos]), digits = 6),
                               ". Steps are halved to keep the risks below 1, which slows the ",
                               "iterations when they approach 1: raise maxit (e.g. 1000). If a ",
                               "risk tends to 1 the estimate lies on the boundary of the ",
                               "parameter space; the modified Poisson model (family = ",
                               "\"poisson\", vcov = \"HC1\") estimates the same risk ratios"),
            call. = FALSE)
  }
  ww <- w * mu_eta_r(eta, lid)^2 / vfun(mu)
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
       theta = if (nb) theta, df.residual = n_pos - k - nb, iterations = it, converged = conv)
}
