#' Two-way (individual + time) demeaning via alternating projections
#'
#' Iteratively demeans by `id` then by `time` (the Gauss-Seidel /
#' method-of-alternating-projections algorithm for two-way fixed-effects
#' demeaning; see Guimaraes & Portugal 2010) until the largest change in
#' one sweep is below `tol`. The iterations converge to the two-way
#' demeaned values for balanced and unbalanced panels (in one sweep for
#' balanced ones); a warning is given if `maxit` sweeps are reached
#' without convergence. The group-sum reduction and the demeaning step are parallelized
#' across observations in C++ (RcppParallel), the same hot-path treatment
#' as the rest of the package's estimators.
#'
#' @keywords internal
#' @noRd
demean_twoway <- function(X, y, group_start, group_size, time, tol = 1e-10, maxit = 10000) {
  # id groups are already contiguous/sorted (from build_panel_index), so no
  # factor()/match() is needed to derive their codes; only the time
  # dimension needs a grouping pass, done via dense_rank0().
  n_id <- length(group_size)
  id_code <- rep.int(seq_len(n_id) - 1L, group_size)
  time_code <- dense_rank0(time)
  n_time <- max(time_code) + 1L

  M <- cbind(y, X)
  dm <- twoway_demean_cpp(M, id_code, time_code, n_id, n_time, maxit, tol)
  twoway_demean_check(dm, maxit)

  list(y = dm$M[, 1], X = dm$M[, -1, drop = FALSE], iterations = dm$iterations, n_time = n_time)
}

# Warn when the alternating projections stop at the iteration limit.
twoway_demean_check <- function(dm, maxit) {
  if (isFALSE(dm$converged)) {
    warning("two-way demeaning (alternating projections) did not converge in ", maxit,
            " iterations; the estimates may be inaccurate", call. = FALSE)
  }
  invisible(dm)
}

fit_within_twoways_gaussian <- function(X, y, group_start, group_size, time, maxit, tol) {
  if (is.null(time)) stop("effect = 'twoways' requires a two-element index = c(id, time)", call. = FALSE)

  dm <- demean_twoway(X, y, group_start, group_size, time, tol = tol)
  res <- irls_fit_cpp(dm$X, dm$y, 0L, 0L, maxit, tol)

  n <- nrow(X); k <- ncol(X)
  n_id <- length(group_size); n_time <- dm$n_time
  # two-way within df: n - k - (n_id - 1) - (n_time - 1) - 1, matching plm's
  # model = "within", effect = "twoways" convention
  df_resid <- n - k - (n_id - 1) - (n_time - 1) - 1

  resid <- dm$y - dm$X %*% res$coefficients
  sigma2 <- sum(resid^2) / max(1, df_resid)
  vcov <- res$vcov_unscaled * sigma2
  bread <- res$vcov_unscaled
  coefs <- as.numeric(res$coefficients)
  names(coefs) <- colnames(X)
  dimnames(vcov) <- dimnames(bread) <- list(colnames(X), colnames(X))
  effects <- recover_twoway_effects(
    y - as.numeric(X %*% coefs), group_size, time, tol = tol
  )
  fitted <- as.numeric(X %*% coefs) + effects$observation_effect

  ## Same Frisch-Waugh-Lovell argument as the one-way within estimator (see
  ## fit_within() in panglm.R): the two-way alternating-projection residuals
  ## are algebraically identical to the full individual+time dummy-variable
  ## regression's residuals, so the standard OLS ML Gaussian log-likelihood
  ## applies at n = total observations (matches panglm_parameter_count()'s
  ## npar = k + (n_id - 1) + (n_time - 1) + 1 for this branch).
  full_resid <- y - fitted
  rss <- sum(full_resid^2)
  loglik <- if (rss > 0) -0.5 * n * (log(2 * pi) + log(rss / n) + 1) else NA_real_

  list(coefficients = coefs, vcov = vcov, bread = bread, fitted.values = fitted,
       individual_effects = effects$individual,
       time_effects = effects$time,
       loglik = loglik, dispersion = sigma2, df.residual = df_resid,
       iterations = res$iterations, demean_iterations = dm$iterations)
}

#' Two-way (individual + time) fixed-effects Poisson via the standard
#' outer-IRLS / inner-weighted-FWL algorithm
#'
#' Poisson fixed effects are multiplicative on the mean scale, not additive
#' like gaussian, so they can't be profiled out by one demeaning pass. This
#' is the algorithm behind `fixest::feglm()`/`alpaca`/`ppmlhdfe`: each IRLS
#' iteration computes the working response `z`, partials out both fixed
#' effects from `z` and `X` via *weighted* (weight = IRLS weight = mu)
#' alternating projections, fits beta by weighted least squares on the
#' partialled-out variables, then reconstructs the full linear predictor as
#' `eta_new = X_tilde %*% beta_new + (z - z_tilde)` (the demeaned-away part
#' of z, `z - z_tilde`, is exactly the fitted two-way effect contribution).
#'
#' @keywords internal
#' @noRd
fit_within_twoways_poisson <- function(X, y, group_start, group_size, time, maxit, tol) {
  if (is.null(time)) stop("effect = 'twoways' requires a two-element index = c(id, time)", call. = FALSE)

  n <- nrow(X); k <- ncol(X)
  n_id <- length(group_size)
  id_code <- rep.int(seq_len(n_id) - 1L, group_size)
  time_code <- dense_rank0(time)
  n_time <- max(time_code) + 1L

  scale <- apply(X, 2, function(col) { s <- stats::sd(col); if (s > 0) s else 1 })
  Xs <- sweep(X, 2, scale, "/")

  eta <- rep(0, n)
  beta <- rep(0, k)
  ll_old <- -Inf
  XtWX <- diag(k)
  demean_iters <- 0L
  iter <- 0L

  for (iter in seq_len(maxit)) {
    eta <- pmin(pmax(eta, -30), 30)
    mu <- exp(eta)
    w <- mu
    z <- eta + (y - mu) / mu

    M <- cbind(z, Xs)
    dm <- twoway_demean_weighted_cpp(M, w, id_code, time_code, n_id, n_time, 10000, tol)
    twoway_demean_check(dm, 10000)
    demean_iters <- dm$iterations
    z_tilde <- dm$M[, 1]
    X_tilde <- dm$M[, -1, drop = FALSE]

    Xw <- X_tilde * w
    XtWX <- crossprod(X_tilde, Xw)
    XtWz <- crossprod(Xw, z_tilde)
    beta_new <- as.numeric(solve(XtWX, XtWz))

    eta_new <- as.numeric(X_tilde %*% beta_new) + (z - z_tilde)
    eta_new <- pmin(pmax(eta_new, -30), 30)
    mu_new <- exp(eta_new)
    ll_new <- sum(y * eta_new - mu_new - lgamma(y + 1))

    beta <- beta_new
    eta <- eta_new
    converged <- abs(ll_new - ll_old) < tol * (abs(ll_old) + 1)
    ll_old <- ll_new
    if (converged) break
  }

  # Refresh the Fisher information at the converged eta/mu -- XtWX above is
  # from the demeaning computed with the *previous* iteration's weights, one
  # step behind the final beta, which understates precision slightly.
  mu <- exp(eta)
  w <- mu
  z <- eta + (y - mu) / mu
  M <- cbind(z, Xs)
  dm <- twoway_demean_weighted_cpp(M, w, id_code, time_code, n_id, n_time, 10000, tol)
  twoway_demean_check(dm, 10000)
  X_tilde <- dm$M[, -1, drop = FALSE]
  XtWX <- crossprod(X_tilde, X_tilde * w)

  coefs <- beta / scale
  names(coefs) <- colnames(X)
  vcov_unscaled <- solve(XtWX)
  vcov <- vcov_unscaled / outer(scale, scale)
  dimnames(vcov) <- list(colnames(X), colnames(X))

  df_resid <- n - k - (n_id - 1) - (n_time - 1) - 1
  effects <- recover_twoway_effects(
    eta - as.numeric(Xs %*% beta), group_size, time, tol = tol
  )

  list(coefficients = coefs, vcov = vcov, bread = vcov, fitted.values = exp(eta),
       individual_effects = effects$individual,
       time_effects = effects$time,
       # regressors net of both effects (weighted by mu): the robust score is
       # (y - mu) * X_tilde, which profiles out the time effects as well
       X_tilde = sweep(X_tilde, 2, scale, "*"),
       loglik = ll_old, dispersion = 1, df.residual = df_resid,
       iterations = iter, demean_iterations = demean_iters)
}

#' Two-way (individual + time) fixed-effects NB2 via the same outer-IRLS /
#' inner-weighted-FWL algorithm used for Poisson twoways
#'
#' NB2's IRLS differs from Poisson's only in the working weight (`mu*theta /
#' (theta+mu)` instead of `mu`, from the NB2 variance function `mu +
#' mu^2/theta`) and the working response formula is unchanged (`eta + (y -
#' mu)/mu`); the outer loop and the two-way weighted-demeaning step carry
#' over unchanged from `fit_within_twoways_poisson()`. `theta` is profiled
#' out via a 1-D optimization of the NB2 log-likelihood at each outer
#' iteration's converged mu (this model has no independent reference
#' implementation to match exactly, unlike the alternating IRLS-for-beta /
#' Newton-for-theta scheme used for pooled/one-way negbin, which matches
#' `MASS::glm.nb`).
#'
#' @keywords internal
#' @noRd
fit_within_twoways_negbin <- function(X, y, group_start, group_size, time, maxit, tol) {
  if (is.null(time)) stop("effect = 'twoways' requires a two-element index = c(id, time)", call. = FALSE)

  n <- nrow(X); k <- ncol(X)
  n_id <- length(group_size)
  id_code <- rep.int(seq_len(n_id) - 1L, group_size)
  time_code <- dense_rank0(time)
  n_time <- max(time_code) + 1L

  scale <- apply(X, 2, function(col) { s <- stats::sd(col); if (s > 0) s else 1 })
  Xs <- sweep(X, 2, scale, "/")

  nb2_loglik <- function(y, mu, theta) {
    sum(lgamma(y + theta) - lgamma(theta) - lgamma(y + 1) +
        theta * log(theta) - theta * log(theta + mu) +
        y * log(mu) - y * log(theta + mu))
  }

  eta <- rep(0, n)
  beta <- rep(0, k)
  ybar <- mean(y); yvar <- stats::var(y)
  theta <- if (yvar > ybar && ybar > 0) ybar^2 / (yvar - ybar) else 100
  if (!is.finite(theta) || theta <= 0) theta <- 100

  ll_old <- -Inf
  XtWX <- diag(k)
  demean_iters <- 0L
  iter <- 0L

  for (iter in seq_len(maxit)) {
    eta <- pmin(pmax(eta, -30), 30)
    mu <- exp(eta)
    w <- (mu * theta) / (theta + mu)
    z <- eta + (y - mu) / mu

    M <- cbind(z, Xs)
    dm <- twoway_demean_weighted_cpp(M, w, id_code, time_code, n_id, n_time, 10000, tol)
    twoway_demean_check(dm, 10000)
    demean_iters <- dm$iterations
    z_tilde <- dm$M[, 1]
    X_tilde <- dm$M[, -1, drop = FALSE]

    Xw <- X_tilde * w
    XtWX <- crossprod(X_tilde, Xw)
    XtWz <- crossprod(Xw, z_tilde)
    beta_new <- as.numeric(solve(XtWX, XtWz))

    eta_new <- as.numeric(X_tilde %*% beta_new) + (z - z_tilde)
    eta_new <- pmin(pmax(eta_new, -30), 30)
    mu_new <- exp(eta_new)

    theta_new <- tryCatch({
      opt <- stats::optimize(function(lk) nb2_loglik(y, mu_new, exp(lk)),
                              interval = c(log(1e-4), log(1e6)), maximum = TRUE, tol = tol)
      exp(opt$maximum)
    }, error = function(e) theta)

    ll_new <- nb2_loglik(y, mu_new, theta_new)

    beta <- beta_new
    eta <- eta_new
    theta <- theta_new
    converged <- abs(ll_new - ll_old) < tol * (abs(ll_old) + 1)
    ll_old <- ll_new
    if (converged) break
  }

  # Refresh the Fisher information at the converged eta/mu/theta, same
  # rationale as fit_within_twoways_poisson()'s post-loop refresh.
  mu <- exp(eta)
  w <- (mu * theta) / (theta + mu)
  z <- eta + (y - mu) / mu
  M <- cbind(z, Xs)
  dm <- twoway_demean_weighted_cpp(M, w, id_code, time_code, n_id, n_time, 10000, tol)
  twoway_demean_check(dm, 10000)
  X_tilde <- dm$M[, -1, drop = FALSE]
  XtWX <- crossprod(X_tilde, X_tilde * w)

  coefs <- beta / scale
  names(coefs) <- colnames(X)
  vcov_unscaled <- solve(XtWX)
  vcov <- vcov_unscaled / outer(scale, scale)
  dimnames(vcov) <- list(colnames(X), colnames(X))

  df_resid <- n - k - (n_id - 1) - (n_time - 1) - 1
  effects <- recover_twoway_effects(
    eta - as.numeric(Xs %*% beta), group_size, time, tol = tol
  )

  list(coefficients = coefs, vcov = vcov, bread = vcov, fitted.values = exp(eta),
       individual_effects = effects$individual,
       time_effects = effects$time,
       X_tilde = sweep(X_tilde, 2, scale, "*"),
       loglik = ll_old, theta = theta, dispersion = 1, df.residual = df_resid,
       iterations = iter, demean_iterations = demean_iters)
}

recover_twoway_effects <- function(effect_observation, group_size, time,
                                   tol = 1e-10, maxit = 10000L) {
  id <- rep.int(seq_along(group_size), group_size)
  time_code <- match(time, sort(unique(time)))
  n_id <- length(group_size)
  n_time <- max(time_code)
  alpha <- numeric(n_id)
  gamma <- numeric(n_time)

  for (iteration in seq_len(maxit)) {
    alpha_new <- as.numeric(rowsum(effect_observation - gamma[time_code], id)) /
      group_size
    gamma_new <- as.numeric(rowsum(effect_observation - alpha_new[id], time_code)) /
      as.numeric(tabulate(time_code, nbins = n_time))
    shift <- mean(gamma_new)
    gamma_new <- gamma_new - shift
    alpha_new <- alpha_new + shift
    change <- max(abs(c(alpha_new - alpha, gamma_new - gamma)))
    alpha <- alpha_new
    gamma <- gamma_new
    if (change < tol) break
  }

  list(
    individual = alpha,
    time = gamma,
    observation_effect = alpha[id] + gamma[time_code]
  )
}

#' Two-way fixed-effects conditional (Chamberlain) logit
#'
#' The individual effects are conditioned out exactly, as in
#' `fit_within_binomial()`, and the time effects enter the conditional
#' likelihood as time dummies (first time level as reference), which is
#' `survival::clogit(y ~ x + factor(time) + strata(id), method = "exact")`.
#' The time dummies are ordinary parameters of the conditional likelihood:
#' the covariance of the slopes is the slope block of the inverse of the
#' full (slopes + time dummies) information, and the robust covariances
#' use the full score and bread (see robust_vcov_within_binomial()).
#'
#' The likelihood is maximized by Newton-Raphson with the analytic gradient
#' and Hessian of the exact conditional likelihood (the recursion of
#' survival::clogit(method = "exact")).
#'
#' @keywords internal
#' @noRd
fit_within_twoways_binomial <- function(X, y, group_start, group_size, time, maxit, tol) {
  if (is.null(time)) stop("effect = 'twoways' requires a two-element index = c(id, time)", call. = FALSE)

  k <- ncol(X)
  time_levels <- sort(unique(time))
  time_code <- match(time, time_levels)
  n_time <- length(time_levels)
  D <- outer(time_code, seq_len(n_time), "==") * 1
  D <- D[, -1L, drop = FALSE]
  time_names <- paste0("time", time_levels[-1L])
  colnames(D) <- time_names

  # Time dummies identified among the informative strata (outcome varies
  # within the individual), net of the individual effects. A dummy for a
  # wave observed only in uninformative strata, for instance, is dropped
  # (its effect is reported as NA). The time dummies come first in the
  # rank check, so that a slope collinear with the time effects among the
  # informative strata is detected rather than a time dummy dropped in its
  # place.
  panel <- rep.int(seq_along(group_size), group_size)
  successes <- as.numeric(rowsum(y, panel))
  informative <- (successes > 0 & successes < group_size)[panel]
  if (!any(informative)) {
    stop("no individual has within-individual outcome variation; the conditional ",
         "logit is not identified", call. = FALSE)
  }
  Z <- cbind(D, X)
  Zd <- within_demean_cpp(Z, y, group_start, group_size)$X[informative, , drop = FALSE]
  candidate <- which(sqrt(colSums(Zd^2)) > 1e-7 * sqrt(colSums(Z[informative, , drop = FALSE]^2)))
  decomposition <- qr(Zd[, candidate, drop = FALSE], tol = 1e-7)
  kept <- candidate[decomposition$pivot[seq_len(decomposition$rank)]]
  kept_time <- sort(kept[kept <= ncol(D)])
  if (!all(ncol(D) + seq_len(k) %in% kept)) {
    stop("a slope is not identified separately from the time effects among the ",
         "individuals whose outcome varies", call. = FALSE)
  }
  Z <- cbind(X, D[, kept_time, drop = FALSE])
  p <- ncol(Z)

  nr <- conditional_logit_newton(Z, y, group_start, group_size, maxit, tol)
  beta <- nr$beta
  H <- -nr$hessian
  opt <- list(value = -nr$loglik, convergence = nr$convergence, counts = c(nr$iterations, NA))
  newton_steps <- nr$iterations

  info <- conditional_logit_loglik_grad_cpp(beta, Z, y, group_start, group_size)
  n_dropped <- length(group_size) - info$n_used_groups
  if (n_dropped > 0) {
    message(n_dropped, " group(s) with no within-group outcome variation dropped ",
            "(all-0 or all-1) -- they carry no information for the conditional likelihood.")
  }

  vcov_full <- tryCatch(solve(H), error = function(e) matrix(NA_real_, p, p))
  dimnames(vcov_full) <- list(colnames(Z), colnames(Z))
  slope <- seq_len(k)
  coefs <- beta[slope]
  names(coefs) <- colnames(X)
  vcov <- vcov_full[slope, slope, drop = FALSE]
  dimnames(vcov) <- list(colnames(X), colnames(X))

  time_effects <- stats::setNames(rep(NA_real_, n_time), as.character(time_levels))
  time_effects[1L] <- 0
  estimated_time <- colnames(Z)[-slope]
  time_effects[match(estimated_time, time_names) + 1L] <- beta[-slope]
  time_se <- stats::setNames(rep(NA_real_, n_time), as.character(time_levels))
  time_se[1L] <- 0
  time_se[match(estimated_time, time_names) + 1L] <- sqrt(diag(vcov_full)[-slope])

  conditional_fitted <- y - as.numeric(info$score_eta)
  list(coefficients = coefs, vcov = vcov, bread = vcov, fitted.values = conditional_fitted,
       time_effects = time_effects, time_effects_se = time_se,
       # full conditional-likelihood design and bread (slopes + time dummies)
       # for the robust covariances
       X_aug = Z, vcov_full_aug = vcov_full,
       loglik = nr$loglik, dispersion = 1, df.residual = NA_real_,
       npar = p,
       n_used_groups = info$n_used_groups, n_dropped_groups = n_dropped,
       iterations = newton_steps, convergence = opt$convergence)
}
