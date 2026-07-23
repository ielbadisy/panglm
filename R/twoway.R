#' Two-way (individual + time) demeaning via alternating projections
#'
#' Iteratively demeans by `id` then by `time` until both group means are
#' simultaneously ~0 (the standard Gauss-Seidel / method-of-alternating-
#' projections algorithm for exact two-way fixed-effects demeaning; see
#' Guimaraes & Portugal 2010). Exact for both balanced and unbalanced
#' panels. The group-sum reduction and the demeaning step are parallelized
#' across observations in C++ (RcppParallel), the same hot-path treatment
#' as the rest of the package's estimators.
#'
#' @keywords internal
#' @noRd
demean_twoway <- function(X, y, group_start, group_size, time, tol = 1e-10, maxit = 10000) {
  # id groups are already contiguous/sorted (from build_panel_index), so no
  # factor()/match() is needed to derive their codes; only the time
  # dimension needs a grouping pass, done via data.table's fast dense rank.
  n_id <- length(group_size)
  id_code <- rep.int(seq_len(n_id) - 1L, group_size)
  time_code <- as.integer(data.table::frankv(time, ties.method = "dense") - 1L)
  n_time <- max(time_code) + 1L

  M <- cbind(y, X)
  dm <- twoway_demean_cpp(M, id_code, time_code, n_id, n_time, maxit, tol)

  list(y = dm$M[, 1], X = dm$M[, -1, drop = FALSE], iterations = dm$iterations, n_time = n_time)
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

  list(coefficients = coefs, vcov = vcov, bread = bread, fitted.values = NULL,
       loglik = NA_real_, dispersion = sigma2, df.residual = df_resid,
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
#' `eta_new = X_tilde \%*\% beta_new + (z - z_tilde)` (the demeaned-away part
#' of z, `z - z_tilde`, is exactly the fitted two-way effect contribution).
#'
#' @keywords internal
#' @noRd
fit_within_twoways_poisson <- function(X, y, group_start, group_size, time, maxit, tol) {
  if (is.null(time)) stop("effect = 'twoways' requires a two-element index = c(id, time)", call. = FALSE)

  n <- nrow(X); k <- ncol(X)
  n_id <- length(group_size)
  id_code <- rep.int(seq_len(n_id) - 1L, group_size)
  time_code <- as.integer(data.table::frankv(time, ties.method = "dense") - 1L)
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
  X_tilde <- dm$M[, -1, drop = FALSE]
  XtWX <- crossprod(X_tilde, X_tilde * w)

  coefs <- beta / scale
  names(coefs) <- colnames(X)
  vcov_unscaled <- solve(XtWX)
  vcov <- vcov_unscaled / outer(scale, scale)
  dimnames(vcov) <- list(colnames(X), colnames(X))

  df_resid <- n - k - (n_id - 1) - (n_time - 1) - 1

  list(coefficients = coefs, vcov = vcov, bread = vcov, fitted.values = exp(eta),
       loglik = ll_old, dispersion = 1, df.residual = df_resid,
       iterations = iter, demean_iterations = demean_iters)
}
