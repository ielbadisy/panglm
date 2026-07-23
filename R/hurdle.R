#' Fixed-effects hurdle model for panel count data
#'
#' A hurdle decomposition for panels where a large fraction of groups report
#' *structural* zeros (never observed to occur), rather than zeros arising
#' from ordinary Poisson/NB sampling variation at a low mean -- a pattern
#' that shows up as extreme, family-invariant overdispersion under either
#' `family = "poisson"` or `family = "negbin"` in [panglm()].
#'
#' The model is fit in two independent parts, both using fixed effects, and
#' reported together rather than folded into `panglm()`'s single-family
#' dispatch, since they operate on different sample sizes (the count part
#' drops the zeros):
#'
#' \enumerate{
#'   \item **Zero vs. positive.** `1(y > 0)` fit via the exact conditional
#'     logistic regression already used for `panglm(..., family =
#'     "binomial")` (Chamberlain 1980) -- no new estimator, just a new
#'     outcome vector.
#'   \item **Count given `y > 0`.** A zero-truncated Poisson fixed-effects
#'     model (Allison-Waterman-style dummy-variable intercepts, the same
#'     unconditional approach used for FE-NB2) on the subsample where
#'     `y > 0`. The log-density is the usual Poisson log-density minus
#'     `log(1 - P(Y = 0))`; see e.g. `pscl::hurdle()` for the reference form
#'     of the truncated-count log-likelihood.
#' }
#'
#' Groups with no remaining within-group variation after truncation (a
#' single surviving observation) contribute no information about the shared
#' count-part slope, the same profiling fact that applies to the other
#' fixed-effects estimators in this package; they are not dropped
#' explicitly since a perfectly-fit singleton doesn't diverge, it simply
#' carries zero weight.
#'
#' @param formula a model formula, e.g. `y ~ x1 + x2`
#' @param data a data.frame containing the variables in `formula` and `index`
#' @param index length-2 character vector identifying the individual and
#'   time columns, e.g. `c("id", "time")`
#' @param maxit maximum IRLS/Newton iterations
#' @param tol convergence tolerance
#' @return an object of class `"panglm_hurdle"`: a list with `zero` (the
#'   `panglm` binomial fit for the zero-vs-positive part) and `count` (the
#'   truncated-Poisson fixed-effects fit for the count part)
#' @export
panglm_hurdle <- function(formula, data, index, maxit = 100, tol = 1e-10) {
  if (missing(index) || length(index) != 2) {
    stop("'index' must be a two-element character vector, e.g. index = c(\"id\", \"time\")", call. = FALSE)
  }

  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  X_full <- stats::model.matrix(formula, data = mf)
  has_intercept <- "(Intercept)" %in% colnames(X_full)
  X <- if (has_intercept) X_full[, setdiff(colnames(X_full), "(Intercept)"), drop = FALSE] else X_full

  panel_cols <- data[index]
  panel_cols <- panel_cols[match(rownames(mf), rownames(data)), , drop = FALSE]
  panel <- build_panel_index(cbind(panel_cols, .panglm_row = seq_len(nrow(mf))), index)
  ord <- panel$data$.panglm_row

  y <- y[ord]
  X <- X[ord, , drop = FALSE]
  group_start <- panel$group_start
  group_size <- panel$group_size

  y_bin <- as.numeric(y > 0)
  zero_fit <- fit_within_binomial(X, y_bin, group_start, group_size, maxit, tol)
  zero_fit$formula <- formula; zero_fit$index <- index
  class(zero_fit) <- "panglm"

  group <- rep(seq_along(group_size), group_size)
  keep <- y > 0
  count_fit <- fit_truncated_poisson_dummy(X[keep, , drop = FALSE], y[keep], group[keep], maxit, tol)

  out <- list(zero = zero_fit, count = count_fit, formula = formula, index = index,
              n = length(y), n_positive = sum(keep))
  class(out) <- "panglm_hurdle"
  out
}

#' Zero-truncated Poisson via the Allison-Waterman-style unconditional
#' dummy-variable estimator
#'
#' Same dummy-augmentation trick as [fit_within_negbin_dummy()] (append one
#' column per surviving group, jointly estimate covariate slopes and
#' per-group intercepts), but for the zero-truncated Poisson log-density
#' instead of NB2, fit by Newton-Raphson with the exact Fisher information
#' (the score of a zero-truncated Poisson w.r.t. the linear predictor has
#' the same GLM-like form as ordinary Poisson, `y - mu_trunc`, where
#' `mu_trunc = mu / (1 - exp(-mu))` is the truncated-Poisson mean).
#'
#' @keywords internal
#' @noRd
fit_truncated_poisson_dummy <- function(X, y, group, maxit, tol) {
  kept_groups <- sort(unique(group))
  group <- match(group, kept_groups)
  n <- length(y); k <- ncol(X); G <- length(kept_groups)
  colnames_X <- colnames(X)

  dummies <- matrix(0, n, G)
  dummies[cbind(seq_len(n), group)] <- 1
  Xa <- cbind(X, dummies)
  p <- k + G

  beta <- c(rep(0, k), log(pmax(as.numeric(rowsum(y, group)) / as.numeric(table(group)), 1e-2)))
  ll_old <- -Inf
  iter <- 0L
  XtVX <- diag(p)

  for (iter in seq_len(maxit)) {
    eta <- as.numeric(Xa %*% beta)
    eta <- pmin(pmax(eta, -30), 30)
    mu <- exp(eta)
    # mu / (1 - exp(-mu)), computed via expm1 to avoid cancellation as mu -> 0
    mu_trunc <- mu / (-expm1(-mu))
    V <- mu_trunc * (1 + mu - mu_trunc)
    V <- pmax(V, 1e-8)

    score <- as.numeric(crossprod(Xa, y - mu_trunc))
    XtVX <- crossprod(Xa, Xa * V) + 1e-10 * diag(p)
    step <- solve(XtVX, score)
    beta <- beta + step

    eta_new <- pmin(pmax(as.numeric(Xa %*% beta), -30), 30)
    mu_new <- exp(eta_new)
    ll_new <- sum(y * eta_new - mu_new - lgamma(y + 1) - log(-expm1(-mu_new)))

    converged <- abs(ll_new - ll_old) < tol * (abs(ll_old) + 1)
    ll_old <- ll_new
    if (converged) break
  }

  coefs <- beta[seq_len(k)]
  alpha_i <- beta[k + seq_len(G)]
  names(coefs) <- colnames_X

  vcov_full <- solve(XtVX)
  vcov <- vcov_full[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(vcov) <- list(colnames_X, colnames_X)

  eta <- pmin(pmax(as.numeric(Xa %*% beta), -30), 30)
  mu <- exp(eta)
  fitted <- mu / (-expm1(-mu))

  list(coefficients = coefs, vcov = vcov, bread = vcov, fitted.values = fitted,
       loglik = ll_old, dispersion = 1,
       individual_effects = stats::setNames(alpha_i, kept_groups),
       n_used_groups = G, df.residual = n - k - G, iterations = iter)
}

#' @export
print.panglm_hurdle <- function(x, ...) {
  cat("\nFixed-effects hurdle model for panel count data\n")
  cat("Zero part (1(y>0), conditional logit):\n")
  print(stats::coef(x$zero))
  cat("\nCount part (zero-truncated Poisson, N =", x$n_positive, "positive obs of", x$n, "):\n")
  print(x$count$coefficients)
  invisible(x)
}
