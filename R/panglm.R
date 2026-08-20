#' Fit a generalized linear model for panel data
#'
#' Fast pooled, fixed-effects ("within") and random-effects panel GLMs,
#' backed by an Rcpp/RcppArmadillo/RcppParallel numerical core. See the
#' package vignette for the estimators currently implemented per family.
#'
#' @param formula a model formula, e.g. `y ~ x1 + x2`
#' @param data a data.frame containing the variables in `formula` and `index`
#' @param index length-1 or length-2 character vector identifying the
#'   individual (and optionally time) columns, e.g. `c("id", "time")`
#' @param model one of `"pooling"`, `"within"`, `"random"`
#' @param family one of `"gaussian"`, `"poisson"`, `"binomial"`, `"negbin"`,
#'   or a family spec from [gaussian_family()] / [poisson_family()] /
#'   [binomial_family()] / [negbin_family()]
#' @param effect one of `"individual"` or `"twoways"` (individual + time
#'   fixed effects). `"twoways"` is currently only implemented for
#'   `model = "within"`, `family = "gaussian"` (exact alternating-projections
#'   demeaning), `"poisson"`, or `"negbin"` (both via outer-IRLS /
#'   inner-weighted-FWL, the algorithm behind `fixest::feglm()`); not yet for
#'   `"binomial"` (no general closed-form two-way conditional logit exists).
#' @param vcov one of `"classical"`, `"HC1"` (heteroskedasticity-robust) or
#'   `"cluster"` (cluster-robust, clustered by the panel individual by
#'   default). Only applies to `model = "pooling"`/`"within"`; random-effects
#'   models always report model-based (information-matrix) standard errors.
#'   See [vcov.panglm()] to recompute a different type after fitting.
#' @param R number of Gauss-Hermite quadrature nodes (random-effects
#'   binomial models only)
#' @param maxit maximum IRLS/Newton iterations
#' @param tol convergence tolerance
#' @return an object of class `"panglm"`
#' @examples
#' data(copd)
#'
#' # pooled OLS on lung function (ignores the panel structure)
#' fit_pool <- panglm(fev1 ~ treatment + age + smoker + crp, data = copd,
#'                     index = c("id", "visit"), model = "pooling", family = "gaussian")
#' summary(fit_pool)
#'
#' # fixed-effects ("within") estimator -- only the time-varying regressor
#' # (crp) survives demeaning; time-invariant covariates like treatment
#' # or age are not identified by a fixed-effects estimator
#' fit_fe <- panglm(fev1 ~ crp, data = copd,
#'                   index = c("id", "visit"), model = "within", family = "gaussian")
#' coef(fit_fe)
#'
#' # random-effects (Swamy-Arora) estimator
#' fit_re <- panglm(fev1 ~ treatment + age + smoker + crp, data = copd,
#'                   index = c("id", "visit"), model = "random", family = "gaussian")
#' coef(fit_re)
#'
#' # fixed-effects Poisson on exacerbation counts
#' fit_pois <- panglm(exacerbations ~ crp, data = copd,
#'                     index = c("id", "visit"), model = "within", family = "poisson")
#' coef(fit_pois)
#' @export
panglm <- function(formula, data, index,
                    model = c("pooling", "within", "random"),
                    family = "gaussian",
                    effect = c("individual", "twoways"),
                    vcov = c("classical", "HC1", "cluster"),
                    R = 21, maxit = 100, tol = 1e-10) {
  model <- match.arg(model)
  effect <- match.arg(effect)
  vcov <- match.arg(vcov)
  family <- resolve_family(family)

  if (missing(index)) stop("'index' is required, e.g. index = c(\"id\", \"time\")", call. = FALSE)
  if (effect == "twoways" && !(model == "within" && family$family %in% c("gaussian", "poisson", "negbin"))) {
    stop("effect = 'twoways' is currently only implemented for ",
         "model = 'within', family = 'gaussian', 'poisson', or 'negbin'", call. = FALSE)
  }

  mf <- stats::model.frame(formula, data = data)
  y <- stats::model.response(mf)
  if (family$family == "binomial") y <- binomial_response_to_numeric(y)
  X_full <- stats::model.matrix(formula, data = mf)

  panel_cols <- data[index]
  panel_cols <- panel_cols[match(rownames(mf), rownames(data)), , drop = FALSE]
  panel <- build_panel_index(cbind(panel_cols, .panglm_row = seq_len(nrow(mf))), index)
  ord <- panel$data$.panglm_row

  y <- y[ord]
  X_full <- X_full[ord, , drop = FALSE]
  group_start <- panel$group_start
  group_size <- panel$group_size
  time_id <- if (length(index) > 1) panel$data[[index[2]]] else NULL

  has_intercept <- "(Intercept)" %in% colnames(X_full)
  X_noint <- if (has_intercept) X_full[, setdiff(colnames(X_full), "(Intercept)"), drop = FALSE] else X_full

  fit <- switch(model,
    pooling = fit_pooled(X_full, y, family, maxit, tol),
    within  = if (effect == "twoways" && family$family == "gaussian") {
      fit_within_twoways_gaussian(X_noint, y, group_start, group_size, time_id, maxit, tol)
    } else if (effect == "twoways" && family$family == "poisson") {
      fit_within_twoways_poisson(X_noint, y, group_start, group_size, time_id, maxit, tol)
    } else if (effect == "twoways" && family$family == "negbin") {
      fit_within_twoways_negbin(X_noint, y, group_start, group_size, time_id, maxit, tol)
    } else {
      fit_within(X_noint, y, family, group_start, group_size, maxit, tol)
    },
    random  = fit_random(X_full, y, family, group_start, group_size, R, maxit, tol)
  )

  fit$X <- if (model == "pooling" || model == "random") X_full else X_noint
  fit$y <- y
  fit$group_start <- group_start
  fit$group_size <- group_size
  fit$time_id <- time_id
  fit$cluster_id <- rep(panel$group_id, group_size)
  fit$call <- match.call()
  fit$formula <- formula
  fit$model <- model
  fit$effect <- effect
  fit$family <- family
  fit$index <- index
  fit$nobs <- length(y)
  fit$n_groups <- panel$n_groups
  fit$vcov_type <- "classical"
  class(fit) <- "panglm"

  if (vcov != "classical") {
    fit$vcov <- vcov.panglm(fit, type = vcov)
    fit$vcov_type <- vcov
  }
  fit
}

fit_pooled <- function(X, y, family, maxit, tol) {
  if (family$family == "negbin") {
    res <- negbin_irls_fit_cpp(X, y, maxit, tol)
    coefs <- as.numeric(res$coefficients)
    names(coefs) <- colnames(X)
    vcov <- res$vcov_unscaled
    dimnames(vcov) <- list(colnames(X), colnames(X))
    fitted <- as.numeric(exp(X %*% res$coefficients))
    return(list(coefficients = coefs, vcov = vcov, bread = vcov, fitted.values = fitted,
                loglik = res$loglik, theta = res$theta, dispersion = 1,
                df.residual = nrow(X) - ncol(X) - 1, iterations = res$iterations))
  }
  res <- irls_fit_cpp(X, y, family$family_id, family$link_id, maxit, tol)
  vcov <- res$vcov_unscaled * res$dispersion
  bread <- res$vcov_unscaled
  coefs <- as.numeric(res$coefficients)
  names(coefs) <- colnames(X)
  dimnames(vcov) <- dimnames(bread) <- list(colnames(X), colnames(X))
  list(coefficients = coefs, vcov = vcov, bread = bread, fitted.values = as.numeric(res$fitted.values),
       loglik = res$loglik, dispersion = res$dispersion, df.residual = nrow(X) - ncol(X),
       iterations = res$iterations)
}

fit_within <- function(X, y, family, group_start, group_size, maxit, tol) {
  if (family$family == "gaussian") {
    dm <- within_demean_cpp(X, y, group_start, group_size)
    res <- irls_fit_cpp(dm$X, dm$y, family$family_id, family$link_id, maxit, tol)
    n <- nrow(X); k <- ncol(X); G <- length(group_start)
    df_resid <- n - G - k
    resid <- dm$y - dm$X %*% res$coefficients
    sigma2 <- sum(resid^2) / max(1, df_resid)
    vcov <- res$vcov_unscaled * sigma2
    bread <- res$vcov_unscaled
    coefs <- as.numeric(res$coefficients)
    names(coefs) <- colnames(X)
    dimnames(vcov) <- dimnames(bread) <- list(colnames(X), colnames(X))
    gm <- group_means_cpp(X, y, group_start, group_size)
    alpha_i <- gm$ybar - as.numeric(gm$Xbar %*% res$coefficients)
    alpha_obs <- rep(alpha_i, group_size)
    fitted <- alpha_obs + as.numeric(X %*% res$coefficients)
    return(list(coefficients = coefs, vcov = vcov, bread = bread, fitted.values = fitted,
                loglik = NA_real_, dispersion = sigma2, df.residual = df_resid,
                iterations = res$iterations))
  }
  if (family$family == "poisson") {
    scale <- apply(X, 2, function(col) { s <- stats::sd(col); if (s > 0) s else 1 })
    Xs <- sweep(X, 2, scale, "/")
    res <- within_poisson_fit_cpp(Xs, y, group_start, group_size, maxit, tol)
    coefs <- as.numeric(res$coefficients) / scale
    names(coefs) <- colnames(X)
    vcov <- res$vcov_unscaled / outer(scale, scale)
    dimnames(vcov) <- list(colnames(X), colnames(X))

    # Closed-form group intercept for the conditional Poisson MLE:
    # alpha_i = log(sum(y_i) / sum(exp(x_it'beta))), so mu_it = (Yi/Li) * exp(x_it'beta).
    eta <- as.numeric(X %*% coefs)
    lit <- exp(eta)
    group <- rep(seq_along(group_size), group_size)
    Li <- as.numeric(rowsum(lit, group))
    Yi <- as.numeric(rowsum(y, group))
    fitted <- (Yi / Li)[group] * lit

    return(list(coefficients = coefs, vcov = vcov, bread = vcov, fitted.values = fitted,
                loglik = res$loglik, dispersion = 1, df.residual = nrow(X) - ncol(X) - length(group_start),
                iterations = res$iterations))
  }
  if (family$family == "negbin") return(fit_within_negbin_dummy(X, y, group_start, group_size, maxit, tol))
  if (family$family == "binomial") return(fit_within_binomial(X, y, group_start, group_size, maxit, tol))
  stop("model = 'within' is not yet implemented for family '", family$family, "'", call. = FALSE)
}

#' Allison-Waterman (2002) unconditional fixed-effects negative binomial
#'
#' Augments the design matrix with one dummy column per panel individual
#' and jointly estimates (covariate slopes, individual intercepts, NB2
#' dispersion) by exact NB2 MLE. Unlike the conditional (Hausman-Hall-
#' Griliches) FE-NB estimator, this does not suffer the "fake fixed
#' effects" critique (Guimaraes 2008): the individual intercepts are
#' estimated directly, not conditioned out via a shared dispersion trick.
#' The incidental-parameters problem that would bias a fixed-effects
#' *logit* dummy-variable estimator does not carry over to NB2 (Allison &
#' Waterman 2002), so this is consistent for moderate N (tractable up to
#' a few hundred individuals; the dummy design matrix grows with N, so it
#' is not intended for very large panels -- see [panglm()]'s vignette).
#'
#' @keywords internal
#' @noRd
fit_within_negbin_dummy <- function(X, y, group_start, group_size, maxit, tol) {
  n0 <- nrow(X); k <- ncol(X); G0 <- length(group_size)
  group <- rep(seq_along(group_size), group_size)
  colnames_X <- colnames(X)

  # Groups with an all-zero outcome have an unbounded dummy-variable MLE
  # (the intercept diverges to -Inf under the log link), the same boundary
  # case fit_within_binomial() screens out for the conditional logit. Left
  # in, the divergent direction can contaminate the shared beta via the
  # joint Newton/IRLS solve, not just that group's own intercept.
  zero_group <- vapply(seq_len(G0), function(g) all(y[group == g] == 0), logical(1))
  n_dropped <- sum(zero_group)
  kept_groups <- which(!zero_group)
  keep_rows <- !zero_group[group]
  if (n_dropped > 0) {
    message(n_dropped, " group(s) with all-zero outcomes dropped -- their fixed effect ",
            "is unbounded under the dummy-variable NB2 MLE and carries no information ",
            "about the shared slope.")
    X <- X[keep_rows, , drop = FALSE]
    y <- y[keep_rows]
    group <- match(group[keep_rows], kept_groups)
  }
  n <- length(y); G <- length(kept_groups)

  dummies <- matrix(0, n, G)
  dummies[cbind(seq_len(n), group)] <- 1
  X_aug <- cbind(X, dummies)

  res <- negbin_irls_fit_cpp(X_aug, y, maxit, tol)
  coefs <- as.numeric(res$coefficients[seq_len(k)])
  alpha_i <- as.numeric(res$coefficients[k + seq_len(G)])
  names(coefs) <- colnames_X

  vcov_full <- res$vcov_unscaled
  vcov <- vcov_full[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(vcov) <- list(colnames_X, colnames_X)

  # Pad back to the original (pre-screening) row count with NA for dropped
  # rows, so fitted.values stays aligned with the full-length y stored on
  # the panglm object (residuals.panglm() does object$y - fitted.values).
  fitted <- rep(NA_real_, n0)
  fitted[keep_rows] <- as.numeric(exp(X_aug %*% res$coefficients))

  list(coefficients = coefs, vcov = vcov, bread = vcov, fitted.values = fitted,
       loglik = res$loglik, theta = res$theta, dispersion = 1,
       individual_effects = stats::setNames(alpha_i, kept_groups),
       n_used_groups = G, n_dropped_groups = n_dropped,
       # Kept for robust/cluster vcov (see robust_vcov_within_negbin() in
       # vcov.R): the full (k+G)-parameter augmented design/bread, and
       # keep_rows to align a full-length cluster/offset vector to it.
       X_aug = X_aug, keep_rows = keep_rows, vcov_full_aug = vcov_full,
       df.residual = n - k - G - 1, iterations = res$iterations)
}

fit_within_binomial <- function(X, y, group_start, group_size, maxit, tol) {
  k <- ncol(X)

  # BFGS calls the objective and gradient separately at (usually) the same
  # beta; cache the last evaluation so the O(sum T_i^2) DP isn't run twice
  # per iteration.
  cache <- new.env(parent = emptyenv())
  eval_at <- function(beta) {
    if (is.null(cache$beta) || !identical(cache$beta, beta)) {
      cache$beta <- beta
      cache$res <- conditional_logit_loglik_grad_cpp(beta, X, y, group_start, group_size)
    }
    cache$res
  }
  objective <- function(beta) -eval_at(beta)$loglik
  gradient <- function(beta) -eval_at(beta)$gradient

  opt <- stats::optim(rep(0, k), objective, gradient, method = "BFGS",
                       control = list(maxit = maxit, reltol = tol), hessian = TRUE)

  info <- conditional_logit_loglik_grad_cpp(opt$par, X, y, group_start, group_size)
  n_dropped <- length(group_size) - info$n_used_groups
  if (n_dropped > 0) {
    message(n_dropped, " group(s) with no within-group outcome variation dropped ",
            "(all-0 or all-1) -- they carry no information for the conditional likelihood.")
  }

  coefs <- as.numeric(opt$par)
  names(coefs) <- colnames(X)
  vcov <- tryCatch(solve(opt$hessian), error = function(e) matrix(NA_real_, k, k))
  bread <- vcov
  dimnames(vcov) <- dimnames(bread) <- list(colnames(X), colnames(X))

  # df.residual isn't well-defined for conditional logistic regression (the
  # per-group fixed effects are never estimated, only conditioned out), so
  # it's left NA; see n_used_groups/n_dropped_groups instead.
  list(coefficients = coefs, vcov = vcov, bread = bread, fitted.values = NULL,
       loglik = -opt$value, dispersion = 1, df.residual = NA_real_,
       n_used_groups = info$n_used_groups, n_dropped_groups = n_dropped,
       iterations = opt$counts[[1]], convergence = opt$convergence)
}

fit_random <- function(X, y, family, group_start, group_size, R, maxit, tol) {
  if (family$family == "gaussian") return(fit_random_gaussian(X, y, group_start, group_size, maxit, tol))
  if (family$family == "poisson") {
    scale <- apply(X, 2, function(col) { s <- stats::sd(col); if (s > 0) s else 1 })
    Xs <- sweep(X, 2, scale, "/")

    # Better starting values than zeros/1: pooled Poisson MLE for beta (on
    # the scaled design, matching the units random_poisson_fit_cpp works
    # in), and a method-of-moments dispersion from the pooled residual
    # over/under-dispersion. Both are cheap and make the Newton path much
    # less likely to wander into a near-singular Hessian region on sparse
    # or skewed panels.
    pooled <- irls_fit_cpp(Xs, y, 1L, 1L, maxit, tol)
    beta_start <- as.numeric(pooled$coefficients)
    mu_pooled <- as.numeric(pooled$fitted.values)
    resid_var <- mean((y - mu_pooled)^2)
    mu_bar <- mean(mu_pooled)
    d_start <- if (resid_var > mu_bar && mu_bar > 0) mu_bar^2 / (resid_var - mu_bar) else 1.0
    if (!is.finite(d_start) || d_start <= 0) d_start <- 1.0

    res <- random_poisson_fit_cpp(Xs, y, group_start, group_size, maxit, tol, beta_start, d_start)
    if (!is.null(res$damping_used) && res$damping_used > 0) {
      warning("random-effects Poisson fit: the Hessian was ill-conditioned and ",
              "required Levenberg-Marquardt damping (max ridge = ",
              format(res$damping_used, digits = 3), "); standard errors may be ",
              "less reliable than usual -- check for sparse or highly skewed panels.",
              call. = FALSE)
    }
    coefs <- as.numeric(res$coefficients) / scale
    names(coefs) <- colnames(X)
    scale_full <- c(scale, 1) # dispersion parameter is unscaled
    vcov_full <- res$vcov_unscaled / outer(scale_full, scale_full)
    vcov <- vcov_full[seq_len(ncol(X)), seq_len(ncol(X)), drop = FALSE]
    dimnames(vcov) <- list(colnames(X), colnames(X))
    return(list(coefficients = coefs, vcov = vcov, fitted.values = NULL,
                loglik = res$loglik, dispersion_param = res$dispersion_param,
                damping_used = res$damping_used,
                df.residual = nrow(X) - ncol(X) - 1, iterations = res$iterations))
  }
  if (family$family == "negbin") {
    scale <- apply(X, 2, function(col) { s <- stats::sd(col); if (s > 0) s else 1 })
    Xs <- sweep(X, 2, scale, "/")
    res <- random_negbin_fit_cpp(Xs, y, group_start, group_size, maxit, tol)
    coefs <- as.numeric(res$coefficients) / scale
    names(coefs) <- colnames(X)
    scale_full <- c(scale, 1, 1) # a, b are unscaled
    vcov_full <- res$vcov_unscaled / outer(scale_full, scale_full)
    vcov <- vcov_full[seq_len(ncol(X)), seq_len(ncol(X)), drop = FALSE]
    dimnames(vcov) <- list(colnames(X), colnames(X))
    return(list(coefficients = coefs, vcov = vcov, fitted.values = NULL,
                loglik = res$loglik, negbin_a = res$a, negbin_b = res$b,
                df.residual = nrow(X) - ncol(X) - 2, iterations = res$iterations))
  }
  if (family$family == "binomial") return(fit_random_binomial(X, y, family, group_start, group_size, R, maxit, tol))
  stop("model = 'random' is not yet implemented for family '", family$family, "'", call. = FALSE)
}

fit_random_gaussian <- function(X, y, group_start, group_size, maxit, tol) {
  has_intercept <- "(Intercept)" %in% colnames(X)
  Xk <- if (has_intercept) X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE] else X
  k <- ncol(Xk); n <- nrow(X); G <- length(group_start)
  Ti <- group_size

  dm <- within_demean_cpp(Xk, y, group_start, group_size)
  fe <- irls_fit_cpp(dm$X, dm$y, 0L, 0L, maxit, tol)
  df_fe <- n - G - k
  resid_fe <- dm$y - dm$X %*% fe$coefficients
  sigma_v2 <- sum(resid_fe^2) / max(1, df_fe)

  gm <- group_means_cpp(X, y, group_start, group_size)
  bw <- irls_fit_cpp(gm$Xbar, gm$ybar, 0L, 0L, maxit, tol)
  df_bw <- G - ncol(X)
  resid_bw <- gm$ybar - gm$Xbar %*% bw$coefficients
  sigma_1_2 <- sum(resid_bw^2) / max(1, df_bw)

  sigma_mu2 <- max(0, sigma_1_2 - sigma_v2 / mean(Ti))
  theta <- 1 - sqrt(sigma_v2 / (Ti * sigma_mu2 + sigma_v2))

  qd <- quasi_demean_cpp(X, y, group_start, group_size, theta)
  res <- irls_fit_cpp(qd$X, qd$y, 0L, 0L, maxit, tol)

  coefs <- as.numeric(res$coefficients)
  names(coefs) <- colnames(X)
  df_resid <- n - ncol(X)
  resid_final <- qd$y - qd$X %*% res$coefficients
  sigma2 <- sum(resid_final^2) / max(1, df_resid)
  vcov <- res$vcov_unscaled * sigma2
  dimnames(vcov) <- list(colnames(X), colnames(X))

  list(coefficients = coefs, vcov = vcov, fitted.values = as.numeric(X %*% res$coefficients),
       loglik = NA_real_, dispersion = sigma2,
       sigma_v2 = sigma_v2, sigma_mu2 = sigma_mu2, theta = theta,
       df.residual = df_resid, iterations = res$iterations)
}

fit_random_binomial <- function(X, y, family, group_start, group_size, R, maxit, tol) {
  gh <- gauss_hermite_quadrature(R)
  k <- ncol(X)

  start <- irls_fit_cpp(X, y, 2L, family$link_id, maxit, tol)
  theta0 <- c(as.numeric(start$coefficients), log(0.5))

  # BFGS calls the objective and gradient separately at (usually) the same
  # theta; cache the last quadrature evaluation to avoid running it twice.
  cache <- new.env(parent = emptyenv())
  eval_at <- function(theta) {
    if (is.null(cache$theta) || !identical(cache$theta, theta)) {
      cache$theta <- theta
      beta <- theta[seq_len(k)]
      sigma <- exp(theta[k + 1])
      cache$res <- random_binomial_loglik_grad_cpp(beta, sigma, X, y, group_start, group_size,
                                                     gh$nodes, gh$weights, family$link_id)
      cache$sigma <- sigma
    }
    cache
  }
  objective <- function(theta) -eval_at(theta)$res$loglik
  gradient <- function(theta) {
    c <- eval_at(theta)
    g <- c$res$gradient
    g[k + 1] <- g[k + 1] * c$sigma # chain rule: d/d(log sigma) = d/dsigma * sigma
    -g
  }

  opt <- stats::optim(theta0, objective, gradient, method = "BFGS",
                       control = list(maxit = maxit, reltol = tol),
                       hessian = TRUE)

  coefs <- opt$par[seq_len(k)]
  names(coefs) <- colnames(X)
  sigma <- exp(opt$par[k + 1])
  vcov_full <- tryCatch(solve(opt$hessian), error = function(e) matrix(NA_real_, k + 1, k + 1))
  vcov <- vcov_full[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(vcov) <- list(colnames(X), colnames(X))

  list(coefficients = coefs, vcov = vcov, fitted.values = NULL,
       loglik = -opt$value, sigma = sigma, dispersion = 1,
       df.residual = nrow(X) - k - 1, iterations = opt$counts[[1]], convergence = opt$convergence)
}
