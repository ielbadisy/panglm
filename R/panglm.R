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
#' @param R number of Gauss-Hermite quadrature nodes (random-effects
#'   binomial models only)
#' @param maxit maximum IRLS/Newton iterations
#' @param tol convergence tolerance
#' @return an object of class `"panglm"`
#' @export
panglm <- function(formula, data, index,
                    model = c("pooling", "within", "random"),
                    family = "gaussian",
                    R = 21, maxit = 100, tol = 1e-10) {
  model <- match.arg(model)
  family <- resolve_family(family)

  if (missing(index)) stop("'index' is required, e.g. index = c(\"id\", \"time\")", call. = FALSE)

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

  has_intercept <- "(Intercept)" %in% colnames(X_full)
  X_noint <- if (has_intercept) X_full[, setdiff(colnames(X_full), "(Intercept)"), drop = FALSE] else X_full

  fit <- switch(model,
    pooling = fit_pooled(X_full, y, family, maxit, tol),
    within  = fit_within(X_noint, y, family, group_start, group_size, maxit, tol),
    random  = fit_random(X_full, y, family, group_start, group_size, R, maxit, tol)
  )

  fit$y <- y
  fit$call <- match.call()
  fit$formula <- formula
  fit$model <- model
  fit$family <- family
  fit$index <- index
  fit$nobs <- length(y)
  fit$n_groups <- panel$n_groups
  class(fit) <- "panglm"
  fit
}

fit_pooled <- function(X, y, family, maxit, tol) {
  res <- irls_fit_cpp(X, y, family$family_id, family$link_id, maxit, tol)
  vcov <- res$vcov_unscaled * res$dispersion
  coefs <- as.numeric(res$coefficients)
  names(coefs) <- colnames(X)
  dimnames(vcov) <- list(colnames(X), colnames(X))
  list(coefficients = coefs, vcov = vcov, fitted.values = as.numeric(res$fitted.values),
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
    coefs <- as.numeric(res$coefficients)
    names(coefs) <- colnames(X)
    dimnames(vcov) <- list(colnames(X), colnames(X))
    gm <- group_means_cpp(X, y, group_start, group_size)
    alpha_i <- gm$ybar - as.numeric(gm$Xbar %*% res$coefficients)
    alpha_obs <- rep(alpha_i, group_size)
    fitted <- alpha_obs + as.numeric(X %*% res$coefficients)
    return(list(coefficients = coefs, vcov = vcov, fitted.values = fitted,
                loglik = NA_real_, dispersion = sigma2, df.residual = df_resid,
                iterations = res$iterations))
  }
  if (family$family %in% c("poisson", "negbin")) {
    scale <- apply(X, 2, function(col) { s <- stats::sd(col); if (s > 0) s else 1 })
    Xs <- sweep(X, 2, scale, "/")
    res <- within_poisson_fit_cpp(Xs, y, group_start, group_size, maxit, tol)
    coefs <- as.numeric(res$coefficients) / scale
    names(coefs) <- colnames(X)
    vcov <- res$vcov_unscaled / outer(scale, scale)
    dimnames(vcov) <- list(colnames(X), colnames(X))
    return(list(coefficients = coefs, vcov = vcov, fitted.values = NULL,
                loglik = res$loglik, dispersion = 1, df.residual = nrow(X) - ncol(X) - length(group_start),
                iterations = res$iterations))
  }
  stop("model = 'within' is not yet implemented for family '", family$family,
       "' (conditional/fixed-effects binomial requires a separate estimator)", call. = FALSE)
}

fit_random <- function(X, y, family, group_start, group_size, R, maxit, tol) {
  if (family$family == "gaussian") return(fit_random_gaussian(X, y, group_start, group_size, maxit, tol))
  if (family$family %in% c("poisson", "negbin")) {
    scale <- apply(X, 2, function(col) { s <- stats::sd(col); if (s > 0) s else 1 })
    Xs <- sweep(X, 2, scale, "/")
    res <- random_poisson_fit_cpp(Xs, y, group_start, group_size, maxit, tol)
    coefs <- as.numeric(res$coefficients) / scale
    names(coefs) <- colnames(X)
    scale_full <- c(scale, 1) # dispersion parameter is unscaled
    vcov_full <- res$vcov_unscaled / outer(scale_full, scale_full)
    vcov <- vcov_full[seq_len(ncol(X)), seq_len(ncol(X)), drop = FALSE]
    dimnames(vcov) <- list(colnames(X), colnames(X))
    return(list(coefficients = coefs, vcov = vcov, fitted.values = NULL,
                loglik = res$loglik, dispersion_param = res$dispersion_param,
                df.residual = nrow(X) - ncol(X) - 1, iterations = res$iterations))
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

  objective <- function(theta) {
    beta <- theta[seq_len(k)]
    sigma <- exp(theta[k + 1])
    res <- random_binomial_loglik_grad_cpp(beta, sigma, X, y, group_start, group_size,
                                            gh$nodes, gh$weights, family$link_id)
    -res$loglik
  }
  gradient <- function(theta) {
    beta <- theta[seq_len(k)]
    sigma <- exp(theta[k + 1])
    res <- random_binomial_loglik_grad_cpp(beta, sigma, X, y, group_start, group_size,
                                            gh$nodes, gh$weights, family$link_id)
    g <- res$gradient
    g[k + 1] <- g[k + 1] * sigma # chain rule: d/d(log sigma) = d/dsigma * sigma
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
