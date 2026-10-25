#' Random slopes: correlated random intercept and slope(s)
#'
#' The model is `g(E[y_it | b_i]) = x_it' beta + z_it' b_i` with
#' `b_i ~ N(0, Sigma)`, `z_it = (1, slope variables)` and `Sigma`
#' unstructured. The Gaussian model is fitted by exact maximum likelihood
#' with the residual variance profiled out (the deviance of
#' `lme4::lmer(REML = FALSE)`, in the relative Cholesky factor of `Sigma`).
#' The Poisson and binomial models (one slope, two-dimensional effect) use
#' adaptive Gauss-Hermite product quadrature centred at the mode of each
#' individual's integrand and scaled by its curvature; one node per
#' dimension is the Laplace approximation.
#'
#' @keywords internal
#' @noRd
NULL

#' Is a random-slope model requested? Validates the `random` formula.
#'
#' @keywords internal
#' @noRd
random_slopes_requested <- function(random) {
  if (!inherits(random, "formula") || length(random) != 2L) {
    stop("'random' must be a one-sided formula such as ~ 1 or ~ time", call. = FALSE)
  }
  tt <- stats::terms(random)
  if (attr(tt, "intercept") != 1L) {
    stop("'random' must keep the random intercept; slopes without an intercept ",
         "are not supported", call. = FALSE)
  }
  length(attr(tt, "term.labels")) > 0L
}

#' Random-effects design matrix Z (intercept first) for a data frame
#'
#' @keywords internal
#' @noRd
slopes_z <- function(random, data) {
  mf <- stats::model.frame(random, data = data, na.action = stats::na.pass)
  Z <- stats::model.matrix(random, data = mf)
  if (anyNA(Z)) {
    stop("the variables in 'random' must not contain missing values", call. = FALSE)
  }
  Z
}

#' Fit dispatcher for random slopes
#'
#' @keywords internal
#' @noRd
fit_random_slopes <- function(X, y, Z, family, group_size, R, maxit, tol) {
  group <- rep(seq_along(group_size), group_size)
  qr_z <- qr(Z, tol = 1e-7)
  if (qr_z$rank < ncol(Z)) {
    stop("the random-effects design (intercept and slopes in 'random') is rank deficient",
         call. = FALSE)
  }
  out <- if (family$family == "gaussian") {
    fit_slopes_gaussian(X, y, Z, group, maxit, tol)
  } else {
    fit_slopes_glm(X, y, Z, group, family, R, maxit, tol)
  }
  r <- ncol(Z)
  out$npar <- ncol(X) + r * (r + 1L) / 2L + as.integer(family$family == "gaussian")
  out
}

# ---------------------------------------------------------------- Gaussian

#' Per-individual cross-products of X, Z and y
#'
#' @keywords internal
#' @noRd
slopes_crossprods <- function(X, y, Z, group) {
  idx <- split(seq_along(group), group)
  lapply(idx, function(rows) {
    Xi <- X[rows, , drop = FALSE]; Zi <- Z[rows, , drop = FALSE]; yi <- y[rows]
    list(ZtZ = crossprod(Zi), ZtX = crossprod(Zi, Xi), Zty = crossprod(Zi, yi))
  })
}

slopes_lambda <- function(par, r) {
  L <- matrix(0, r, r)
  L[lower.tri(L, diag = TRUE)] <- par
  L
}

#' Profiled Gaussian quantities at a relative Cholesky factor
#'
#' `Sigma = sigma2 * Lambda Lambda'`; `M_i = I + Lambda' Z_i'Z_i Lambda`.
#'
#' @keywords internal
#' @noRd
slopes_gaussian_profile <- function(par, cp, XtX, Xty, yty, n, r) {
  L <- slopes_lambda(par, r)
  XVX <- XtX; XVy <- Xty; yVy <- yty; logdet <- 0
  Ir <- diag(r)
  for (cpi in cp) {
    M <- Ir + crossprod(L, cpi$ZtZ %*% L)
    Rm <- chol(M)
    logdet <- logdet + 2 * sum(log(diag(Rm)))
    A <- backsolve(Rm, crossprod(L, cpi$ZtX), transpose = TRUE)
    a <- backsolve(Rm, crossprod(L, cpi$Zty), transpose = TRUE)
    XVX <- XVX - crossprod(A)
    XVy <- XVy - crossprod(A, a)
    yVy <- yVy - sum(a^2)
  }
  beta <- solve(XVX, XVy)
  r2 <- as.numeric(yVy - crossprod(XVy, beta))
  sigma2 <- r2 / n
  list(dev = n * log(2 * pi * sigma2) + logdet + n, beta = as.numeric(beta),
       sigma2 = sigma2, XVX = XVX, L = L)
}

fit_slopes_gaussian <- function(X, y, Z, group, maxit, tol) {
  n <- length(y); r <- ncol(Z); k <- ncol(X)
  cp <- slopes_crossprods(X, y, Z, group)
  XtX <- crossprod(X); Xty <- crossprod(X, y); yty <- sum(y^2)
  lt <- lower.tri(diag(r), diag = TRUE)
  diag_pos <- which((row(lt) == col(lt))[lt])
  par0 <- diag(r)[lower.tri(diag(r), diag = TRUE)]
  lower <- rep(-Inf, length(par0)); lower[diag_pos] <- 0
  dev <- function(par) {
    v <- tryCatch(slopes_gaussian_profile(par, cp, XtX, Xty, yty, n, r)$dev,
                  error = function(e) NA_real_)
    if (is.finite(v)) v else 1e100
  }
  opt <- stats::nlminb(par0, dev, lower = lower,
                       control = list(iter.max = max(200L, maxit), eval.max = 2000L,
                                      rel.tol = 1e-14, x.tol = 1e-12))
  # polish with a bounded quasi-Newton step from the nlminb optimum
  opt2 <- stats::optim(opt$par, dev, method = "L-BFGS-B", lower = lower,
                       control = list(factr = 1, pgtol = 0, maxit = 500L))
  par <- if (opt2$value < opt$objective) opt2$par else opt$par
  pr <- slopes_gaussian_profile(par, cp, XtX, Xty, yty, n, r)

  coefs <- pr$beta; names(coefs) <- colnames(X)
  vcov <- pr$sigma2 * solve(pr$XVX)
  vcov <- (vcov + t(vcov)) / 2
  dimnames(vcov) <- list(colnames(X), colnames(X))
  varcomp <- pr$sigma2 * tcrossprod(pr$L)
  dimnames(varcomp) <- list(colnames(Z), colnames(Z))

  # BLUPs: b_i = Lambda M_i^{-1} Lambda' Z_i'(y_i - X_i beta)
  ranef <- t(vapply(cp, function(cpi) {
    M <- diag(r) + crossprod(pr$L, cpi$ZtZ %*% pr$L)
    as.numeric(pr$L %*% solve(M, crossprod(pr$L, cpi$Zty - cpi$ZtX %*% coefs)))
  }, numeric(r)))
  if (r == 1L) ranef <- matrix(ranef, ncol = 1L)
  colnames(ranef) <- colnames(Z)

  eta <- as.numeric(X %*% coefs)
  list(coefficients = coefs, vcov = vcov, fitted.values = eta,
       loglik = -pr$dev / 2, dispersion = pr$sigma2, sigma_v2 = pr$sigma2,
       varcomp = varcomp, ranef = ranef, slopes_lambda = pr$L,
       df.residual = n - k, iterations = opt$iterations,
       convergence = opt$convergence)
}

#' Cluster-robust covariance of a Gaussian random-slope fit
#'
#' Scores `X_i' V_i^{-1} e_i`, bread the model-based covariance, factor
#' G / (G - 1).
#'
#' @keywords internal
#' @noRd
slopes_gaussian_scores <- function(object) {
  group <- rep(seq_along(object$group_size), object$group_size)
  e <- object$y - as.numeric(object$X %*% object$coefficients)
  L <- object$slopes_lambda; r <- ncol(L)
  idx <- split(seq_along(group), group)
  t(vapply(idx, function(rows) {
    Xi <- object$X[rows, , drop = FALSE]; Zi <- object$Z[rows, , drop = FALSE]
    ZL <- Zi %*% L
    M <- diag(r) + crossprod(ZL)
    ei <- e[rows]
    as.numeric(crossprod(Xi, ei) - crossprod(Xi, ZL %*% solve(M, crossprod(ZL, ei)))) /
      object$sigma_v2
  }, numeric(ncol(object$X))))
}

robust_vcov_slopes <- function(object, cluster) {
  cl_obs <- panglm_cluster_vector(object, cluster)
  cl <- individual_clusters(cl_obs, object$group_size)
  k <- length(object$coefficients)
  if (object$family$family == "gaussian") {
    score <- slopes_gaussian_scores(object)
    if (k == 1L) score <- matrix(score, ncol = 1L)
    bread <- object$vcov_classical %||% object$vcov
  } else {
    score <- object$theta_score
    bread <- object$theta_vcov
  }
  v <- ml_sandwich(bread, score, rep(1, nrow(score)), cl, clusters_only = TRUE)
  v <- v[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(v) <- list(names(object$coefficients), names(object$coefficients))
  v
}

# ----------------------------------------------------- Poisson and binomial

#' Observation log-density and its first two eta-derivatives
#'
#' Returns `ll`, `d1` and `w = -d2` (nonnegative, the observed information).
#'
#' @keywords internal
#' @noRd
slopes_obs_terms <- function(family, y) {
  if (family$family == "poisson") {
    lfy <- lgamma(y + 1)
    return(function(eta, deriv = TRUE) {
      mu <- exp(eta)
      if (!deriv) return(y * eta - mu - lfy)
      list(ll = y * eta - mu - lfy, d1 = y - mu, w = mu)
    })
  }
  s <- 2 * y - 1
  if (family$link == "logit") {
    return(function(eta, deriv = TRUE) {
      ll <- stats::plogis(s * eta, log.p = TRUE)
      if (!deriv) return(ll)
      p <- stats::plogis(eta)
      list(ll = ll, d1 = y - p, w = p * (1 - p))
    })
  }
  function(eta, deriv = TRUE) {
    se <- s * eta
    ll <- stats::pnorm(se, log.p = TRUE)
    if (!deriv) return(ll)
    lam <- exp(stats::dnorm(se, log = TRUE) - ll)
    list(ll = ll, d1 = s * lam, w = lam * (lam + se))
  }
}

#' Sigma from the unconstrained parameters (log l11, l21, log l22)
#'
#' @keywords internal
#' @noRd
slopes_sigma2d <- function(par) {
  L <- matrix(c(exp(par[1]), par[2], 0, exp(par[3])), 2, 2)
  tcrossprod(L)
}

#' Individual integrated log-likelihoods of the two-dimensional GLMM
#'
#' `theta = (beta, log l11, l21, log l22)`, `Sigma = L L'`. The modes of
#' the last call are kept in `environment(contrib)$state$modes` (and warm
#' start the next call).
#'
#' @keywords internal
#' @noRd
slopes_glm_contrib <- function(family, X, y, Z, group, R) {
  k <- ncol(X); G <- max(group)
  z1 <- Z[, 1]; z2 <- Z[, 2]
  obs <- slopes_obs_terms(family, y)
  gh <- gauss_hermite_quadrature(R)
  grid <- expand.grid(a = seq_len(R), b = seq_len(R))
  zn1 <- gh$nodes[grid$a]; zn2 <- gh$nodes[grid$b]
  log_wk <- log(gh$weights[grid$a]) + log(gh$weights[grid$b]) + zn1^2 + zn2^2
  rs <- function(v) as.numeric(rowsum(v, group, reorder = FALSE))
  state <- new.env(parent = emptyenv())
  state$modes <- matrix(0, G, 2)
  function(theta, score) {
    beta <- theta[seq_len(k)]
    S <- slopes_sigma2d(theta[k + 1:3])
    detS <- S[1, 1] * S[2, 2] - S[1, 2]^2
    if (!is.finite(detS) || detS <= 0) return(list(ll = rep(-Inf, G)))
    Si <- matrix(c(S[2, 2], -S[1, 2], -S[1, 2], S[1, 1]), 2, 2) / detS
    eta0 <- as.numeric(X %*% beta)
    log_prior_const <- -log(2 * pi) - 0.5 * log(detS)
    hfun <- function(bb1, bb2) {
      rs(obs(eta0 + z1 * bb1[group] + z2 * bb2[group], deriv = FALSE)) + log_prior_const -
        0.5 * (Si[1, 1] * bb1^2 + 2 * Si[1, 2] * bb1 * bb2 + Si[2, 2] * bb2^2)
    }
    # Newton ascent on each individual's (concave) integrand, with
    # per-individual step halving
    newton <- function(b) {
      hb <- hfun(b[, 1], b[, 2])
      if (!all(is.finite(hb))) return(NULL)
      for (step in seq_len(500L)) {
        o <- obs(eta0 + z1 * b[group, 1] + z2 * b[group, 2])
        g1 <- rs(o$d1 * z1) - (Si[1, 1] * b[, 1] + Si[1, 2] * b[, 2])
        g2 <- rs(o$d1 * z2) - (Si[2, 1] * b[, 1] + Si[2, 2] * b[, 2])
        h11 <- rs(o$w * z1 * z1) + Si[1, 1]
        h12 <- rs(o$w * z1 * z2) + Si[1, 2]
        h22 <- rs(o$w * z2 * z2) + Si[2, 2]
        dh <- h11 * h22 - h12^2
        d1 <- (h22 * g1 - h12 * g2) / dh
        d2 <- (h11 * g2 - h12 * g1) / dh
        if (!all(is.finite(d1)) || !all(is.finite(d2))) return(NULL)
        if (max(abs(d1), abs(d2)) < 1e-10) {
          b[, 1] <- b[, 1] + d1; b[, 2] <- b[, 2] + d2
          return(b)
        }
        big <- pmax(1, abs(d1) / 10, abs(d2) / 10)
        d1 <- d1 / big; d2 <- d2 / big
        if (max(abs(d1), abs(d2)) < 1e-3) {
          # quadratic convergence region: full Newton steps
          b[, 1] <- b[, 1] + d1; b[, 2] <- b[, 2] + d2
          hb <- hfun(b[, 1], b[, 2])
          next
        }
        for (halving in 0:40) {
          hn <- hfun(b[, 1] + d1, b[, 2] + d2)
          bad <- !is.finite(hn) | hn < hb - 1e-12 * pmax(1, abs(hb))
          if (!any(bad)) break
          d1[bad] <- d1[bad] / 2; d2[bad] <- d2[bad] / 2
        }
        if (all(bad)) return(b)
        d1[bad] <- 0; d2[bad] <- 0
        b[, 1] <- b[, 1] + d1; b[, 2] <- b[, 2] + d2
        hb[!bad] <- hn[!bad]
      }
      NULL
    }
    b <- newton(state$modes)
    if (is.null(b)) b <- newton(matrix(0, G, 2))
    if (is.null(b)) return(list(ll = rep(-Inf, G)))
    state$modes <- b
    # curvature at the mode
    o <- obs(eta0 + z1 * b[group, 1] + z2 * b[group, 2])
    h11 <- rs(o$w * z1 * z1) + Si[1, 1]
    h12 <- rs(o$w * z1 * z2) + Si[1, 2]
    h22 <- rs(o$w * z2 * z2) + Si[2, 2]
    dh <- h11 * h22 - h12^2
    # C C' = H^{-1}
    c11 <- sqrt(h22 / dh)
    c21 <- (-h12 / dh) / c11
    c22 <- sqrt(pmax(h11 / dh - c21^2, 0))
    terms <- vapply(seq_along(zn1), function(q) {
      bb1 <- b[, 1] + sqrt(2) * c11 * zn1[q]
      bb2 <- b[, 2] + sqrt(2) * (c21 * zn1[q] + c22 * zn2[q])
      log_wk[q] + hfun(bb1, bb2)
    }, numeric(G))
    if (!is.matrix(terms)) terms <- matrix(terms, nrow = G)
    list(ll = log(2 * pi) + log(c11 * c22) + log_row_sum_exp(terms))
  }
}

fit_slopes_glm <- function(X, y, Z, group, family, R, maxit, tol) {
  if (!family$family %in% c("poisson", "binomial")) {
    stop("random slopes are available for family = 'gaussian', 'poisson' and 'binomial' ",
         "(logit or probit) only", call. = FALSE)
  }
  if (ncol(Z) != 2L) {
    stop("random slopes for family = '", family$family, "' support one slope ",
         "(random = ~ x, a two-dimensional random effect); use family = 'gaussian' ",
         "for more slopes", call. = FALSE)
  }
  if (family$family == "binomial" && any(!y %in% c(0, 1))) {
    stop("random slopes for family = 'binomial' need a 0/1 response", call. = FALSE)
  }
  k <- ncol(X); G <- max(group)
  fam_id <- if (family$family == "poisson") 1L else 2L
  link_id <- if (fam_id == 2L) family$link_id else 1L
  b0 <- as.numeric(irls_fit_cpp(X, y, fam_id, link_id, maxit, tol)$coefficients)
  sz <- stats::sd(Z[, 2]); if (!is.finite(sz) || sz <= 0) sz <- 1
  start <- c(b0, log(0.5), 0, log(0.2 / sz))
  contrib <- slopes_glm_contrib(family, X, y, Z, group, R)
  res <- ml_optimize(start, contrib, rep(1, G), maxit, tol)
  # recover the modes at the optimum
  contrib(res$theta, FALSE)
  modes <- environment(contrib)$state$modes
  colnames(modes) <- colnames(Z)

  coefs <- res$theta[seq_len(k)]; names(coefs) <- colnames(X)
  vcov <- res$vcov[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(vcov) <- list(colnames(X), colnames(X))
  varcomp <- slopes_sigma2d(res$theta[k + 1:3])
  dimnames(varcomp) <- list(colnames(Z), colnames(Z))
  eta <- as.numeric(X %*% coefs)
  out <- list(coefficients = coefs, vcov = vcov, loglik = res$loglik,
              theta = res$theta, theta_vcov = res$vcov, theta_score = res$score,
              varcomp = varcomp, ranef = modes, dispersion = 1,
              quadrature_nodes = R, df.residual = length(y) - length(res$theta),
              iterations = res$iterations, convergence = res$convergence)
  out$fitted.values <- slopes_marginal_mean(family, eta, Z, varcomp)
  out
}

#' Population-averaged mean: E over b of g^{-1}(eta + z'b)
#'
#' @keywords internal
#' @noRd
slopes_marginal_mean <- function(family, eta, Z, varcomp, nodes = 21L) {
  s2 <- as.numeric(rowSums((Z %*% varcomp) * Z))
  switch(family$family,
    gaussian = eta,
    poisson = exp(eta + s2 / 2),
    binomial = {
      gh <- gauss_hermite_quadrature(nodes)
      shifts <- sqrt(2) * outer(sqrt(s2), gh$nodes)
      as.numeric(linkinv_r(eta + shifts, family$link_id) %*% gh$weights)
    }
  )
}
