#' Random slopes: correlated random intercept and slope(s)
#'
#' The model is `g(E[y_it | b_i]) = x_it' beta + z_it' b_i` with
#' `b_i ~ N(0, Sigma)`, `z_it = (1, slope variables)` and `Sigma`
#' unstructured. The Gaussian model is fitted by exact maximum likelihood
#' with the residual variance profiled out (the deviance of
#' `lme4::lmer(REML = FALSE)`, in the relative Cholesky factor of `Sigma`).
#' The Poisson and binomial models (any number of slopes) use adaptive
#' Gauss-Hermite product quadrature centred at the mode of each
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
  out$npar <- ncol(X) + r * (r + 1L) / 2L +
    as.integer(family$family %in% c("gaussian", "negbin", "gamma"))
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

#' Observation log-density and its eta-derivatives
#'
#' Returns a function of `eta` (a vector or a matrix of columns of `eta`)
#' and of `alpha`, the log of the dispersion parameter of the negative
#' binomial (NB2 `theta`) and gamma (`shape`) families. With `deriv = 0`
#' it returns the log density `ll`; with `deriv = 1` also `d1` (its
#' derivative in eta), `w = -d2` (nonnegative, the observed information)
#' and `w1 = dw / deta`; with `deriv = 2` also the derivatives of `ll`,
#' `d1` and `w` in `alpha` (`da_ll`, `da_d1`, `da_w`).
#'
#' @keywords internal
#' @noRd
slopes_obs_terms <- function(family, y) {
  if (family$family == "poisson") {
    lfy <- lgamma(y + 1)
    return(function(eta, deriv = 1L, alpha = NULL) {
      mu <- exp(eta)
      ll <- y * eta - mu - lfy
      if (!deriv) return(ll)
      list(ll = ll, d1 = y - mu, w = mu, w1 = mu)
    })
  }
  if (family$family == "negbin") {
    lfy <- lgamma(y + 1)
    return(function(eta, deriv = 1L, alpha) {
      th <- exp(alpha)
      mu <- exp(eta)
      tm <- th + mu
      ll <- lgamma(y + th) - lgamma(th) - lfy + th * alpha + y * eta - (y + th) * log(tm)
      if (!deriv) return(ll)
      w <- th * mu * (th + y) / tm^2
      out <- list(ll = ll, d1 = th * (y - mu) / tm, w = w, w1 = w * (th - mu) / tm)
      if (deriv >= 2L) {
        out$da_ll <- th * (digamma(y + th) - digamma(th) + alpha + 1 - log(tm) - (y + th) / tm)
        out$da_d1 <- th * (y - mu) * mu / tm^2
        out$da_w <- th * mu * (2 * th * mu + y * mu - th * y) / tm^3
      }
      out
    })
  }
  if (family$family == "gamma") {
    ly <- log(y)
    return(function(eta, deriv = 1L, alpha) {
      nu <- exp(alpha)
      r <- y * exp(-eta)
      ll <- nu * alpha - lgamma(nu) + (nu - 1) * ly - nu * eta - nu * r
      if (!deriv) return(ll)
      out <- list(ll = ll, d1 = nu * (r - 1), w = nu * r, w1 = -nu * r)
      if (deriv >= 2L) {
        out$da_ll <- nu * (alpha + 1 - digamma(nu) + ly - eta - r)
        out$da_d1 <- out$d1
        out$da_w <- out$w
      }
      out
    })
  }
  s <- 2 * y - 1
  if (family$link == "logit") {
    return(function(eta, deriv = 1L, alpha = NULL) {
      ll <- stats::plogis(s * eta, log.p = TRUE)
      if (!deriv) return(ll)
      p <- stats::plogis(eta)
      w <- p * (1 - p)
      list(ll = ll, d1 = y - p, w = w, w1 = w * (1 - 2 * p))
    })
  }
  function(eta, deriv = 1L, alpha = NULL) {
    se <- s * eta
    ll <- stats::pnorm(se, log.p = TRUE)
    if (!deriv) return(ll)
    lam <- exp(stats::dnorm(se, log = TRUE) - ll)
    w <- lam * (lam + se)
    list(ll = ll, d1 = s * lam, w = w, w1 = s * (lam - w * (2 * lam + se)))
  }
}

# ------------------------------------- batched small-matrix linear algebra
# A batch of G symmetric q x q matrices is an array of dimension c(G, q, q);
# A[, a, c] is the vector of the (a, c) entries over the G individuals.

#' Lower Cholesky factors of a batch of positive definite matrices
#'
#' @keywords internal
#' @noRd
bchol <- function(A) {
  q <- dim(A)[2]
  L <- array(0, dim(A))
  for (c in seq_len(q)) {
    s <- A[, c, c]
    for (m in seq_len(c - 1L)) s <- s - L[, c, m]^2
    L[, c, c] <- sqrt(ifelse(s >= 0, s, NaN))
    for (a in seq_len(q - c) + c) {
      t <- A[, a, c]
      for (m in seq_len(c - 1L)) t <- t - L[, a, m] * L[, c, m]
      L[, a, c] <- t / L[, c, c]
    }
  }
  L
}

#' Inverses of a batch of lower triangular matrices
#'
#' @keywords internal
#' @noRd
blowinv <- function(L) {
  q <- dim(L)[2]
  M <- array(0, dim(L))
  for (c in seq_len(q)) {
    M[, c, c] <- 1 / L[, c, c]
    for (a in seq_len(q - c) + c) {
      s <- 0
      for (m in c:(a - 1L)) s <- s + L[, a, m] * M[, m, c]
      M[, a, c] <- -s / L[, a, a]
    }
  }
  M
}

#' Inverses of a batch of positive definite matrices
#'
#' @keywords internal
#' @noRd
bspdinv <- function(A) {
  q <- dim(A)[2]
  M <- blowinv(bchol(A))
  P <- array(0, dim(A))
  for (a in seq_len(q)) for (c in seq_len(a)) {
    s <- 0
    for (m in a:q) s <- s + M[, m, a] * M[, m, c]
    P[, a, c] <- s; P[, c, a] <- s
  }
  P
}

#' Batch matrix product: slice i of the result is A_i B_i (tA: A_i' B_i)
#'
#' @keywords internal
#' @noRd
bmatmul <- function(A, B, tA = FALSE) {
  q <- dim(A)[2]
  out <- array(0, dim(A))
  for (a in seq_len(q)) for (c in seq_len(q)) for (m in seq_len(q)) {
    out[, a, c] <- out[, a, c] + (if (tA) A[, m, a] else A[, a, m]) * B[, m, c]
  }
  out
}

#' Batch matrix-vector product: row i of the result is A_i v_i
#'
#' @keywords internal
#' @noRd
bmatvec <- function(A, v) {
  q <- dim(A)[2]
  out <- matrix(0, nrow(v), q)
  for (a in seq_len(q)) for (c in seq_len(q)) out[, a] <- out[, a] + A[, a, c] * v[, c]
  out
}

# ------------------------------------------- GLM families with random slopes

#' Positions of the vech of a q x q lower triangle (column-major)
#'
#' @keywords internal
#' @noRd
slopes_vech_index <- function(q) {
  lt <- lower.tri(diag(q), diag = TRUE)
  data.frame(row = row(lt)[lt], col = col(lt)[lt])
}

#' Cholesky factor L of Sigma from its parameters
#'
#' The parameters are the vech of L, unconstrained: changing the sign of a
#' column of L leaves Sigma unchanged, and a zero diagonal entry is a
#' singular Sigma (a variance on the boundary). Unlike a log-scale
#' diagonal, this parameterisation stays regular at the boundary (the
#' likelihood is even in a zero diagonal entry), so the information matrix
#' and the fixed-effect standard errors remain finite there, as in lme4.
#'
#' @keywords internal
#' @noRd
slopes_lfactor <- function(par, q) {
  vi <- slopes_vech_index(q)
  L <- matrix(0, q, q)
  L[cbind(vi$row, vi$col)] <- par
  L
}

#' Individual integrated log-likelihoods of a GLMM with q random effects
#'
#' `b_i = L u_i` with `u_i ~ N(0, I_q)`, `Sigma = L L'` and
#' `theta = (beta, vech L)`. The integral over `u_i` is computed by
#' adaptive Gauss-Hermite product quadrature: `R^q` nodes centred at the
#' mode of each individual's integrand and scaled by the Cholesky factor
#' of the inverse curvature there (one node per dimension is the Laplace
#' approximation). The rule is affine invariant, so it gives the same
#' value as adaptive quadrature over `b_i` directly. The modes of the last
#' call are kept in `environment(contrib)$state$modes` (on the `u` scale)
#' and warm start the next call.
#'
#' @keywords internal
#' @noRd
slopes_glm_contrib <- function(family, X, y, Z, group, R) {
  k <- ncol(X); G <- max(group); q <- ncol(Z)
  obs_alpha <- slopes_obs_terms(family, y)
  has_alpha <- family$family %in% c("negbin", "gamma")
  gh <- gauss_hermite_quadrature(R)
  grid <- as.matrix(expand.grid(rep(list(seq_len(R)), q)))
  zn <- matrix(gh$nodes[grid], ncol = q)
  log_wk <- rowSums(matrix(log(gh$weights[grid]), ncol = q)) + rowSums(zn^2)
  K <- nrow(zn)
  vi <- slopes_vech_index(q)
  pairs <- vi[vi$row >= vi$col, ]
  half_log_2pi <- 0.5 * q * log(2 * pi)
  rs <- function(v) rowsum(v, group, reorder = FALSE)
  state <- new.env(parent = emptyenv())
  state$modes <- matrix(0, G, q)
  # an infeasible point: zero scores spare the numerical derivatives
  fail <- function(p) list(ll = rep(-Inf, G), score = matrix(0, G, p))
  function(theta, score) {
    beta <- theta[seq_len(k)]
    L <- slopes_lfactor(theta[k + seq_len(nrow(vi))], q)
    if (any(!is.finite(L))) return(fail(length(theta)))
    alpha <- if (has_alpha) theta[length(theta)]
    if (has_alpha && (!is.finite(alpha) || abs(alpha) > 30)) return(fail(length(theta)))
    obs <- function(eta, deriv = 1L) obs_alpha(eta, deriv, alpha)
    Zt <- Z %*% L
    ZZ <- Zt[, pairs$row, drop = FALSE] * Zt[, pairs$col, drop = FALSE]
    eta0 <- as.numeric(X %*% beta)
    lin <- function(u) eta0 + rowSums(Zt * u[group, , drop = FALSE])
    hfun <- function(u) {
      as.numeric(rs(obs(lin(u), deriv = 0L))) - 0.5 * rowSums(u^2) - half_log_2pi
    }
    curvature <- function(w) {
      S <- rs(w * ZZ)
      H <- array(0, c(G, q, q))
      for (j in seq_len(nrow(pairs))) {
        H[, pairs$row[j], pairs$col[j]] <- S[, j]
        H[, pairs$col[j], pairs$row[j]] <- S[, j]
      }
      for (a in seq_len(q)) H[, a, a] <- H[, a, a] + 1
      H
    }
    # Newton ascent on each individual's (concave) integrand, with
    # per-individual step halving
    newton <- function(u) {
      hb <- hfun(u)
      if (!all(is.finite(hb))) return(NULL)
      for (step in seq_len(500L)) {
        o <- obs(lin(u))
        g <- rs(o$d1 * Zt) - u
        d <- bmatvec(bspdinv(curvature(o$w)), g)
        if (!all(is.finite(d))) return(NULL)
        if (max(abs(d)) < 1e-10) return(u + d)
        d <- d / pmax(1, apply(abs(d), 1L, max) / 10)
        if (max(abs(d)) < 1e-3) {
          # quadratic convergence region: full Newton steps
          u <- u + d
          hb <- hfun(u)
          next
        }
        for (halving in 0:40) {
          hn <- hfun(u + d)
          bad <- !is.finite(hn) | hn < hb - 1e-12 * pmax(1, abs(hb))
          if (!any(bad)) break
          d[bad, ] <- d[bad, ] / 2
        }
        if (all(bad)) return(u)
        d[bad, ] <- 0
        u <- u + d
        hb[!bad] <- hn[!bad]
      }
      NULL
    }
    u <- newton(state$modes)
    if (is.null(u)) u <- newton(matrix(0, G, q))
    if (is.null(u)) return(fail(length(theta)))
    state$modes <- u
    # C C' = H^{-1} at the mode; nodes u + sqrt(2) C z_k
    o <- obs(lin(u), if (score && has_alpha) 2L else 1L)
    P <- bspdinv(curvature(o$w))
    C <- bchol(P)
    log_det_C <- 0
    for (a in seq_len(q)) log_det_C <- log_det_C + log(C[, a, a])
    U <- lapply(seq_len(q), function(a) {
      Ua <- matrix(u[, a], G, K)
      for (c in seq_len(a)) Ua <- Ua + sqrt(2) * outer(C[, a, c], zn[, c])
      Ua
    })
    eta <- matrix(eta0, length(eta0), K)
    for (a in seq_len(q)) eta <- eta + Zt[, a] * U[[a]][group, , drop = FALSE]
    on <- obs(eta, if (!score) 0L else if (has_alpha) 2L else 1L)
    terms <- rs(if (score) on$ll else on) - half_log_2pi +
      matrix(log_wk, G, K, byrow = TRUE)
    for (a in seq_len(q)) terms <- terms - 0.5 * U[[a]]^2
    lse <- log_row_sum_exp(terms)
    ll <- half_log_2pi + log_det_C + lse
    if (!score) return(list(ll = ll))
    list(ll = ll, score = slopes_glm_score(X, Z, Zt, group, u, P, C, U, zn,
                                           exp(terms - lse), on, o, vi, has_alpha))
  }
}

#' Analytic scores of the adaptive quadrature log-likelihood
#'
#' Total derivatives of `log|C| + log sum_k w_k exp(h(u_k))` with respect to
#' `theta = (beta, vech L)`, where the nodes `u_k = u + sqrt(2) C z_k` move
#' with the mode `u` (implicit function theorem: `du = H^{-1} g`) and with
#' the Cholesky factor `C` of `H^{-1}` (`dC = -C Phi(C' dH C)`, `Phi` the
#' lower triangle with half the diagonal). `pk` holds the normalised node
#' weights of each individual, `on` the observation terms at the nodes and
#' `o` those at the mode. With `has_alpha` the last parameter is the log
#' dispersion parameter, which enters the log density only.
#'
#' @keywords internal
#' @noRd
slopes_glm_score <- function(X, Z, Zt, group, u, P, C, U, zn, pk, on, o, vi,
                             has_alpha = FALSE) {
  k <- ncol(X); q <- ncol(Z); G <- nrow(u)
  rs <- function(v) rowsum(v, group, reorder = FALSE)
  E1 <- pk[group, , drop = FALSE] * on$d1
  e1 <- rowSums(E1)
  EU <- vapply(U, function(Ua) rowSums(pk * Ua), numeric(G))
  if (G == 1L) EU <- matrix(EU, 1L)
  # expected gradient of h over the nodes, and M[m, c] = E[grad_m h z_c]
  Egrad <- rs(e1 * Zt) - EU
  Ez <- E1 %*% zn
  M <- array(0, c(G, q, q))
  for (m in seq_len(q)) {
    PU <- (pk * U[[m]]) %*% zn
    for (c in seq_len(q)) M[, m, c] <- rs(Zt[, m] * Ez[, c]) - PU[, c]
  }
  W <- o$w; W1 <- o$w1
  p <- k + nrow(vi) + has_alpha
  out <- matrix(0, G, p)
  ZZ <- Zt[, rep(seq_len(q), q), drop = FALSE] * Zt[, rep(seq_len(q), each = q), drop = FALSE]
  # derivatives of the mode (du), of the curvature (dH) and the expected
  # direct derivative of h for parameter j
  pieces <- function(j) {
    if (j > k + nrow(vi)) {
      # the log dispersion: no direct effect on eta
      du <- bmatvec(P, rs(o$da_d1 * Zt))
      deta <- rowSums(Zt * du[group, , drop = FALSE])
      return(list(du = du, dH = array(rs((W1 * deta + o$da_w) * ZZ), c(G, q, q)),
                  Edirect = rs(rowSums(pk[group, , drop = FALSE] * on$da_ll))))
    }
    if (j <= k) {
      direct <- X[, j]
      Edirect <- rs(e1 * direct)
    } else {
      a <- vi$row[j - k]; cc <- vi$col[j - k]
      direct <- Z[, a] * u[group, cc]
      Edirect <- rs(Z[, a] * rowSums(E1 * U[[cc]][group, , drop = FALSE]))
    }
    g <- -rs(W * direct * Zt)
    if (j > k) g[, cc] <- g[, cc] + rs(o$d1 * Z[, a])
    du <- bmatvec(P, g)
    deta <- direct + rowSums(Zt * du[group, , drop = FALSE])
    dH <- array(rs((W1 * deta) * ZZ), c(G, q, q))
    if (j > k) {
      v <- rs(W * Z[, a] * Zt)
      for (m in seq_len(q)) {
        dH[, m, cc] <- dH[, m, cc] + v[, m]
        dH[, cc, m] <- dH[, cc, m] + v[, m]
      }
    }
    list(du = du, dH = dH, Edirect = Edirect)
  }
  for (j in seq_len(p)) {
    pj <- pieces(j)
    du <- pj$du; dH <- pj$dH
    dlogC <- 0
    for (m in seq_len(q)) for (c in seq_len(q)) dlogC <- dlogC - 0.5 * P[, m, c] * dH[, c, m]
    Phi <- bmatmul(C, bmatmul(dH, C), tA = TRUE)
    for (m in seq_len(q)) {
      Phi[, m, m] <- Phi[, m, m] / 2
      for (c in seq_len(q)) if (c > m) Phi[, m, c] <- 0
    }
    dC <- -bmatmul(C, Phi)
    out[, j] <- dlogC + pj$Edirect + rowSums(Egrad * du) + sqrt(2) * rowSums(M * dC, dims = 1L)
  }
  out
}

#' Default number of quadrature nodes per dimension for q random effects
#'
#' @keywords internal
#' @noRd
slopes_default_nodes <- function(q) {
  if (q <= 2L) 7L else if (q == 3L) 5L else if (q <= 6L) 3L else 1L
}

#' Number of nodes per dimension: the default for q effects, or R with a
#' warning when the product rule is large
#'
#' @keywords internal
#' @noRd
slopes_nodes <- function(R, q) {
  if (is.null(R)) return(slopes_default_nodes(q))
  if (R^q > 2000) {
    warning("adaptive quadrature with R = ", R, " nodes per dimension and ", q,
            " random effects uses ", R^q, " nodes per individual and may be slow; ",
            "consider a smaller R (R = 1 is the Laplace approximation)", call. = FALSE)
  }
  R
}

fit_slopes_glm <- function(X, y, Z, group, family, R, maxit, tol) {
  if (!family$family %in% c("poisson", "binomial", "negbin", "gamma")) {
    stop("random slopes are available for family = 'gaussian', 'poisson', 'negbin', ",
         "'gamma' and 'binomial' (logit or probit) only", call. = FALSE)
  }
  if (family$family == "binomial" && any(!y %in% c(0, 1))) {
    stop("random slopes for family = 'binomial' need a 0/1 response", call. = FALSE)
  }
  k <- ncol(X); G <- max(group); q <- ncol(Z)
  R <- slopes_nodes(R, q)
  alpha0 <- NULL
  if (family$family == "negbin") {
    pooled <- negbin_irls_fit_cpp(X, y, maxit, tol, NULL)
    b0 <- as.numeric(pooled$coefficients)
    alpha0 <- log(pooled$theta)
  } else {
    fam_id <- switch(family$family, poisson = 1L, binomial = 2L, gamma = 3L)
    link_id <- if (fam_id == 2L) family$link_id else 1L
    pooled <- irls_fit_cpp(X, y, fam_id, link_id, maxit, tol)
    b0 <- as.numeric(pooled$coefficients)
    if (family$family == "gamma") alpha0 <- -log(max(1e-3, pooled$dispersion))
  }
  vi <- slopes_vech_index(q)
  sz <- vapply(seq_len(q), function(a) {
    s <- stats::sd(Z[, a]); if (is.finite(s) && s > 0) s else 1
  }, numeric(1))
  l0 <- ifelse(vi$row == vi$col, ifelse(vi$row == 1L, 0.5, 0.2 / sz[vi$row]), 0)
  start <- c(b0, l0, alpha0)
  contrib <- slopes_glm_contrib(family, X, y, Z, group, R)
  res <- ml_optimize(start, contrib, rep(1, G), maxit, tol)
  # recover the modes at the optimum
  contrib(res$theta, FALSE)
  L <- slopes_lfactor(res$theta[k + seq_len(nrow(vi))], q)
  modes <- environment(contrib)$state$modes %*% t(L)
  colnames(modes) <- colnames(Z)

  slopes_boundary_warning(L, sz, colnames(Z))
  theta_vcov <- slopes_theta_vcov(res$hessian)

  coefs <- res$theta[seq_len(k)]; names(coefs) <- colnames(X)
  vcov <- theta_vcov[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(vcov) <- list(colnames(X), colnames(X))
  varcomp <- tcrossprod(L)
  dimnames(varcomp) <- list(colnames(Z), colnames(Z))
  eta <- as.numeric(X %*% coefs)
  dispersion <- 1
  extra <- list()
  if (!is.null(alpha0)) {
    alpha <- exp(res$theta[length(res$theta)])
    if (family$family == "negbin") {
      extra$negbin_theta <- alpha
    } else {
      extra$shape <- alpha
      dispersion <- 1 / alpha
    }
  }
  out <- list(coefficients = coefs, vcov = vcov, loglik = res$loglik,
              theta = res$theta, theta_vcov = theta_vcov, theta_score = res$score,
              varcomp = varcomp, ranef = modes, dispersion = dispersion,
              quadrature_nodes = R, df.residual = length(y) - length(res$theta),
              iterations = res$iterations, convergence = res$convergence)
  out$fitted.values <- slopes_marginal_mean(family, eta, Z, varcomp)
  c(out, extra)
}

#' Warn when the random-effects covariance is singular
#'
#' `|L_aa|` is the standard deviation of effect `a` given the preceding
#' ones; it is compared with 1e-4 on the scale of the slope variable (the
#' singularity tolerance of `lme4::isSingular()`).
#'
#' @keywords internal
#' @noRd
slopes_boundary_warning <- function(L, scale, names) {
  at_zero <- abs(diag(L)) * scale < 1e-4
  if (any(at_zero)) {
    warning("the random-effects covariance matrix is singular (on the boundary of the ",
            "parameter space): the variance of ", paste0("'", names[at_zero], "'", collapse = ", "),
            " given the preceding random effects is estimated at zero; consider a simpler ",
            "'random' formula. Standard errors use the Cholesky parameterisation, which ",
            "stays regular at the boundary", call. = FALSE)
  }
  invisible(at_zero)
}

#' Covariance of the parameters from the Hessian of the log-likelihood
#'
#' The inverse of the observed information; when it is singular (a
#' parameter not identified at the boundary), its Moore-Penrose generalized
#' inverse, which leaves the identified directions, among them the fixed
#' effects, finite.
#'
#' @keywords internal
#' @noRd
slopes_theta_vcov <- function(H) {
  V <- tryCatch(solve(-H), error = function(e) NULL)
  if (!is.null(V) && all(is.finite(V)) && all(diag(V) > 0)) return((V + t(V)) / 2)
  e <- eigen(-(H + t(H)) / 2, symmetric = TRUE)
  keep <- e$values > 1e-8 * max(e$values)
  V <- e$vectors[, keep, drop = FALSE] %*% (t(e$vectors[, keep, drop = FALSE]) / e$values[keep])
  (V + t(V)) / 2
}

#' Population-averaged mean: E over b of g^{-1}(eta + z'b)
#'
#' @keywords internal
#' @noRd
slopes_marginal_mean <- function(family, eta, Z, varcomp, nodes = 21L) {
  s2 <- as.numeric(rowSums((Z %*% varcomp) * Z))
  switch(family$family,
    gaussian = eta,
    poisson = , negbin = , gamma = exp(eta + s2 / 2),
    binomial = {
      gh <- gauss_hermite_quadrature(nodes)
      shifts <- sqrt(2) * outer(sqrt(s2), gh$nodes)
      as.numeric(linkinv_r(eta + shifts, family$link_id) %*% gh$weights)
    }
  )
}
