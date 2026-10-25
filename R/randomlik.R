#' Individual-level log-likelihoods of the random-effects models
#'
#' One contribution per individual, written in R with the parameterization
#' of the compiled estimators: `theta = (beta, log d)` for Poisson-gamma,
#' `(beta, log a, log b)` for the beta-negative-binomial, `(beta, log sigma)`
#' for the binomial, `(beta, log sigma, log shape)` for the gamma and
#' `(beta, log sigma_e, log sigma_u)` for the Gaussian model. They drive
#' weighted random-effects fits (sum over individuals of `w_i * ll_i`) and
#' cluster-robust covariances (individual scores by central differences).
#'
#' @keywords internal
#' @noRd
random_contrib <- function(family, X, y, group, R = 21L) {
  k <- ncol(X)
  Ti <- tabulate(group)
  switch(family$family,
    poisson = {
      base <- y * 0 - lgamma(y + 1)
      Yi <- as.numeric(rowsum(y, group))
      function(theta, score) {
        eta <- as.numeric(X %*% theta[seq_len(k)])
        d <- exp(theta[k + 1L])
        Li <- as.numeric(rowsum(exp(eta), group))
        lnA <- as.numeric(rowsum(y * eta + base, group))
        list(ll = lnA + d * log(d) - (d + Yi) * log(Li + d) + lgamma(d + Yi) - lgamma(d))
      }
    },
    negbin = {
      Yi <- as.numeric(rowsum(y, group))
      lfy <- lgamma(y + 1)
      function(theta, score) {
        l <- exp(as.numeric(X %*% theta[seq_len(k)]))
        a <- exp(theta[k + 1L]); b <- exp(theta[k + 2L])
        Li <- as.numeric(rowsum(l, group))
        lnA <- as.numeric(rowsum(lgamma(l + y) - lgamma(l) - lfy, group))
        # lbeta stays accurate when a and b are large (near the Poisson limit)
        list(ll = lnA + lbeta(a + Li, b + Yi) - lbeta(a, b))
      }
    },
    binomial = {
      gh <- gauss_hermite_quadrature(R)
      logcdf <- if (family$link == "logit") {
        function(z) stats::plogis(z, log.p = TRUE)
      } else {
        function(z) stats::pnorm(z, log.p = TRUE)
      }
      sgn <- 2 * y - 1
      function(theta, score) {
        eta <- as.numeric(X %*% theta[seq_len(k)])
        sigma <- exp(theta[k + 1L])
        L <- vapply(seq_along(gh$nodes), function(q) {
          log(gh$weights[q]) +
            as.numeric(rowsum(logcdf(sgn * (eta + sqrt(2) * sigma * gh$nodes[q])), group))
        }, numeric(length(Ti)))
        if (!is.matrix(L)) L <- matrix(L, nrow = length(Ti))
        list(ll = log_row_sum_exp(L))
      }
    },
    gamma = random_gamma_contrib(X, y, group, R),
    gaussian = {
      function(theta, score) {
        e <- y - as.numeric(X %*% theta[seq_len(k)])
        s2e <- exp(2 * theta[k + 1L]); s2u <- exp(2 * theta[k + 2L])
        se <- as.numeric(rowsum(e, group)); see <- as.numeric(rowsum(e^2, group))
        quad <- (see - s2u / (s2e + Ti * s2u) * se^2) / s2e
        logdet <- (Ti - 1) * log(s2e) + log(s2e + Ti * s2u)
        list(ll = -0.5 * (Ti * log(2 * pi) + logdet + quad))
      }
    },
    stop("no individual likelihood for family '", family$family, "'", call. = FALSE)
  )
}

log_row_sum_exp <- function(L) {
  top <- apply(L, 1, max)
  top + log(rowSums(exp(L - top)))
}

#' Gamma random-intercept contributions by adaptive Gauss-Hermite quadrature
#'
#' @keywords internal
#' @noRd
random_gamma_contrib <- function(X, y, group, R) {
  k <- ncol(X); G <- max(group)
  gh <- gauss_hermite_quadrature(R)
  log_w <- log(gh$weights * sqrt(pi)) + gh$nodes^2
  Ti <- tabulate(group)
  slogy <- as.numeric(rowsum(log(y), group))
  function(theta, score) {
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
    u <- rep(0, G)
    for (step in seq_len(100L)) {
      g1 <- -nu * Ti + nu * sy * exp(-u) - u / sigma^2
      g2 <- -nu * sy * exp(-u) - 1 / sigma^2
      delta <- pmax(-5, pmin(5, g1 / g2))
      if (anyNA(delta)) return(list(ll = rep(-Inf, G)))
      u <- u - delta
      if (max(abs(delta)) < 1e-10) break
    }
    tau <- 1 / sqrt(nu * sy * exp(-u) + 1 / sigma^2)
    terms <- vapply(seq_along(gh$nodes), function(q) {
      log_w[q] + h(u + sqrt(2) * tau * gh$nodes[q])
    }, numeric(G))
    if (!is.matrix(terms)) terms <- matrix(terms, nrow = G)
    list(ll = log(sqrt(2) * tau) + log_row_sum_exp(terms))
  }
}

#' Natural-scale parameter vector of a fitted random-effects model
#'
#' @keywords internal
#' @noRd
random_theta <- function(fit, family) {
  b <- unname(fit$coefficients)
  switch(family$family,
    poisson = c(b, log(fit$dispersion_param)),
    negbin = c(b, log(fit$negbin_a), log(fit$negbin_b)),
    binomial = c(b, log(fit$sigma)),
    gamma = c(b, log(fit$sigma), log(fit$shape)),
    gaussian = c(b, 0.5 * log(fit$sigma_v2), 0.5 * log(max(fit$sigma_mu2, 1e-8)))
  )
}

#' Generic starting values: pooled GLM slopes and moderate variance components
#'
#' @keywords internal
#' @noRd
random_default_theta <- function(X, y, family, maxit, tol) {
  fam_id <- switch(family$family, poisson = 1L, negbin = 1L, binomial = 2L, gamma = 3L, 0L)
  link_id <- if (fam_id == 0L) 0L else if (fam_id == 2L) family$link_id else 1L
  b <- as.numeric(irls_fit_cpp(X, y, fam_id, link_id, maxit, tol)$coefficients)
  switch(family$family,
    poisson = c(b, 0),
    negbin = c(b, log(2), 0),
    binomial = c(b, log(0.5)),
    gamma = c(b, log(0.5), 0),
    gaussian = c(b, rep(0.5 * log(stats::var(y) / 2), 2))
  )
}

#' Weighted random-effects fit: maximize sum_i w_i ll_i
#'
#' Starts from the unweighted estimate. The Gaussian model is estimated by
#' maximum likelihood here (the unweighted fit uses Swamy-Arora FGLS).
#'
#' @keywords internal
#' @noRd
fit_random_weighted <- function(X, y, family, group_start, group_size, R, maxit, tol, weights,
                                compiled_start = TRUE) {
  w <- individual_weights(weights, group_size)
  # individuals with zero weight carry no information: estimate without them
  keep_ind <- w > 0
  keep_obs <- rep(keep_ind, group_size)
  Xk <- X[keep_obs, , drop = FALSE]; yk <- y[keep_obs]
  size_k <- group_size[keep_ind]
  start_k <- c(0L, cumsum(size_k)[-length(size_k)])
  group <- rep(seq_along(size_k), size_k)
  start <- if (compiled_start) {
    tryCatch(random_theta(fit_random(Xk, yk, family, start_k, size_k, R, maxit, tol), family),
             error = function(e) NULL)
  }
  if (is.null(start) || any(!is.finite(start))) start <- random_default_theta(Xk, yk, family, maxit, tol)
  contrib <- random_contrib(family, Xk, yk, group, R)
  res <- ml_optimize(start, contrib, w[keep_ind], maxit, tol)
  score <- matrix(0, length(w), length(res$theta))
  score[keep_ind, ] <- res$score
  res$score <- score
  k <- ncol(X)
  coefs <- res$theta[seq_len(k)]
  names(coefs) <- colnames(X)
  vcov <- res$vcov[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(vcov) <- list(colnames(X), colnames(X))
  eta <- as.numeric(X %*% coefs)
  extra <- res$theta[-seq_len(k)]
  out <- list(coefficients = coefs, vcov = vcov, loglik = res$loglik, theta = res$theta,
              theta_vcov = res$vcov, theta_score = res$score, unit_weights = w,
              df.residual = sum(w > 0) - length(res$theta), iterations = res$iterations,
              convergence = res$convergence, quadrature_nodes = R)
  switch(family$family,
    poisson = {
      out$dispersion_param <- exp(extra[1])
      out$fitted.values <- exp(eta)
    },
    negbin = {
      a <- exp(extra[1]); b <- exp(extra[2])
      out$negbin_a <- a; out$negbin_b <- b
      out$marginal_mean_multiplier <- if (a > 1) b / (a - 1) else NA_real_
      out$fitted.values <- out$marginal_mean_multiplier * exp(eta)
    },
    binomial = {
      out$sigma <- exp(extra[1]); out$dispersion <- 1
      out$fitted.values <- random_binomial_marginal_mean(eta, out$sigma, family, R)
    },
    gamma = {
      out$sigma <- exp(extra[1]); out$shape <- exp(extra[2]); out$dispersion <- 1 / out$shape
      out$fitted.values <- exp(eta + out$sigma^2 / 2)
    },
    gaussian = {
      out$sigma_v2 <- exp(2 * extra[1]); out$sigma_mu2 <- exp(2 * extra[2])
      out$dispersion <- out$sigma_v2
      out$fitted.values <- eta
    }
  )
  out
}

#' Cluster-robust covariance for a random-effects panglm fit
#'
#' Likelihood fits: sandwich of the individual scores of the integrated
#' likelihood (clusters must contain whole individuals), with the small-sample
#' factor G / (G - 1). The unweighted Gaussian (Swamy-Arora) fit uses the
#' cluster sandwich of the quasi-demeaned regression, as in Stata's
#' `xtreg, re vce(cluster)`.
#'
#' @keywords internal
#' @noRd
robust_vcov_random <- function(object, cluster) {
  cl_obs <- panglm_cluster_vector(object, cluster)
  cl <- individual_clusters(cl_obs, object$group_size)
  k <- length(object$coefficients)
  weighted <- !is.null(object$unit_weights)
  if (object$family$family == "gaussian" && !weighted) {
    qd <- quasi_demean_cpp(object$X, object$y, object$group_start, object$group_size,
                           object$theta)
    bread <- solve(crossprod(qd$X))
    resid <- as.numeric(qd$y - qd$X %*% object$coefficients)
    group <- rep(seq_along(object$group_size), object$group_size)
    score <- rowsum(resid * qd$X, group)
    v <- ml_sandwich(bread, score, rep(1, nrow(score)), cl, clusters_only = TRUE)
    dimnames(v) <- dimnames(object$vcov)
    return(v)
  }
  group <- rep(seq_along(object$group_size), object$group_size)
  if (weighted) {
    w <- object$unit_weights
    score <- object$theta_score
    bread <- object$theta_vcov
  } else {
    w <- rep(1, length(object$group_size))
    contrib <- random_contrib(object$family, object$X, object$y, group,
                              object$quadrature_nodes %||% 21L)
    theta <- random_theta(object, object$family)
    score <- ml_scores(contrib, theta)
    bread <- solve(-ml_hessian(contrib, theta, w))
  }
  v <- ml_sandwich(bread, score, w, cl, clusters_only = TRUE)
  v <- v[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(v) <- dimnames(object$vcov)
  v
}
