#' Individual effects of a panel model
#'
#' `fixef()` returns the estimated individual effects `alpha_i` of a
#' fixed-effects (`model = "within"`) fit, and `ranef()` the empirical
#' Bayes predictions of the random effects of a `model = "random"` fit.
#'
#' @section Fixed effects:
#' * Gaussian within fits: `alpha_i = ybar_i - xbar_i' beta`, the
#'   intercepts of the least-squares dummy-variable regression. They are
#'   not normalized: `fixest::fixef()` gives the same values for a model
#'   without an intercept.
#' * Poisson within fits: the conditional (Poisson fixed-effects) estimator
#'   does not estimate `alpha_i`, but the individual intercept that
#'   maximizes the unconditional likelihood at the conditional `beta` has
#'   the closed form `alpha_i = log(sum_t y_it / sum_t exp(x_it' beta))`,
#'   the values reported by `fixest::fepois()`. Individuals with no
#'   positive outcome (`alpha_i = -Inf`) are omitted, as in fixest.
#' * Negative binomial and gamma within fits: the dummy-variable
#'   intercepts of the joint maximum likelihood fit. The negative binomial
#'   fit drops the individuals whose outcome is always zero.
#' * Two-way fits (Gaussian, Poisson, negative binomial): a list with the
#'   `individual` and `time` effects, normalized so that the time effects
#'   average zero.
#' * The conditional logit (`family = "binomial"`, one-way or two-way) and
#'   the fixed-effects ordered logit condition the individual effects out
#'   of the likelihood and do not estimate them: `fixef()` stops with an
#'   error. The individual intercepts are not consistently estimable with
#'   a fixed number of periods anyway (incidental parameters).
#' * [panhurdle()] fits: the dummy-variable intercepts of the zero-truncated
#'   count part (individuals with at least one positive count); the binary
#'   part is a conditional logit.
#'
#' @section Random effects:
#' `ranef()` returns the posterior mean (`type = "mean"`, default) or the
#' posterior mode (`type = "mode"`) of each individual's random effects
#' given the individual's outcomes, at the estimated parameters (empirical
#' Bayes).
#' * Gaussian random intercept: the BLUP
#'   `sigma_mu^2 T_i / (sigma_mu^2 T_i + sigma_e^2) * mean_t(y_it - x_it' beta)`
#'   with the variance components of the fit (Swamy-Arora, or maximum
#'   likelihood for weighted fits). Mean and mode coincide.
#' * Gaussian random slopes: the BLUPs stored by the fit.
#' * Binomial, gamma, beta, Tobit and ordered random intercepts, which are
#'   normal on the link scale: the posterior mean is computed by adaptive
#'   Gauss-Hermite quadrature with the `R` nodes of the fit, centred at the
#'   posterior mode and scaled by the curvature there, the rule used for the
#'   likelihood. With `R = 1` (Laplace) the mean is the mode.
#' * Poisson and negative binomial random intercepts (Hausman, Hall and
#'   Griliches), whose mixing distribution is not normal: the effect is
#'   reported on the log scale of the multiplicative effect on the mean.
#'   For the Poisson-gamma model `nu_i ~ Gamma(d, d)` and the posterior is
#'   `Gamma(d + sum_t y_it, d + sum_t exp(x_it' beta))`, so the posterior
#'   mean of `log(nu_i)` is `digamma(d + Y_i) - log(d + L_i)` and its mode
#'   `log((d + Y_i) / (d + L_i))`. For the beta-mixed negative binomial
#'   model the effect is `log(delta_i)`, with `delta_i / (1 + delta_i)`
#'   distributed as `Beta(b, a)`, giving `digamma(b + Y_i) - digamma(a + L_i)`
#'   and `log((b + Y_i) / (a + L_i))`.
#' * Random slopes of the Poisson, negative binomial, gamma and binomial
#'   families: posterior modes are those stored by the fit; posterior means
#'   use the adaptive product rule of the fit (`R^q` nodes).
#'
#' Empirical Bayes predictions are shrunk towards zero, the more so the
#' fewer observations an individual has: their spread is smaller than the
#' estimated standard deviation of the random effect (for the Gaussian
#' model, `var(b_i) = sigma_mu^2 * shrinkage_i`), so their histogram is
#' narrower than the assumed normal density (see [plot.panglm()],
#' `which = "effects"`).
#'
#' @section Generics:
#' The generics `fixef()` and `ranef()` are also defined by nlme, whose
#' generics lme4 and fixest re-export. panglm exports its own generics so
#' that it does not depend on nlme, and when nlme is loaded it also
#' registers its methods for nlme's generics: whichever package is
#' attached last, `fixef(fit)` and `ranef(fit)` dispatch to the panglm
#' methods.
#'
#' @param object a fitted [panglm()] model (or a [panhurdle()] fit for
#'   `fixef()`)
#' @param type for `ranef()`, `"mean"` (posterior mean, default) or
#'   `"mode"` (posterior mode)
#' @param ... unused
#' @return `fixef()`: a named numeric vector of individual effects (names
#'   are the individual identifiers), or for two-way fits a list with
#'   components `individual` and `time`. `ranef()`: a matrix with one row
#'   per individual (row names are the identifiers) and one column per
#'   random effect, with attribute `"type"`.
#' @family panglm methods
#' @seealso [plot.panglm()] with `which = "effects"` to plot them.
#' @examples
#' data(copd)
#' fe <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
#'              model = "within")
#' head(fixef(fe))
#'
#' re <- panglm(fev1 ~ crp + treatment, data = copd,
#'              index = c("id", "visit"), model = "random")
#' head(ranef(re))
#'
#' hosp <- panglm(hospitalized ~ crp + treatment, data = copd,
#'                index = c("id", "visit"), model = "random",
#'                family = "binomial")
#' head(ranef(hosp))
#' @name fixef
#' @export
fixef <- function(object, ...) UseMethod("fixef")

#' @rdname fixef
#' @export
ranef <- function(object, ...) UseMethod("ranef")

#' @rdname fixef
#' @export
fixef.panglm <- function(object, ...) {
  if (!identical(object$model, "within")) {
    stop("fixef() needs a fixed-effects fit (model = 'within')",
         if (identical(object$model, "random")) "; use ranef() for random effects",
         call. = FALSE)
  }
  if (inherits(object, "panordered")) {
    stop("the fixed-effects ordered logit conditions the individual effects out of ",
         "the likelihood and does not estimate them", call. = FALSE)
  }
  if (identical(object$family$family, "binomial")) {
    stop("the conditional logit conditions the individual effects out of the ",
         "likelihood and does not estimate them",
         if (identical(object$effect, "twoways")) "; the time effects are in fit$time_effects",
         call. = FALSE)
  }
  alpha <- object$individual_effects
  if (is.null(alpha)) stop("the fit stores no individual effects", call. = FALSE)
  alpha <- alpha[is.finite(alpha)]
  if (identical(object$effect, "twoways")) {
    return(list(individual = alpha, time = object$time_effects))
  }
  alpha
}

#' @rdname fixef
#' @export
fixef.panglm_hurdle <- function(object, ...) {
  alpha <- object$count$individual_effects
  if (is.null(alpha)) stop("the fit stores no individual effects", call. = FALSE)
  if (!is.null(object$group_levels)) {
    names(alpha) <- as.character(object$group_levels[as.integer(names(alpha))])
  }
  alpha[is.finite(alpha)]
}

#' @rdname fixef
#' @export
ranef.panglm <- function(object, type = c("mean", "mode"), ...) {
  type <- match.arg(type)
  if (!identical(object$model, "random")) {
    stop("ranef() needs a random-effects fit (model = 'random')",
         if (identical(object$model, "within")) "; use fixef() for fixed effects",
         call. = FALSE)
  }
  ids <- as.character(object$group_levels)
  out <- panglm_random_effects(object, type)
  if (!is.matrix(out)) out <- matrix(out, ncol = 1L, dimnames = list(NULL, "(Intercept)"))
  rownames(out) <- ids
  attr(out, "type") <- type
  out
}

#' Empirical Bayes random effects (in panel order of the individuals)
#'
#' @keywords internal
#' @noRd
panglm_random_effects <- function(object, type) {
  fam <- object$family$family
  y <- object$y
  eta <- object$linear.predictors
  group_size <- object$group_size
  group <- rep(seq_along(group_size), group_size)
  if (!is.null(object$random)) {
    if (fam == "gaussian") return(object$ranef)
    if (type == "mode") return(object$ranef)
    logf <- panglm_conditional_logf(object)
    return(posterior_effects(logf, eta, group, object$Z, object$varcomp,
                             object$quadrature_nodes %||% 1L, start = object$ranef)$mean)
  }
  switch(fam,
    gaussian = {
      r <- as.numeric(rowsum(y - eta, group, reorder = FALSE)) / group_size
      s_mu2 <- object$sigma_mu2; s_v2 <- object$sigma_v2
      s_mu2 * group_size / (s_mu2 * group_size + s_v2) * r
    },
    poisson = {
      d <- object$dispersion_param
      A <- d + as.numeric(rowsum(y, group, reorder = FALSE))
      B <- d + as.numeric(rowsum(exp(eta), group, reorder = FALSE))
      if (type == "mean") digamma(A) - log(B) else log(A / B)
    },
    negbin = {
      A <- object$negbin_a + as.numeric(rowsum(exp(eta), group, reorder = FALSE))
      B <- object$negbin_b + as.numeric(rowsum(y, group, reorder = FALSE))
      if (type == "mean") digamma(B) - digamma(A) else log(B / A)
    },
    {
      logf <- panglm_conditional_logf(object)
      Z <- matrix(1, length(y), 1L)
      post <- posterior_effects(logf, eta, group, Z, matrix(object$sigma^2),
                                object$quadrature_nodes %||% 21L)
      post[[type]]
    }
  )
}

#' Conditional log-density of the outcome given the linear predictor
#'
#' Returns `function(eta)`: the vector of `log f(y_it | eta_it)` at the
#' fitted parameters, where `eta` includes the random effect.
#'
#' @keywords internal
#' @noRd
panglm_conditional_logf <- function(object) {
  fam <- object$family$family
  y <- object$y
  if (inherits(object, "panordered")) {
    lf <- ordered_link(object$family$link)
    cuts <- c(-Inf, object$cutpoints, Inf)
    upper <- cuts[y + 1L]; lower <- cuts[y]
    return(function(eta) {
      log(pmax(lf$F(upper - eta) - lf$F(lower - eta), 1e-300))
    })
  }
  switch(fam,
    binomial = {
      if (identical(object$family$link, "probit")) {
        function(eta) ifelse(y > 0.5, stats::pnorm(eta, log.p = TRUE),
                             stats::pnorm(-eta, log.p = TRUE))
      } else {
        function(eta) ifelse(y > 0.5, stats::plogis(eta, log.p = TRUE),
                             stats::plogis(-eta, log.p = TRUE))
      }
    },
    poisson = function(eta) stats::dpois(y, exp(eta), log = TRUE),
    negbin = {
      theta <- object$negbin_theta
      function(eta) stats::dnbinom(y, size = theta, mu = exp(eta), log = TRUE)
    },
    gamma = {
      nu <- object$shape
      function(eta) stats::dgamma(y, shape = nu, rate = nu * exp(-eta), log = TRUE)
    },
    beta = {
      lf <- beta_link(object$family$link)
      phi <- object$precision
      function(eta) {
        mu <- pmin(pmax(lf$mu(eta), 1e-12), 1 - 1e-12)
        stats::dbeta(y, mu * phi, (1 - mu) * phi, log = TRUE)
      }
    },
    tobit = {
      left <- object$family$left; right <- object$family$right
      se <- object$sigma_e
      status <- tobit_status(y, left, right)
      function(eta) {
        out <- stats::dnorm(y, eta, se, log = TRUE)
        lo <- status == -1L; hi <- status == 1L
        out[lo] <- stats::pnorm((left - eta[lo]) / se, log.p = TRUE)
        out[hi] <- stats::pnorm((eta[hi] - right) / se, log.p = TRUE)
        out
      }
    },
    stop("no conditional density for family '", fam, "'", call. = FALSE)
  )
}

#' Posterior modes and means of normal random effects
#'
#' `b_i ~ N(0, Sigma)` enters the linear predictor as `eta_it + z_it' b_i`.
#' The modes are found by Newton's method on
#' `h_i(b) = sum_t log f(y_it | eta_it + z_it' b) - b' Sigma^-1 b / 2`
#' (numerical derivatives of `logf` in the linear predictor). The means use
#' adaptive Gauss-Hermite product quadrature with `R` nodes per dimension,
#' centred at the mode and scaled by the Cholesky factor of the inverse
#' curvature there.
#'
#' @keywords internal
#' @noRd
posterior_effects <- function(logf, eta, group, Z, Sigma, R, start = NULL,
                              maxit = 100L, tol = 1e-10) {
  G <- max(group); q <- ncol(Z)
  Z <- unname(as.matrix(Z))
  Sigma <- as.matrix(Sigma)
  # a (near) singular Sigma: regularize the prior precision
  Prec <- solve(Sigma + diag(1e-10 * max(1, max(diag(Sigma))), q))
  b <- if (is.null(start)) matrix(0, G, q) else unname(as.matrix(start))
  rs <- function(v) as.numeric(rowsum(v, group, reorder = FALSE))
  derivs <- function(b) {
    e <- eta + rowSums(Z * b[group, , drop = FALSE])
    h <- 1e-4
    lp <- logf(e + h); l0 <- logf(e); lm <- logf(e - h)
    d1 <- (lp - lm) / (2 * h)
    d2 <- pmin((lp - 2 * l0 + lm) / h^2, 0)
    grad <- vapply(seq_len(q), function(a) rs(d1 * Z[, a]), numeric(G)) -
      b %*% Prec
    if (!is.matrix(grad)) grad <- matrix(grad, nrow = G)
    A <- array(0, c(G, q, q))
    for (a in seq_len(q)) for (c in seq_len(a)) {
      v <- -rs(d2 * Z[, a] * Z[, c]) + Prec[a, c]
      A[, a, c] <- v; A[, c, a] <- v
    }
    list(grad = grad, A = A)
  }
  for (iter in seq_len(maxit)) {
    dv <- derivs(b)
    step <- bmatvec(bspdinv(dv$A), dv$grad)
    step[!is.finite(step)] <- 0
    step <- pmax(-2, pmin(2, step))
    b <- b + step
    if (max(abs(step)) < tol) break
  }
  dv <- derivs(b)
  # Cholesky factor of the posterior covariance approximation A^-1
  C <- bchol(bspdinv(dv$A))
  gh <- gauss_hermite_quadrature(max(1L, as.integer(R)))
  grid <- as.matrix(expand.grid(rep(list(seq_along(gh$nodes)), q)))
  x <- matrix(gh$nodes[grid], ncol = q)
  logw <- rowSums(matrix(log(gh$weights[grid]), ncol = q)) + rowSums(x^2)
  h_at <- function(bb) {
    e <- eta + rowSums(Z * bb[group, , drop = FALSE])
    rs(logf(e)) - 0.5 * rowSums((bb %*% Prec) * bb)
  }
  h0 <- h_at(b)
  nodes <- nrow(x)
  logterm <- matrix(0, G, nodes)
  sums <- matrix(0, G, q)
  for (j in seq_len(nodes)) {
    bj <- b + sqrt(2) * bmatvec(C, matrix(x[j, ], G, q, byrow = TRUE))
    logterm[, j] <- logw[j] + h_at(bj) - h0
  }
  top <- apply(logterm, 1L, max)
  wts <- exp(logterm - top)
  wts <- wts / rowSums(wts)
  for (j in seq_len(nodes)) {
    bj <- b + sqrt(2) * bmatvec(C, matrix(x[j, ], G, q, byrow = TRUE))
    sums <- sums + wts[, j] * bj
  }
  list(mode = b, mean = sums)
}

#' Register the fixef() and ranef() methods for nlme's generics
#'
#' nlme (re-exported by lme4 and fixest) defines generics of the same
#' names; registering the methods there keeps `fixef(fit)` and `ranef(fit)`
#' working when one of these packages is attached after panglm.
#'
#' @keywords internal
#' @noRd
register_nlme_methods <- function(...) {
  ns <- asNamespace("nlme")
  registerS3method("fixef", "panglm", fixef.panglm, envir = ns)
  registerS3method("fixef", "panglm_hurdle", fixef.panglm_hurdle, envir = ns)
  registerS3method("ranef", "panglm", ranef.panglm, envir = ns)
  invisible(NULL)
}

.onLoad <- function(libname, pkgname) {
  if (isNamespaceLoaded("nlme")) register_nlme_methods()
  setHook(packageEvent("nlme", "onLoad"), register_nlme_methods)
  invisible(NULL)
}
