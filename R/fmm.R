#' Latent class (finite mixture) models for panel data
#'
#' Finite mixture regression in which every individual belongs to one of
#' `k` latent classes for all of its observations, the panel latent class
#' model of Deb and Trivedi (1997, 2002) and Bago d'Uva (2006). Within class
#' `c` the outcome follows a GLM with class-specific coefficients `beta_c`
#' (and dispersion), and the individual likelihood is
#' `L_i = sum_c pi_ic prod_t f_c(y_it | x_it)`. Class probabilities are
#' constant (`pi_c`) or depend on individual-level concomitant variables
#' through a multinomial logit. In health economics the classes typically
#' separate a large group of low users of care from a smaller group of high
#' users with a different response to income, insurance or health status.
#'
#' @section Estimation:
#' The EM algorithm is run from `nstart` starting partitions of the
#' individuals: the first splits them by the quantiles of their mean
#' outcome, the others are random. The E step computes the posterior class
#' probabilities `tau_ic` of each individual from the whole history
#' `y_i1, ..., y_iT`; the M step fits one weighted GLM per class (weights
#' `tau_ic` repeated over the individual's observations, iteratively
#' reweighted least squares warm-started at the previous iteration's
#' coefficients; NB2 by IRLS and a Newton step for `theta`) and the class
#' probabilities (sample means of `tau_ic`, or a weighted multinomial logit
#' on the concomitant variables). The run with the largest log-likelihood is
#' then polished by BFGS and Newton steps on the observed-data
#' log-likelihood `sum_i log L_i`, whose numerical Hessian (from analytic
#' scores) gives the classical covariance `solve(-H)`. The cluster
#' covariance is the sandwich `B (sum_g s_g s_g') B` with `B = solve(-H)`,
#' the individual scores `s_i` summed within clusters (the individual by
#' default) and the factor `G / (G - 1)`, the convention used in the rest of
#' the package for likelihoods whose units are individuals.
#'
#' @section Weights:
#' Class membership is a property of the individual, so `weights` are
#' individual (sampling) weights: they must be constant within individuals
#' (one weight per individual, repeated over its rows) and `panfmm()` stops
#' otherwise. The fit maximizes the weighted pseudo log-likelihood
#' `sum_i w_i log L_i`, the convention of flexmix (whose integer weights must
#' also be identical within groups). In the EM algorithm the M step weights
#' become `w_i tau_ic`, and the class probabilities are weighted means of
#' the posteriors. The classical covariance is `solve(-H_w)` from the
#' weighted Hessian, which is correct for frequency weights; for sampling
#' weights use `vcov = "cluster"`, the sandwich `B (sum_g S_g S_g') B` with
#' the weighted individual scores `S_i = w_i s_i` summed within clusters.
#' Individuals with zero weight do not contribute. With weights the
#' reported `loglik` (hence AIC and BIC) is the weighted log-likelihood, and
#' BIC keeps `n` the number of observations (flexmix uses the sum of the
#' weights, which differs for sampling weights scaled to population totals).
#'
#' Ancillary parameters are estimated on the log scale: `log(sigma)` for the
#' Gaussian family and `log(theta)` for NB2 (variance
#' `mu + mu^2 / theta`, i.e. `alpha = 1 / theta`). The class probabilities
#' are parameterized by the multinomial logit coefficients of classes
#' `2, ..., k` against class 1 (an intercept only when there are no
#' concomitant variables), and these are the parameters reported in
#' `coef()` and `vcov()`; the implied probabilities are in `fit$prior`.
#'
#' @section Label switching:
#' Mixture likelihoods are invariant to permutations of the class labels.
#' `panfmm()` orders the classes by ascending mean outcome: the mean over
#' all observations of the class-specific fitted mean `mu_c(x_it)`, so that
#' class 1 is the low-use class and class `k` the high-use class. The
#' multinomial logit of the class probabilities uses class 1 as reference.
#'
#' @param formula a model formula `y ~ x1 + x2` for the class-specific
#'   regressions (the same regressors in every class). It may contain
#'   `offset()` terms, e.g. `offset(log(exposure))` for counts over
#'   different exposure periods; the offset enters the linear predictor of
#'   every class with coefficient 1.
#'   `lag(x, k)`, `lead(x, k)` and `diff(x)` are taken within each
#'   individual on the grid of sorted distinct time values (see
#'   "Lags and differences" in [panglm()]; needs
#'   `index = c(id, time)`); the first waves of each individual get `NA`
#'   and leave the sample, and `lag(x, 0:2)` gives several lags
#' @param data a data.frame in long format
#' @param index length-1 or length-2 character vector giving the individual
#'   (and time) columns; class membership is constant within individual. Use
#'   a row identifier for cross-section data
#' @param k number of latent classes
#' @param family `"poisson"`, `"negbin"` (NB2), `"gaussian"`, `"binomial"`
#'   (logit) or `"probit"` (binomial with probit link), or the corresponding
#'   panglm family spec, e.g. `binomialfamily("probit")`
#' @param concomitant optional one-sided formula `~ z1 + z2` of
#'   individual-level (time-invariant) variables for the class
#'   probabilities
#' @param nstart number of EM starts (the first deterministic, the others
#'   random)
#' @param seed optional random seed for the random starts; the global
#'   random number stream is restored afterwards
#' @param vcov `"classical"` (inverse negative Hessian) or `"cluster"`
#'   (sandwich clustered by individual); see [vcov.panfmm()]
#' @param maxit maximum number of EM iterations per start
#' @param tol convergence tolerance on the relative change of the EM
#'   log-likelihood
#' @param na.action missing-data action passed to [stats::model.frame()]
#' @param weights optional nonnegative individual (sampling) weights: the
#'   name of a column of `data` or a numeric vector with one value per row
#'   of `data`, constant within individuals; see the Weights section
#' @param offset optional numeric offset (the name of a column of `data` or
#'   a vector with one value per row of `data`), added
#'   to the linear predictor of every class (in addition to any `offset()`
#'   terms of `formula`), as in [stats::glm()]
#' @return an object of class `"panfmm"`, a list with components
#'   `coefficients` (all parameters, named `class<c>:<term>` and
#'   `prior<c>:<term>`), `vcov`, `vcov_type`, `classcoef` (a matrix of
#'   class-specific regression coefficients, one column per class), `sigma`
#'   (Gaussian) or `theta` (NB2) by class, `prior` (the class probabilities
#'   of each individual, one row per individual; with no concomitant
#'   variables all rows are equal), `prior_mean` (their averages),
#'   `posterior` (posterior class probabilities, one row per individual,
#'   rownames the individual identifiers), `class` (the modal posterior
#'   class), `classmean` (the mean outcome of each class, used for the label
#'   order), `loglik`, `df`, `starts` (log-likelihood of each EM start),
#'   `nobs`, `n_groups`, `weights` (the individual weights, `NULL` without
#'   weights), `offset` (the total offset of each observation in panel
#'   order, `NULL` without offsets), `k`, `family`, `call` and `formula`. Methods:
#'   `print`, `summary`, `coef`, `vcov`, `confint`, `predict`, `logLik`
#'   (hence `AIC` and `BIC`, with `n` the number of observations), `nobs`,
#'   `tidy` and `glance`.
#' @references
#' Deb, P. and Trivedi, P. K. (1997). Demand for medical care by the
#' elderly: a finite mixture approach. *Journal of Applied Econometrics*,
#' 12, 313-336.
#'
#' Deb, P. and Trivedi, P. K. (2002). The structure of demand for health
#' care: latent class versus two-part models. *Journal of Health
#' Economics*, 21, 601-625. \doi{10.1016/S0167-6296(02)00008-5}
#'
#' Bago d'Uva, T. (2006). Latent class models for utilisation of health
#' care. *Health Economics*, 15, 329-343. \doi{10.1002/hec.1112}
#'
#' Gruen, B. and Leisch, F. (2008). FlexMix version 2: finite mixtures with
#' concomitant variables and varying and constant parameters. *Journal of
#' Statistical Software*, 28(4). \doi{10.18637/jss.v028.i04}
#' @seealso [summary.panfmm()], [vcov.panfmm()] and [predict.panfmm()] for
#'   the fitted model, [tidy-methods] for `tidy()` and `glance()`.
#' @family panel models
#' @section Validation:
#' Checked against flexmix in the package's validation suite on GitHub:
#' `flexmix::flexmix()` with the `| id` grouping, started from the same
#' posterior, gives the same log-likelihood, class and concomitant
#' coefficients, class probabilities, posteriors and BIC, and
#' `flexmix::refit()` the same standard errors, for Gaussian, Poisson, NB2 and
#' binomial components; flexmix from its own random starts finds no better
#' optimum.
#' @examples
#' # Doctor visits of 400 patients over 4 years: 70% low users and 30%
#' # high users with a stronger effect of chronic illness
#' set.seed(1)
#' N <- 400; TT <- 4
#' d <- data.frame(id = rep(seq_len(N), each = TT), year = rep(seq_len(TT), N))
#' high <- rep(rbinom(N, 1, 0.3), each = TT)
#' d$chronic <- rbinom(N * TT, 1, 0.4)
#' d$age <- rep(round(runif(N, 40, 80)), each = TT)
#' mu <- exp(ifelse(high == 1, 1.5 + 0.6 * d$chronic, 0 + 0.3 * d$chronic) +
#'           0.01 * (d$age - 60))
#' d$visits <- rpois(N * TT, mu)
#'
#' fit <- panfmm(visits ~ chronic + age, data = d, index = c("id", "year"),
#'               k = 2, family = "poisson", nstart = 3, seed = 1)
#' summary(fit)
#' head(fit$posterior)
#' table(predict(fit, type = "class"))
#' @export
panfmm <- function(formula, data, index, k = 2, family = "poisson",
                   concomitant = NULL, nstart = 10, seed = NULL,
                   vcov = c("classical", "cluster"), maxit = 1000, tol = 1e-6,
                   na.action = stats::na.omit, weights = NULL, offset = NULL) {
  vcov <- match.arg(vcov)
  family <- resolve_family(family)
  weights <- data_column(weights, data, "weights")
  offset <- data_column(offset, data, "offset")
  if (!family$family %in% c("gaussian", "poisson", "negbin", "binomial")) {
    stop("panfmm() supports family = 'gaussian', 'poisson', 'negbin' and 'binomial'",
         call. = FALSE)
  }
  link <- if (family$family == "binomial") family$link %||% "logit" else ""
  if (family$family == "binomial" && !link %in% c("logit", "probit")) {
    stop("panfmm() supports the logit and probit links for family = 'binomial'",
         call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || k < 1 || k != round(k)) {
    stop("'k' must be a positive integer", call. = FALSE)
  }
  k <- as.integer(k)
  if (!is.numeric(nstart) || length(nstart) != 1L || nstart < 1) {
    stop("'nstart' must be a positive integer", call. = FALSE)
  }
  if (missing(index)) stop("'index' is required, e.g. index = c(\"id\", \"time\")", call. = FALSE)
  if (!is.character(index) || !length(index) %in% 1:2 || anyDuplicated(index)) {
    stop("'index' must contain one or two distinct column names", call. = FALSE)
  }
  missing_index <- setdiff(index, names(data))
  if (length(missing_index)) {
    stop("index column(s) not found in data: ", paste(missing_index, collapse = ", "),
         call. = FALSE)
  }

  stop_panel_calls(concomitant, what = "the 'concomitant' formula")
  formula <- panel_guard_formula(formula)
  parts <- fmm_formulas(formula, concomitant)
  mf <- panel_model_frame(parts$combined, data, index, na.action = na.action)
  model_rows <- match(rownames(mf), rownames(data))
  if (anyNA(model_rows)) stop("could not align model-frame rows with data", call. = FALSE)
  y <- stats::model.response(mf)
  if (family$family == "binomial") y <- binomial_response_to_numeric(y)
  if (!is.numeric(y) || any(!is.finite(y))) {
    stop("the response must be finite and numeric", call. = FALSE)
  }
  if (family$family %in% c("poisson", "negbin") && (any(y < 0) || any(y != round(y)))) {
    stop("count families need a non-negative integer response", call. = FALSE)
  }
  X <- stats::model.matrix(parts$regressors, data = mf)
  off <- stats::model.offset(mf)
  if (!is.null(offset)) {
    if (!is.numeric(offset) || length(offset) != nrow(data) || any(!is.finite(offset))) {
      stop("'offset' must be a column name or a finite numeric vector with nrow(data) values", call. = FALSE)
    }
    off <- (off %||% 0) + offset[model_rows]
  }
  if (!is.null(off) && any(!is.finite(off))) stop("the offset must be finite", call. = FALSE)
  Zobs <- stats::model.matrix(parts$concomitant, data = mf)

  panel_cols <- data[model_rows, index, drop = FALSE]
  if (anyNA(panel_cols)) stop("panel index columns must not contain missing values", call. = FALSE)
  panel <- build_panel_index(cbind(panel_cols, .panfmm_row = seq_len(nrow(mf))), index)
  ord <- panel$data$.panfmm_row
  y <- as.numeric(y[ord])
  X <- X[ord, , drop = FALSE]
  Zobs <- Zobs[ord, , drop = FALSE]
  if (!is.null(off)) {
    off <- rep_len(as.numeric(off), nrow(mf))[ord]
    if (all(off == 0)) off <- NULL
  }
  group <- rep.int(seq_len(panel$n_groups), panel$group_size)
  Z <- Zobs[panel$group_start + 1L, , drop = FALSE]
  if (max(abs(Zobs - Z[group, , drop = FALSE])) > 1e-10 * max(1, abs(Zobs))) {
    stop("the concomitant variables must be constant within individuals", call. = FALSE)
  }
  rownames(Z) <- NULL
  if (qr(X)$rank < ncol(X)) stop("the regressors are collinear", call. = FALSE)
  if (qr(Z)$rank < ncol(Z)) stop("the concomitant variables are collinear", call. = FALSE)
  if (k > 1L && panel$n_groups <= k) {
    stop("panfmm() needs more individuals than classes", call. = FALSE)
  }
  w <- NULL
  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != nrow(data) ||
        any(!is.finite(weights)) || any(weights < 0)) {
      stop("'weights' must be a column name or a finite nonnegative numeric vector with nrow(data) values",
           call. = FALSE)
    }
    w <- fmm_individual_weights(weights[model_rows][ord], panel$group_size)
    if (!any(w > 0)) stop("all weights are zero", call. = FALSE)
    if (all(w == 1)) w <- NULL
  }

  spec <- fmm_spec(X, Z, y, group, family$family, k, w, off, link)
  em <- fmm_em_starts(spec, nstart, seed, maxit, tol)
  par <- em$par
  for (pass in 1:3) {
    par <- fmm_reorder(spec, par)
    contrib <- function(theta, score) fmm_contrib(spec, theta, score)
    opt <- ml_optimize(par, contrib, spec$w, maxit = 500, tol = 1e-12)
    par <- opt$theta
    if (identical(fmm_class_order(spec, par), seq_len(k))) break
  }
  names(par) <- spec$names
  vc <- opt$vcov
  dimnames(vc) <- list(spec$names, spec$names)
  state <- fmm_contrib(spec, par, TRUE)
  score <- state$score
  colnames(score) <- spec$names

  fit <- fmm_unpack(spec, par)
  fit$coefficients <- par
  fit$vcov <- vc
  fit$vcov_classical <- vc
  fit$vcov_type <- "classical"
  fit$hessian <- opt$hessian
  fit$score <- score
  post <- state$tau
  ids <- as.character(panel$group_id)
  dimnames(post) <- list(ids, paste0("class", seq_len(k)))
  fit$posterior <- post
  fit$class <- stats::setNames(max.col(post, ties.method = "first"), ids)
  rownames(fit$prior) <- ids
  fit$loglik <- opt$loglik
  fit$df <- length(par)
  fit$starts <- em$starts
  fit$em_iterations <- em$iterations
  fit$convergence <- opt$convergence
  fit$k <- k
  fit$family <- family
  fit$call <- match.call()
  fit$formula <- formula
  fit$concomitant <- concomitant
  fit$terms <- stats::terms(parts$regressors)
  fit$terms_concomitant <- stats::terms(parts$concomitant)
  fit$xlevels <- stats::.getXlevels(stats::terms(parts$combined), mf)
  fit$index <- index
  fit$nobs <- length(y)
  fit$n_groups <- panel$n_groups
  fit$group_size <- panel$group_size
  fit$weights <- w
  fit$offset <- off
  fit$offset_argument <- !is.null(offset)
  fit$y <- y
  fit$X <- X
  fit$Z <- Z
  fit$cluster_id <- rep(panel$group_id, panel$group_size)
  fit$panel_order <- ord
  fit$na.action <- attr(mf, "na.action")
  class(fit) <- "panfmm"
  if (vcov == "cluster") {
    fit$vcov <- vcov.panfmm(fit, type = "cluster")
    fit$vcov_type <- "cluster"
  }
  fit
}

#' Regression, concomitant and combined formulas
#'
#' @keywords internal
#' @noRd
fmm_formulas <- function(formula, concomitant) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("'formula' must be a two-sided formula 'y ~ x1 + x2'", call. = FALSE)
  }
  env <- environment(formula) %||% parent.frame()
  rhs <- formula[[3L]]
  if (is.null(concomitant)) {
    conc <- stats::as.formula(~ 1, env = env)
    combined <- formula
  } else {
    if (!inherits(concomitant, "formula") || length(concomitant) != 2L) {
      stop("'concomitant' must be a one-sided formula '~ z1 + z2'", call. = FALSE)
    }
    if (!is.null(attr(stats::terms(concomitant), "offset"))) {
      stop("offsets are not allowed in 'concomitant'", call. = FALSE)
    }
    conc <- concomitant
    combined <- stats::as.formula(call("~", formula[[2L]], call("+", rhs, concomitant[[2L]])),
                                  env = env)
  }
  list(regressors = stats::as.formula(call("~", rhs), env = env),
       concomitant = conc, combined = combined)
}

#' Model description shared by the EM and the likelihood
#'
#' @keywords internal
#' @noRd
fmm_spec <- function(X, Z, y, group, family, k, w = NULL, offset = NULL, link = "") {
  p <- ncol(X)
  q <- ncol(Z)
  aux <- family %in% c("gaussian", "negbin")
  aux_name <- if (family == "gaussian") "log(sigma)" else "log(theta)"
  m <- p + aux
  cls <- unlist(lapply(seq_len(k), function(c) {
    paste0("class", c, ":", c(colnames(X), if (aux) aux_name))
  }))
  pri <- if (k > 1L) {
    unlist(lapply(2:k, function(c) paste0("prior", c, ":", colnames(Z))))
  } else character(0)
  list(X = X, Z = Z, y = y, group = group, family = family, k = k, p = p, q = q,
       aux = aux, m = m, N = nrow(Z), names = c(cls, pri),
       w = if (is.null(w)) rep(1, nrow(Z)) else w,
       offset = if (is.null(offset)) rep(0, length(y)) else offset,
       has_offset = !is.null(offset), probit = identical(link, "probit"),
       lgy = if (family == "negbin") lgamma(y + 1) else NULL)
}

#' Individual weights from observation weights constant within individuals
#'
#' @keywords internal
#' @noRd
fmm_individual_weights <- function(weights, group_size) {
  group <- rep.int(seq_along(group_size), group_size)
  w <- as.numeric(rowsum(weights, group, reorder = FALSE)) / group_size
  if (any(abs(weights - w[group]) > 1e-10 * pmax(1, abs(weights)))) {
    stop("panfmm() weights must be constant within individuals ",
         "(one weight per individual, e.g. a sampling weight)", call. = FALSE)
  }
  w
}

#' Split the parameter vector into class blocks and the prior matrix
#'
#' @keywords internal
#' @noRd
fmm_split <- function(spec, theta) {
  k <- spec$k; m <- spec$m
  blocks <- lapply(seq_len(k), function(c) theta[(c - 1L) * m + seq_len(m)])
  G <- matrix(0, spec$q, k)
  if (k > 1L) G[, 2:k] <- theta[k * m + seq_len(spec$q * (k - 1L))]
  list(blocks = blocks, G = G)
}

#' Class probabilities (multinomial logit), one row per individual
#'
#' @keywords internal
#' @noRd
fmm_logprior <- function(Z, G) {
  eta <- Z %*% G
  mx <- fmm_rowmax(eta)
  eta - (mx + log(rowSums(exp(eta - mx))))
}

#' Row maxima of a matrix
#'
#' @keywords internal
#' @noRd
fmm_rowmax <- function(a) a[cbind(seq_len(nrow(a)), max.col(a, ties.method = "first"))]

#' Observation log-density and score of one class
#'
#' @keywords internal
#' @noRd
fmm_component <- function(spec, block, score = FALSE) {
  X <- spec$X; y <- spec$y; p <- spec$p
  beta <- block[seq_len(p)]
  eta <- as.numeric(X %*% beta) + spec$offset
  fam <- spec$family
  if (fam == "gaussian") {
    s <- exp(block[p + 1L])
    r <- y - eta
    ll <- stats::dnorm(y, eta, s, log = TRUE)
    if (score) sc <- cbind(X * (r / s^2), r^2 / s^2 - 1)
    mu <- eta
  } else if (fam == "poisson") {
    mu <- exp(eta)
    ll <- stats::dpois(y, mu, log = TRUE)
    if (score) sc <- X * (y - mu)
  } else if (fam == "binomial" && spec$probit) {
    mu <- stats::pnorm(eta)
    lp1 <- stats::pnorm(eta, log.p = TRUE)
    lp0 <- stats::pnorm(eta, lower.tail = FALSE, log.p = TRUE)
    ll <- y * lp1 + (1 - y) * lp0
    if (score) {
      ld <- stats::dnorm(eta, log = TRUE)
      sc <- X * (y * exp(ld - lp1) - (1 - y) * exp(ld - lp0))
    }
  } else if (fam == "binomial") {
    mu <- stats::plogis(eta)
    ll <- stats::dbinom(y, 1, mu, log = TRUE)
    if (score) sc <- X * (y - mu)
  } else {
    th <- exp(block[p + 1L])
    mu <- exp(eta)
    ll <- lgamma(y + th) - lgamma(th) - spec$lgy + th * log(th) + y * eta -
      (y + th) * log(th + mu)
    if (score) {
      sc <- cbind(X * ((y - mu) * th / (th + mu)),
                  th * (digamma(y + th) - digamma(th) + log(th) + 1 -
                          log(th + mu) - (y + th) / (th + mu)))
    }
  }
  list(ll = ll, mu = mu, score = if (score) sc else NULL)
}

#' Individual log-likelihood contributions, scores and posteriors
#'
#' @keywords internal
#' @noRd
fmm_contrib <- function(spec, theta, score = FALSE) {
  sp <- fmm_split(spec, theta)
  k <- spec$k
  comps <- lapply(sp$blocks, function(b) fmm_component(spec, b, score))
  lf <- vapply(comps, function(cm) rowsum(cm$ll, spec$group, reorder = FALSE)[, 1L],
               numeric(spec$N))
  lf <- matrix(lf, spec$N, k)
  lp <- fmm_logprior(spec$Z, sp$G)
  a <- lp + lf
  mx <- fmm_rowmax(a)
  ll <- mx + log(rowSums(exp(a - mx)))
  tau <- exp(a - ll)
  out <- list(ll = ll, tau = tau, logprior = lp)
  if (score) {
    s <- lapply(seq_len(k), function(c) {
      tau[, c] * rowsum(comps[[c]]$score, spec$group, reorder = FALSE)
    })
    if (k > 1L) {
      pr <- exp(lp)
      s <- c(s, lapply(2:k, function(c) spec$Z * (tau[, c] - pr[, c])))
    }
    out$score <- do.call(cbind, s)
    if (any(!is.finite(out$score))) stop("non-finite score", call. = FALSE)
  }
  if (!all(is.finite(ll))) stop("non-finite log-likelihood", call. = FALSE)
  out
}

#' Weighted M step for one class
#'
#' @keywords internal
#' @noRd
fmm_mstep_class <- function(spec, w, block) {
  X <- spec$X; y <- spec$y; p <- spec$p; off <- spec$offset
  fam <- spec$family
  start <- if (is.null(block)) NULL else block[seq_len(p)]
  if (fam == "gaussian") {
    fit <- stats::lm.wfit(X, y - off, w)
    beta <- fit$coefficients
    beta[is.na(beta)] <- 0
    r <- y - off - as.numeric(X %*% beta)
    return(c(beta, log(sqrt(max(sum(w * r^2) / sum(w), 1e-12)))))
  }
  if (fam == "negbin") {
    if (is.null(block) && !spec$has_offset) {
      res <- negbin_irls_fit_cpp(X, y, 100L, 1e-8, w)
      block <- c(as.numeric(res$coefficients), log(min(max(res$theta, 1e-4), 1e6)))
    } else if (is.null(block)) {
      fit <- suppressWarnings(stats::glm.fit(X, y, weights = w, offset = off,
                                             family = stats::poisson()))
      beta <- fit$coefficients
      beta[is.na(beta)] <- 0
      return(fmm_mstep_negbin(X, y, w, c(beta, 0), off, maxit = 100L))
    }
    return(fmm_mstep_negbin(X, y, w, block, off))
  }
  famobj <- if (fam == "poisson") {
    stats::poisson()
  } else {
    stats::binomial(link = if (spec$probit) "probit" else "logit")
  }
  if (!is.null(start)) {
    beta <- fmm_irls(X, y, w, off, famobj, start)
    if (!is.null(beta)) return(beta)
  }
  fit <- suppressWarnings(stats::glm.fit(X, y, weights = w, offset = off, family = famobj,
                                         start = start))
  beta <- fit$coefficients
  beta[is.na(beta)] <- 0
  beta
}

#' Warm-started weighted IRLS for the Poisson and binomial M steps
#'
#' The IRLS iterations of [stats::glm.fit()] from the previous EM
#' iteration's coefficients, with its convergence rule on the deviance
#' (relative change below `eps`) but without its per-call overhead
#' (deviance residual checks, AIC, QR with pivoting): the weighted normal
#' equations of the small `p x p` system are solved directly. Returns
#' `NULL` (the caller then falls back to `glm.fit()`) when the system is
#' singular or the deviance is not finite.
#'
#' @keywords internal
#' @noRd
fmm_irls <- function(X, y, w, offset, famobj, start, maxit = 25L, eps = 1e-8) {
  keep <- w > 0
  if (!all(keep)) {
    X <- X[keep, , drop = FALSE]; y <- y[keep]; w <- w[keep]; offset <- offset[keep]
  }
  beta <- start
  eta <- as.numeric(X %*% beta) + offset
  mu <- famobj$linkinv(eta)
  dev <- sum(famobj$dev.resids(y, mu, w))
  if (!is.finite(dev)) return(NULL)
  for (it in seq_len(maxit)) {
    me <- famobj$mu.eta(eta)
    z <- eta - offset + (y - mu) / me
    ww <- w * me^2 / famobj$variance(mu)
    if (any(!is.finite(z)) || any(!is.finite(ww))) return(NULL)
    bnew <- tryCatch(as.numeric(solve(crossprod(X, X * ww), crossprod(X, ww * z))),
                     error = function(e) NULL)
    if (is.null(bnew) || any(!is.finite(bnew))) return(NULL)
    eta <- as.numeric(X %*% bnew) + offset
    mu <- famobj$linkinv(eta)
    devnew <- sum(famobj$dev.resids(y, mu, w))
    if (!is.finite(devnew)) return(NULL)
    beta <- bnew
    if (abs(devnew - dev) / (abs(devnew) + 0.1) < eps) break
    dev <- devnew
  }
  beta
}

#' Warm-started weighted NB2 M step
#'
#' Alternates one IRLS step for the coefficients given `theta` with one
#' Newton step for `log(theta)` given the coefficients.
#'
#' @keywords internal
#' @noRd
fmm_mstep_negbin <- function(X, y, w, block, offset = 0, maxit = 5L, tol = 1e-8) {
  p <- ncol(X)
  beta <- block[seq_len(p)]
  phi <- block[p + 1L]
  for (it in seq_len(maxit)) {
    th <- exp(phi)
    eta <- as.numeric(X %*% beta) + offset
    mu <- exp(eta)
    ww <- w * mu * th / (th + mu)
    z <- eta - offset + (y - mu) / mu
    fit <- stats::lm.wfit(X, z, ww)
    bnew <- fit$coefficients
    if (anyNA(bnew)) break
    mu <- exp(as.numeric(X %*% bnew) + offset)
    s_th <- digamma(y + th) - digamma(th) + log(th) + 1 - log(th + mu) - (y + th) / (th + mu)
    h_th <- trigamma(y + th) - trigamma(th) + 1 / th - 2 / (th + mu) + (y + th) / (th + mu)^2
    g <- sum(w * th * s_th)
    h <- sum(w * (th * s_th + th^2 * h_th))
    step <- if (h < 0) -g / h else sign(g) * 0.5
    step <- max(min(step, 2), -2)
    phinew <- min(max(phi + step, -10), 15)
    change <- max(abs(bnew - beta), abs(phinew - phi))
    beta <- bnew
    phi <- phinew
    if (change < tol) break
  }
  c(beta, phi)
}

#' Weighted multinomial logit M step for the class probabilities
#'
#' Maximizes `sum_i w_i sum_c tau_ic log pi_ic(z_i)` by Newton steps.
#'
#' @keywords internal
#' @noRd
fmm_mstep_prior <- function(Z, tau, G, w) {
  k <- ncol(tau); q <- ncol(Z)
  if (q == 1L && all(Z == 1)) {
    pm <- pmax(colSums(w * tau) / sum(w), 1e-300)
    G[1L, ] <- log(pm) - log(pm[1L])
    return(G)
  }
  obj <- function(G) sum(w * tau * fmm_logprior(Z, G))
  for (it in seq_len(25L)) {
    pr <- exp(fmm_logprior(Z, G))
    g <- as.numeric(crossprod(Z, w * (tau[, -1L, drop = FALSE] - pr[, -1L, drop = FALSE])))
    H <- matrix(0, q * (k - 1L), q * (k - 1L))
    for (a in 2:k) for (b in 2:k) {
      wab <- w * pr[, a] * ((a == b) - pr[, b])
      ia <- (a - 2L) * q + seq_len(q); ib <- (b - 2L) * q + seq_len(q)
      H[ia, ib] <- -crossprod(Z * wab, Z)
    }
    step <- tryCatch(solve(-H, g), error = function(e) NULL)
    if (is.null(step)) break
    old <- obj(G)
    for (h in 0:20) {
      Gn <- G
      Gn[, 2:k] <- G[, 2:k] + matrix(step / 2^h, q)
      if (obj(Gn) >= old - 1e-12) break
    }
    G <- Gn
    if (max(abs(g)) < 1e-9) break
  }
  G
}

#' Parameter vector from class blocks and the prior matrix
#'
#' @keywords internal
#' @noRd
fmm_pack <- function(spec, blocks, G) {
  c(unlist(blocks), if (spec$k > 1L) as.numeric(G[, -1L, drop = FALSE]))
}

#' One EM run from an initial posterior matrix
#'
#' @keywords internal
#' @noRd
fmm_em <- function(spec, tau, maxit, tol) {
  k <- spec$k
  G <- matrix(0, spec$q, k)
  blocks <- vector("list", k)
  old <- -Inf
  ll <- -Inf
  for (it in seq_len(maxit)) {
    blocks <- lapply(seq_len(k), function(c) {
      fmm_mstep_class(spec, (spec$w * tau[, c])[spec$group], blocks[[c]])
    })
    if (k > 1L) G <- fmm_mstep_prior(spec$Z, tau, G, spec$w)
    par <- fmm_pack(spec, blocks, G)
    st <- fmm_contrib(spec, par, FALSE)
    ll <- sum(spec$w * st$ll)
    tau <- st$tau
    if (k > 1L && any(colSums(spec$w * tau) < 1e-6 * sum(spec$w))) {
      stop("an empty class", call. = FALSE)
    }
    if (k == 1L || abs(ll - old) <= tol * (abs(ll) + tol)) break
    old <- ll
  }
  list(par = par, loglik = ll, iterations = it)
}

#' EM from several starts; the first splits individuals by mean outcome
#'
#' @keywords internal
#' @noRd
fmm_em_starts <- function(spec, nstart, seed, maxit, tol) {
  k <- spec$k; N <- spec$N
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else NULL
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = globalenv())
      } else {
        assign(".Random.seed", old_seed, envir = globalenv())
      }
    }, add = TRUE)
    set.seed(seed)
  }
  hard <- function(cl) {
    tau <- matrix(0.05 / max(k - 1L, 1L), N, k)
    tau[cbind(seq_len(N), cl)] <- 0.95
    if (k == 1L) tau[] <- 1
    tau
  }
  ybar <- as.numeric(rowsum(spec$y, spec$group, reorder = FALSE)) / tabulate(spec$group)
  first <- if (k == 1L) rep(1L, N) else {
    r <- rank(ybar, ties.method = "first")
    as.integer(ceiling(r * k / N))
  }
  n_runs <- if (k == 1L) 1L else as.integer(nstart)
  runs <- vector("list", n_runs)
  for (s in seq_len(n_runs)) {
    cl <- if (s == 1L) first else c(seq_len(k), sample.int(k, N - k, replace = TRUE))[sample.int(N)]
    runs[[s]] <- tryCatch(fmm_em(spec, hard(cl), maxit, tol), error = function(e) NULL)
  }
  ll <- vapply(runs, function(r) if (is.null(r)) NA_real_ else r$loglik, numeric(1))
  if (all(is.na(ll))) stop("the EM algorithm failed from every start", call. = FALSE)
  best <- which.max(ll)
  list(par = runs[[best]]$par, starts = ll, iterations = runs[[best]]$iterations)
}

#' Order of the classes by ascending mean outcome
#'
#' @keywords internal
#' @noRd
fmm_class_means <- function(spec, theta) {
  sp <- fmm_split(spec, theta)
  vapply(sp$blocks, function(b) mean(fmm_component(spec, b)$mu), numeric(1))
}

fmm_class_order <- function(spec, theta) order(fmm_class_means(spec, theta))

#' Relabel the classes by ascending mean outcome, class 1 the reference
#'
#' @keywords internal
#' @noRd
fmm_reorder <- function(spec, theta) {
  o <- fmm_class_order(spec, theta)
  sp <- fmm_split(spec, theta)
  G <- sp$G[, o, drop = FALSE]
  G <- G - G[, 1L]
  fmm_pack(spec, sp$blocks[o], G)
}

#' Readable pieces of the parameter vector
#'
#' @keywords internal
#' @noRd
fmm_unpack <- function(spec, theta) {
  sp <- fmm_split(spec, theta)
  k <- spec$k; p <- spec$p
  cc <- vapply(sp$blocks, function(b) b[seq_len(p)], numeric(p))
  cc <- matrix(cc, p, k, dimnames = list(colnames(spec$X), paste0("class", seq_len(k))))
  out <- list(classcoef = cc)
  if (spec$aux) {
    a <- stats::setNames(exp(vapply(sp$blocks, function(b) b[p + 1L], numeric(1))),
                         colnames(cc))
    if (spec$family == "gaussian") out$sigma <- a else out$theta <- a
  }
  pr <- exp(fmm_logprior(spec$Z, sp$G))
  colnames(pr) <- colnames(cc)
  out$prior <- pr
  out$prior_mean <- colSums(spec$w * pr) / sum(spec$w)
  out$prior_coef <- matrix(sp$G, spec$q, k, dimnames = list(colnames(spec$Z), colnames(cc)))
  out$classmean <- stats::setNames(fmm_class_means(spec, theta), colnames(cc))
  out
}

#' Covariance of a panfmm fit
#'
#' @param object a `"panfmm"` fit
#' @param type `NULL` (the covariance stored on the fit), `"classical"`
#'   (inverse of the negative Hessian of the observed-data log-likelihood)
#'   or `"cluster"` (sandwich with individual scores summed within
#'   clusters, factor `G / (G - 1)`; with weights the individual scores are
#'   multiplied by the individual weights)
#' @param cluster optional clustering vector with one value per observation
#'   used in the fit (in the row order of the data), constant within
#'   individuals; defaults to the individual
#' @param ... unused
#' @return the covariance matrix of all parameters
#' @family panfmm methods
#' @examples
#' set.seed(3)
#' d <- data.frame(id = rep(1:150, each = 3), t = rep(1:3, 150), x = rnorm(450))
#' d$y <- rpois(450, exp(ifelse(rep(rbinom(150, 1, 0.4), each = 3) == 1, 1.5, 0) + 0.3 * d$x))
#' fit <- panfmm(y ~ x, data = d, index = c("id", "t"), nstart = 2, seed = 1)
#' sqrt(diag(vcov(fit, type = "cluster")))
#' @export
vcov.panfmm <- function(object, type = NULL, cluster = NULL, ...) {
  if (is.null(type)) return(object$vcov)
  type <- match.arg(type, c("classical", "cluster"))
  if (type == "classical") return(object$vcov_classical)
  cl <- if (is.null(cluster)) {
    seq_len(object$n_groups)
  } else {
    individual_clusters(panglm_cluster_vector(object, cluster), object$group_size)
  }
  v <- ml_sandwich(object$vcov_classical, object$score,
                   object$weights %||% rep(1, object$n_groups),
                   cluster = cl, clusters_only = TRUE)
  dimnames(v) <- dimnames(object$vcov_classical)
  v
}

#' @export
coef.panfmm <- function(object, ...) object$coefficients

#' @export
nobs.panfmm <- function(object, ...) object$nobs

#' @export
logLik.panfmm <- function(object, ...) {
  structure(object$loglik, df = object$df, nobs = object$nobs, class = "logLik")
}

#' @export
confint.panfmm <- function(object, parm, level = 0.95, ...) {
  confint.panglm(object, parm, level = level, ...)
}

#' Predictions from a latent class panel model
#'
#' @param object a `"panfmm"` fit
#' @param newdata optional data.frame; by default the estimation sample (in
#'   its original row order). Offsets are included in the class means:
#'   the estimation-sample offsets without `newdata`, the `offset()` terms
#'   of the formula evaluated in `newdata` otherwise (an `offset` argument
#'   of [panfmm()] cannot be evaluated in `newdata` and is not applied,
#'   with a warning)
#' @param type `"response"` for the marginal mean
#'   `sum_c w_ic mu_c(x_it)`; `"class"` for the modal class of each
#'   observation's individual (an integer, the class with the largest
#'   posterior probability in the estimation sample, or the largest prior
#'   probability, from the concomitant variables when there are some, with
#'   `newdata`); `"posterior"` / `"prior"` for the matrix of class
#'   probabilities of each observation's individual; `"components"` for the
#'   matrix of class-specific means `mu_c(x_it)` (one column per class).
#'   In panglm 1.5.0 and earlier, `type = "class"` returned the
#'   class-specific means now given by `"components"`; a message recalls
#'   this once per session.
#' @param weights for `type = "response"`: `"prior"` (the class
#'   probabilities, from the concomitant variables when there are some) or
#'   `"posterior"` (the in-sample posterior probabilities; only without
#'   `newdata`)
#' @param ... unused
#' @return a numeric vector (`"response"`), an integer vector of class
#'   labels `1, ..., k` (`"class"`), or a matrix with one column per class
#'   (`"posterior"`, `"prior"`, `"components"`)
#' @family panfmm methods
#' @examples
#' set.seed(3)
#' d <- data.frame(id = rep(1:150, each = 3), t = rep(1:3, 150), x = rnorm(450))
#' d$y <- rpois(450, exp(ifelse(rep(rbinom(150, 1, 0.4), each = 3) == 1, 1.5, 0) + 0.3 * d$x))
#' fit <- panfmm(y ~ x, data = d, index = c("id", "t"), nstart = 2, seed = 1)
#' head(predict(fit))
#' table(predict(fit, type = "class"))              # modal class per row
#' head(predict(fit, type = "posterior"))           # class probabilities
#' predict(fit, newdata = data.frame(x = c(-1, 0, 1)), type = "components")
#' @export
predict.panfmm <- function(object, newdata = NULL,
                           type = c("response", "class", "prior", "posterior",
                                    "components"),
                           weights = c("prior", "posterior"), ...) {
  type <- match.arg(type)
  if (type == "class") fmm_class_notice()
  weights <- match.arg(weights)
  k <- object$k
  if (is.null(newdata)) {
    X <- object$X
    off <- object$offset %||% 0
    prior <- object$prior[rep(seq_len(object$n_groups), object$group_size), , drop = FALSE]
    post <- object$posterior[rep(seq_len(object$n_groups), object$group_size), , drop = FALSE]
    back <- order(object$panel_order)
  } else {
    if (type == "posterior" || (type == "response" && weights == "posterior")) {
      stop("posterior probabilities are only available for the estimation sample",
           call. = FALSE)
    }
    tt <- stats::delete.response(object$terms)
    mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass,
                             xlev = object$xlevels)
    X <- stats::model.matrix(tt, mf)
    off <- stats::model.offset(mf) %||% 0
    if (isTRUE(object$offset_argument)) {
      warning("the 'offset' argument of panfmm() is not applied to 'newdata'; ",
              "put the offset in the formula as offset()", call. = FALSE)
    }
    tz <- object$terms_concomitant
    mz <- stats::model.frame(tz, newdata, na.action = stats::na.pass, xlev = object$xlevels)
    Zn <- stats::model.matrix(tz, mz)
    if (nrow(Zn) != nrow(X)) Zn <- matrix(1, nrow(X), 1L)
    prior <- exp(fmm_logprior(Zn, object$prior_coef))
    post <- NULL
    back <- seq_len(nrow(X))
  }
  cn <- paste0("class", seq_len(k))
  eta <- X %*% object$classcoef + off
  probit <- identical(object$family$link, "probit")
  mu <- switch(object$family$family,
               gaussian = eta,
               binomial = if (probit) stats::pnorm(eta) else stats::plogis(eta),
               exp(eta))
  mu <- matrix(mu, nrow(X), k, dimnames = list(NULL, cn))
  if (type == "class") {
    probs <- if (is.null(post)) prior else post
    return(as.integer(max.col(probs, ties.method = "first"))[back])
  }
  out <- switch(type,
    components = mu,
    prior = prior,
    posterior = post,
    response = rowSums(mu * if (weights == "prior") prior else post)
  )
  if (is.matrix(out)) {
    out <- out[back, , drop = FALSE]
    dimnames(out) <- list(NULL, cn)
    out
  } else {
    as.numeric(out[back])
  }
}

# Tell users of the old predict(type = "class") (class-specific means)
# once per session that it now returns class labels.
fmm_class_notice <- function() {
  key <- "panfmm predict class"
  if (isTRUE(.panglm_deprecated[[key]])) return(invisible(FALSE))
  assign(key, TRUE, envir = .panglm_deprecated)
  msg <- paste0("predict.panfmm(type = \"class\") now returns the modal class label; ",
                "the class-specific means it returned in panglm 1.5.0 and earlier are ",
                "type = \"components\". This message is shown once per session.\n")
  message(structure(list(message = msg, call = NULL),
                    class = c("panglm_behaviour_change", "message", "condition")))
  invisible(TRUE)
}

#' @export
print.panfmm <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("panfmm(", fmm_label(x), ")\n\n", sep = "")
  cat("Class-specific coefficients:\n")
  print(x$classcoef, digits = digits)
  fmm_print_aux(x, digits)
  cat("\nClass probabilities:\n")
  print(x$prior_mean, digits = digits)
  cat("\nlogLik =", format(x$loglik, digits = digits + 3), " AIC =",
      format(stats::AIC(x), digits = digits + 3), " BIC =",
      format(stats::BIC(x), digits = digits + 3), "\n")
  cat("N =", x$nobs, " groups =", x$n_groups, "\n")
  invisible(x)
}

fmm_label <- function(x) {
  fam <- switch(x$family$family, negbin = "NB2", poisson = "Poisson",
                gaussian = "Gaussian",
                binomial = if (identical(x$family$link, "probit")) "probit" else "logit")
  paste0(x$k, "-class latent class ", fam, " model")
}

fmm_print_aux <- function(x, digits) {
  if (!is.null(x$sigma)) {
    cat("\nsigma:\n"); print(x$sigma, digits = digits)
  }
  if (!is.null(x$theta)) {
    cat("\ntheta (NB2, variance mu + mu^2 / theta):\n"); print(x$theta, digits = digits)
  }
}

#' Summarize a panfmm fit
#'
#' @param object a `"panfmm"` fit
#' @param ... unused
#' @return an object of class `"summary.panfmm"`: a list with one
#'   coefficient table per class (`classes`: estimate, standard error, z
#'   value and p-value, ancillary parameter on the log scale), the table of
#'   the class-probability (multinomial logit) parameters (`prior`), the mean
#'   class probabilities, the class sizes by modal posterior class, and the
#'   fit information (`call`, `loglik`, `AIC`, `BIC`, `nobs`, `n_groups`,
#'   `vcov_type`).
#' @family panfmm methods
#' @examples
#' set.seed(3)
#' d <- data.frame(id = rep(1:150, each = 3), t = rep(1:3, 150), x = rnorm(450))
#' d$y <- rpois(450, exp(ifelse(rep(rbinom(150, 1, 0.4), each = 3) == 1, 1.5, 0) + 0.3 * d$x))
#' summary(panfmm(y ~ x, data = d, index = c("id", "t"), nstart = 2, seed = 1))
#' @export
summary.panfmm <- function(object, ...) {
  est <- object$coefficients
  se <- sqrt(diag(object$vcov))
  z <- est / se
  tab <- cbind(Estimate = est, `Std. Error` = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  nm <- names(est)
  classes <- lapply(seq_len(object$k), function(c) {
    pre <- paste0("class", c, ":")
    t <- tab[startsWith(nm, pre), , drop = FALSE]
    rownames(t) <- substring(rownames(t), nchar(pre) + 1L)
    t
  })
  prior <- tab[startsWith(nm, "prior"), , drop = FALSE]
  out <- list(call = object$call, label = fmm_label(object), classes = classes,
              prior = prior, prior_mean = object$prior_mean,
              size = tabulate(object$class, object$k), classmean = object$classmean,
              sigma = object$sigma, theta = object$theta, loglik = object$loglik,
              AIC = stats::AIC(object), BIC = stats::BIC(object), df = object$df,
              nobs = object$nobs, n_groups = object$n_groups, vcov_type = object$vcov_type)
  class(out) <- "summary.panfmm"
  out
}

#' @export
print.summary.panfmm <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nModel:", x$label, "\n")
  cat("N =", x$nobs, " groups =", x$n_groups, " vcov:", x$vcov_type, "\n")
  for (c in seq_along(x$classes)) {
    cat(sprintf("\nClass %d (prior %.3f, %d individuals by modal posterior, mean outcome %s):\n",
                c, x$prior_mean[c], x$size[c], format(x$classmean[c], digits = digits)))
    stats::printCoefmat(x$classes[[c]], digits = digits, signif.stars = TRUE,
                        has.Pvalue = TRUE)
  }
  fmm_print_aux(x, digits)
  if (nrow(x$prior)) {
    cat("\nClass probabilities (multinomial logit, class 1 = reference):\n")
    stats::printCoefmat(x$prior, digits = digits, signif.stars = TRUE, has.Pvalue = TRUE)
  }
  cat("\nlogLik =", format(x$loglik, digits = digits + 3), " df =", x$df,
      " AIC =", format(x$AIC, digits = digits + 3),
      " BIC =", format(x$BIC, digits = digits + 3), "\n\n")
  invisible(x)
}

#' Compare latent class panel models with different numbers of classes
#'
#' Fits [panfmm()] for each number of classes in `k` and tabulates the
#' log-likelihood, number of parameters, AIC and BIC (with `n` the number of
#' observations), the usual criteria for choosing the number of classes
#' (Deb and Trivedi 2002 prefer BIC). Mixtures with too many classes may
#' fail to converge or return near-empty classes; such fits are reported
#' with `NA`.
#'
#' @param formula,data,index,family,concomitant,... passed to [panfmm()]
#' @param k integer vector of numbers of classes
#' @return a data.frame with columns `k`, `logLik`, `df`, `AIC`, `BIC` and
#'   `smallest` (the smallest mean class probability), with the fits in the
#'   attribute `"fits"`
#' @family panel models
#' @section Validation:
#' Each fit is a [panfmm()] fit, whose log-likelihood and BIC are checked
#' against flexmix in the package's validation suite on GitHub.
#' @examples
#' set.seed(3)
#' d <- data.frame(id = rep(1:150, each = 3), t = rep(1:3, 150), x = rnorm(450))
#' d$y <- rpois(450, exp(ifelse(rep(rbinom(150, 1, 0.4), each = 3) == 1, 1.5, 0) + 0.3 * d$x))
#' panfmmselect(y ~ x, data = d, index = c("id", "t"), k = 1:2, nstart = 2, seed = 1)
#' \donttest{
#' # up to three classes (a few seconds)
#' panfmmselect(y ~ x, data = d, index = c("id", "t"), k = 1:3, nstart = 2, seed = 1)
#' }
#' @export
panfmmselect <- function(formula, data, index, k = 1:3, family = "poisson",
                         concomitant = NULL, ...) {
  fits <- lapply(k, function(kk) {
    tryCatch(panfmm(formula, data, index, k = kk, family = family,
                    concomitant = concomitant, ...),
             error = function(e) NULL)
  })
  row <- function(f) {
    if (is.null(f)) return(c(NA, NA, NA, NA, NA))
    c(f$loglik, f$df, stats::AIC(f), stats::BIC(f), min(f$prior_mean))
  }
  tab <- do.call(rbind, lapply(fits, row))
  out <- data.frame(k = k, logLik = tab[, 1], df = tab[, 2], AIC = tab[, 3],
                    BIC = tab[, 4], smallest = tab[, 5])
  names(fits) <- paste0("k", k)
  attr(out, "fits") <- fits
  out
}
