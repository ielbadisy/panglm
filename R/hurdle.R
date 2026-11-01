#' Fixed-effects hurdle model for panel count data
#'
#' A hurdle decomposition for panels where a large fraction of groups report
#' *structural* zeros (never observed to occur), rather than zeros arising
#' from ordinary Poisson/NB sampling variation at a low mean, a pattern
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
#'     "binomial")` (Chamberlain 1980): no new estimator, just a new
#'     outcome vector.
#'   \item **Count given `y > 0`.** A zero-truncated Poisson or NB2
#'     fixed-effects model (Allison-Waterman-style dummy-variable
#'     intercepts, the same unconditional approach used for FE-NB2) on the
#'     subsample where `y > 0`. The log-density is the ordinary count
#'     log-density minus `log(1 - P(Y = 0))`; see `pscl::hurdle()` for the
#'     reference form of the truncated-count likelihood.
#' }
#'
#' Groups with no remaining within-group variation after truncation (a
#' single surviving observation) contribute no information about the shared
#' count-part slope, the same profiling fact that applies to the other
#' fixed-effects estimators in this package; they are not dropped
#' explicitly since a perfectly-fit singleton doesn't diverge, it simply
#' carries zero weight. Groups whose positive counts are all equal to 1 are
#' dropped from the count part with a message: their intercept has no
#' finite estimate (the truncated likelihood approaches its maximum as the
#' group mean goes to 0, whatever the slopes), so they carry no
#' information either.
#'
#' Both parts remove one effect per individual, so regressors that do not
#' vary within individuals (sex, baseline covariates) are not identified:
#' they are dropped with a warning, as in `panglm(model = "within")`. A
#' regressor that varies only in individuals a part does not use (for
#' example individuals with a single positive count, for the count part) is
#' dropped from that part alone, with a warning.
#'
#' The count part is fit by Newton-Raphson with the analytic gradient and
#' Hessian of the joint likelihood. The group intercepts enter through a
#' diagonal Hessian block that is eliminated at each step, so the cost is
#' linear in the number of groups and the intercepts are never stored as
#' dummy columns. Standard errors are the slope block of the inverse
#' observed information of the joint likelihood (slopes, intercepts and,
#' for NB2, `log(theta)`). When the NB2 shape parameter runs off to
#' infinity (no overdispersion), it is stopped at `1e8` with a warning and
#' the fit equals the zero-truncated Poisson one.
#'
#' @param formula a model formula, e.g. `y ~ x1 + x2`
#' @param data a data.frame containing the variables in `formula` and `index`
#' @param index length-2 character vector identifying the individual and
#'   time columns, e.g. `c("id", "time")`
#' @param family count distribution for positive observations, either
#'   `"poisson"` or `"negbin"` (NB2 parameterization). It is stored on the
#'   fit as `count_family`.
#' @param maxit maximum IRLS/Newton iterations
#' @param tol convergence tolerance
#' @return an object of class `"panglm_hurdle"`: a list with `zero` (the
#'   [panglm()] binomial fit for the zero-vs-positive part) and `count` (the
#'   zero-truncated fixed-effects fit for the count part). Methods: `print`,
#'   `summary`, `coef`, `vcov`, `confint`, `nobs`, `plot`, `tidy` and `glance`.
#' @family panel models
#' @seealso [hurdle-methods] and [plot.panglm_hurdle()] for the fitted model,
#'   [pandispersiontest()] to check the count part, [twopart()] for continuous
#'   costs.
#' @section Validation:
#' Numerical accuracy is checked in the package tests and in the validation
#' suite; see the validation article on the package website.
#' @examples
#' data(copd)
#' fit <- panhurdle(exacerbations ~ crp, data = copd, index = c("id", "visit"))
#' fit
#' @export
panhurdle <- function(formula, data, index,
                      family = c("poisson", "negbin"),
                      maxit = 100, tol = 1e-10) {
  count_family <- match.arg(family)
  if (missing(index) || length(index) != 2) {
    stop("'index' must be a two-element character vector, e.g. index = c(\"id\", \"time\")", call. = FALSE)
  }

  stop_panel_calls(formula, what = "panhurdle()")
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
  group <- rep(seq_along(group_size), group_size)

  # Both parts condition on (or estimate) one effect per individual, so a
  # regressor constant within every individual (sex, a baseline covariate)
  # is not identified: it is dropped with a warning, as in
  # panglm(model = "within"). Each part is then checked on the individuals
  # it actually uses.
  if (!ncol(X)) stop("the model contains no regressors", call. = FALSE)
  varying <- hurdle_identified(X, group, rep(TRUE, length(y)))
  invariant <- setdiff(colnames(X), varying)
  if (length(invariant)) {
    warning("time-invariant regressor(s) dropped (not identified with fixed effects): ",
            paste(invariant, collapse = ", "), call. = FALSE)
    if (!length(varying)) stop("no time-varying regressors left", call. = FALSE)
    X <- X[, varying, drop = FALSE]
  }

  y_bin <- as.numeric(y > 0)
  informative <- rowsum(y_bin, group, reorder = TRUE)[, 1]
  informative <- (informative > 0 & informative < group_size)[group]
  X_zero <- hurdle_part_columns(X, group, informative, "zero")
  zero_fit <- fit_within_binomial(X_zero, y_bin, group_start, group_size, maxit, tol)
  zero_fit$formula <- formula; zero_fit$index <- index
  class(zero_fit) <- "panglm"

  keep <- y > 0
  not_ones <- rowsum(as.numeric(keep & y != 1), group, reorder = TRUE)[, 1] > 0
  X_count <- hurdle_part_columns(X, group, keep & not_ones[group], "count")
  count_fit <- fit_truncated_count_fe(
    X_count[keep, , drop = FALSE], y[keep], group[keep], count_family, maxit, tol
  )

  out <- list(zero = zero_fit, count = count_fit, formula = formula, index = index,
              count_family = count_family, n = length(y), n_positive = sum(keep),
              group_levels = panel$group_id)
  class(out) <- "panglm_hurdle"
  out
}

# Columns of X identified once each individual's effect is removed,
# judged on the rows in `rows` (the individuals a part actually uses).
hurdle_identified <- function(X, group, rows) {
  if (!any(rows)) return(character(0))
  g <- group[rows]
  size <- tabulate(match(g, unique(g)))
  start <- as.integer(c(0L, cumsum(size)[-length(size)]))
  cols <- tryCatch(
    panglm_estimable_columns(X[rows, , drop = FALSE], "within", "individual", start, size, NULL),
    error = function(e) integer(0)
  )
  colnames(X)[cols]
}

# The regressors of one hurdle part: those not identified on its rows are
# dropped with a warning.
hurdle_part_columns <- function(X, group, rows, part) {
  ok <- hurdle_identified(X, group, rows)
  if (!length(ok)) {
    stop("no regressor is identified in the ", part, " part", call. = FALSE)
  }
  if (length(ok) < ncol(X)) {
    warning("regressor(s) not identified in the ", part, " part dropped from it: ",
            paste(setdiff(colnames(X), ok), collapse = ", "), call. = FALSE)
  }
  X[, ok, drop = FALSE]
}

#' Zero-truncated Poisson or NB2 count part with group intercepts
#'
#' Joint maximum likelihood for the covariate slopes, one intercept per
#' group (Allison-Waterman-style dummy-variable estimator) and, for NB2,
#' `s = log(theta)`, by Newton-Raphson with the analytic gradient and
#' observed Hessian. The group intercepts are never formed as dummy
#' columns: their Hessian block is diagonal, so each Newton step solves
#' only the small Schur complement for the slopes (and `s`), at a cost
#' linear in the number of observations and groups.
#'
#' Groups whose positive counts are all equal to 1 have no finite
#' intercept (the truncated likelihood tends to its supremum 1 as the
#' group mean goes to 0, whatever the slopes), so they carry no
#' information and are set aside with a message; their fitted value is
#' the limit 1.
#'
#' The returned `vcov` is the slope block of the inverse observed
#' information of the joint likelihood (slopes, intercepts and, for NB2,
#' `log(theta)`).
#'
#' @keywords internal
#' @noRd
fit_truncated_count_fe <- function(X, y, group, family = c("poisson", "negbin"),
                                   maxit = 100, tol = 1e-10, quiet = FALSE) {
  family <- match.arg(family)
  colnames_X <- colnames(X)
  X <- as.matrix(X)
  n_all <- length(y)
  kept_groups <- sort(unique(group))
  g_all <- match(group, kept_groups)

  ones <- as.logical(rowsum(as.numeric(y != 1), g_all, reorder = TRUE) == 0)
  if (all(ones)) {
    stop("every group's positive counts are all 1: the count part is not identified",
         call. = FALSE)
  }
  if (any(ones)) {
    if (!quiet) {
      message(sum(ones), " group(s) whose positive counts are all 1 dropped from the ",
              "count part -- their intercept is not finite and they carry no information.")
    }
  }
  use <- !ones[g_all]
  fit_groups <- kept_groups[!ones]
  X <- X[use, , drop = FALSE]
  y <- y[use]
  group <- match(group[use], fit_groups)
  n <- length(y); k <- ncol(X); G <- length(fit_groups)
  nb <- family == "negbin"
  m <- k + nb

  terms_at <- function(beta, alpha, s, derivs = TRUE) {
    eta <- as.numeric(X %*% beta) + alpha[group]
    truncated_count_terms(eta, y, s, family, derivs)
  }
  loglik_at <- function(beta, alpha, s) {
    v <- sum(terms_at(beta, alpha, s, derivs = FALSE)$ll)
    if (is.finite(v)) v else -Inf
  }
  # gradient, Schur pieces of the Hessian over (beta, s) | alpha
  newton_parts <- function(beta, alpha, s) {
    tt <- terms_at(beta, alpha, s)
    gg <- as.numeric(crossprod(X, tt$e))
    Hgg <- crossprod(X, X * tt$ee)
    R <- rowsum(X * tt$ee, group, reorder = TRUE)        # G x k
    if (nb) {
      gg <- c(gg, sum(tt$s))
      hbs <- as.numeric(crossprod(X, tt$es))
      Hgg <- rbind(cbind(Hgg, hbs), c(hbs, sum(tt$ss)))
      R <- cbind(R, rowsum(tt$es, group, reorder = TRUE))
    }
    list(ll = sum(tt$ll), gg = gg, ga = as.numeric(rowsum(tt$e, group, reorder = TRUE)),
         Hgg = Hgg, R = R, D = as.numeric(rowsum(tt$ee, group, reorder = TRUE)))
  }
  newton_step <- function(p, lambda = 0) {
    # Marquardt damping: the diagonal is scaled by (1 + lambda)
    D <- p$D - lambda * abs(p$D)
    Hgg <- p$Hgg - diag(lambda * abs(diag(p$Hgg)), m)
    S <- Hgg - crossprod(p$R, p$R / D)
    rhs <- p$gg - as.numeric(crossprod(p$R, p$ga / D))
    dg <- tryCatch(-solve(S, rhs), error = function(e) NULL)
    if (is.null(dg) || any(!is.finite(dg))) return(NULL)
    da <- -(p$ga + as.numeric(p$R %*% dg)) / D
    list(dg = dg, da = da)
  }

  # intercepts alone at fixed slopes and s: the likelihood separates over
  # groups, so each group takes its own Newton step, halved where it fails
  update_alpha <- function(beta, alpha, s, steps = 25L) {
    xb <- as.numeric(X %*% beta)
    group_ll <- function(a) {
      v <- rowsum(truncated_count_terms(xb + a[group], y, s, family, FALSE)$ll,
                  group, reorder = TRUE)[, 1]
      v[!is.finite(v)] <- -Inf
      v
    }
    for (i in seq_len(steps)) {
      tt <- truncated_count_terms(xb + alpha[group], y, s, family)
      ga <- rowsum(tt$e, group, reorder = TRUE)[, 1]
      D <- rowsum(tt$ee, group, reorder = TRUE)[, 1]
      if (max(abs(ga)) < 1e-8) break
      da <- ga / pmax(-D, 1e-8 * (1 + abs(D)))
      ll0 <- group_ll(alpha)
      h <- rep(1, G)
      for (j in 1:20) {
        cand <- alpha + h * da
        bad <- group_ll(cand) < ll0 - 1e-12 * (abs(ll0) + 1)
        if (!any(bad)) break
        h[bad] <- h[bad] / 2
      }
      alpha <- ifelse(bad, alpha, cand)
    }
    alpha
  }

  # start: slopes 0, intercepts at log group means (the truncated-Poisson
  # moment start); for NB2 the truncated-Poisson fit, then log(theta) from
  # a one-dimensional search with the intercepts re-fit at each value
  beta <- rep(0, k)
  alpha <- log(pmax(as.numeric(rowsum(y, group, reorder = TRUE)) / tabulate(group, G) - 0.5, 0.1))
  s <- NULL
  s_max <- log(1e8)
  if (nb) {
    pf <- fit_truncated_count_fe(X, y, group, "poisson", maxit, tol, quiet = TRUE)
    beta <- unname(pf$coefficients)
    alpha <- unname(pf$individual_effects)
    profile_s <- function(sv) loglik_at(beta, update_alpha(beta, alpha, sv, 10L), sv)
    s <- stats::optimize(profile_s, c(log(1e-3), s_max), maximum = TRUE, tol = 0.05)$maximum
    alpha <- update_alpha(beta, alpha, s)
  }

  p <- newton_parts(beta, alpha, s)
  convergence <- 1L
  iter <- 0L
  lambda <- 0
  for (iter in seq_len(maxit)) {
    accepted <- FALSE
    for (attempt in 1:30) {
      st <- newton_step(p, lambda)
      if (!is.null(st)) {
        h <- 1
        repeat {
          b_new <- beta + h * st$dg[seq_len(k)]
          a_new <- alpha + h * st$da
          s_new <- if (nb) min(s + h * st$dg[m], s_max) else NULL
          ll_new <- loglik_at(b_new, a_new, s_new)
          if (ll_new >= p$ll - 1e-12 * (abs(p$ll) + 1)) { accepted <- TRUE; break }
          h <- h / 2
          if (h < 1e-4) break
        }
      }
      if (accepted) break
      lambda <- if (lambda == 0) 1e-3 else 10 * lambda
    }
    if (!accepted) break
    size <- max(abs(c(b_new - beta, a_new - alpha, if (nb) s_new - s)))
    scale <- 1 + max(abs(c(b_new, if (nb) s_new)))
    beta <- b_new; alpha <- a_new; s <- s_new
    ll_old <- p$ll
    p <- newton_parts(beta, alpha, s)
    if (size <= 1e-10 * scale ||
        (lambda == 0 && h == 1 && abs(p$ll - ll_old) <= tol * (abs(ll_old) + 1) &&
         size <= 1e-6 * scale)) {
      convergence <- 0L
      break
    }
    lambda <- if (lambda < 1e-6) 0 else lambda / 10
  }
  if (convergence != 0L) {
    warning("zero-truncated ", if (nb) "NB2" else "Poisson",
            " count fit did not converge", call. = FALSE)
  }
  theta <- if (nb) exp(s) else NULL
  if (nb && s >= s_max - 1e-8) {
    warning("NB2 theta reached its upper bound (", format(theta),
            "): no overdispersion beyond Poisson in the count part", call. = FALSE)
  }

  S <- p$Hgg - crossprod(p$R, p$R / p$D)
  vfull <- tryCatch(solve(-S), error = function(e) matrix(NA_real_, m, m))
  vcov <- vfull[seq_len(k), seq_len(k), drop = FALSE]
  dimnames(vcov) <- list(colnames_X, colnames_X)
  coefs <- stats::setNames(beta, colnames_X)

  eta <- as.numeric(X %*% beta) + alpha[group]
  mu <- exp(eta)
  p0 <- if (nb) exp(-exp(s) * log1p(mu / exp(s))) else exp(-mu)
  fitted <- rep(1, n_all)
  fitted[use] <- mu / (1 - p0)

  out <- list(coefficients = coefs, vcov = vcov, bread = vcov,
              fitted.values = fitted, loglik = p$ll, dispersion = 1,
              individual_effects = stats::setNames(alpha, fit_groups),
              n_used_groups = G, n_dropped_groups = sum(ones),
              df.residual = n - k - G - nb, iterations = iter,
              convergence = convergence)
  if (nb) out$theta <- theta
  out
}

fit_truncated_poisson_dummy <- function(X, y, group, maxit, tol) {
  fit_truncated_count_fe(X, y, group, "poisson", maxit, tol)
}

fit_truncated_negbin_dummy <- function(X, y, group, maxit, tol) {
  fit_truncated_count_fe(X, y, group, "negbin", maxit, tol)
}

# Per-observation log-likelihood and derivatives of the zero-truncated count
# density with respect to the linear predictor eta and (NB2) s = log(theta).
truncated_count_terms <- function(eta, y, s = NULL, family = "poisson", derivs = TRUE) {
  mu <- exp(eta)
  if (family == "poisson") {
    l1 <- log(-expm1(-mu))
    ll <- y * eta - mu - lgamma(y + 1) - l1
    if (!derivs) return(list(ll = ll))
    mu_t <- mu / (-expm1(-mu))
    return(list(ll = ll, e = y - mu_t, ee = -mu_t * (1 + mu - mu_t)))
  }
  th <- exp(s)
  a <- th + mu
  la <- log1p(mu / th)             # log(a / theta)
  q <- -th * la                    # log P(Y = 0)
  l1 <- log(-expm1(q))
  lf <- lgamma(y + th) - lgamma(th) - lgamma(y + 1) - th * la + y * (eta - log(a))
  ll <- lf - l1
  if (!derivs) return(list(ll = ll))
  r <- exp(q) / (-expm1(q))        # P0 / (1 - P0)
  q_e <- -th * mu / a
  q_ee <- -th^2 * mu / a^2
  q_t <- -la + mu / a
  q_tt <- mu^2 / (th * a^2)
  q_et <- -mu^2 / a^2
  f_e <- th * (y - mu) / a
  f_ee <- -th * mu * (th + y) / a^2
  f_t <- digamma(y + th) - digamma(th) - la + (mu - y) / a
  f_tt <- trigamma(y + th) - trigamma(th) + mu / (th * a) + (y - mu) / a^2
  f_et <- (y - mu) * mu / a^2
  rr <- r * (1 + r)
  l_e <- f_e + r * q_e
  l_ee <- f_ee + rr * q_e^2 + r * q_ee
  l_t <- f_t + r * q_t
  l_tt <- f_tt + rr * q_t^2 + r * q_tt
  l_et <- f_et + rr * q_e * q_t + r * q_et
  list(ll = ll, e = l_e, ee = l_ee, s = th * l_t, ss = th * l_t + th^2 * l_tt,
       es = th * l_et)
}

#' @export
print.panglm_hurdle <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  family_label <- if (identical(x$count_family, "negbin")) "NB2" else "Poisson"
  cat("Fixed-effects hurdle model (zero part: conditional logit; count part: ",
      "zero-truncated ", family_label, ")\n\n", sep = "")
  cat("Zero part, 1(y > 0):\n")
  print(stats::coef(x$zero), digits = digits)
  cat("\nCount part, y | y > 0:\n")
  print(x$count$coefficients, digits = digits)
  if (identical(x$count_family, "negbin")) {
    cat("Theta:", format(x$count$theta, digits = digits), "\n")
  }
  cat("\nn =", x$n, " positive =", x$n_positive, "\n")
  invisible(x)
}
