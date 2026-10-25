#' Extended two-way fixed-effects difference-in-differences
#'
#' The extended two-way fixed-effects (ETWFE) estimator of Wooldridge
#' (2021, 2023) for staggered adoption, with linear, Poisson or logit
#' outcome models. The model interacts the treatment with every cohort by
#' period cell, so that each cohort can have its own effect in each period:
#' with `g` the first treated period of the cohort and `t` the period,
#' \deqn{E[y_{it}] = F(\alpha + \gamma_g + \lambda_t + \sum_{(g, t)}
#' \tau_{gt} D_{gt,it} + x_{it}'\beta + \ldots)}{E[y] = F(alpha + gamma_g
#' + lambda_t + sum tau_gt D_gt + x'beta + ...)}
#' where `D_gt` indicates cohort `g` in period `t`, `gamma_g` and
#' `lambda_t` are cohort and period effects and `F` is the identity, the
#' exponential or the logistic function. Because the cell effects are
#' unrestricted, the estimator does not suffer from the negative weights
#' of the single-coefficient TWFE regression of [did()]; in the linear case
#' with `control = "notyet"` it is numerically equivalent to an imputation
#' estimator. For nonlinear models, parallel trends is imposed on the index
#' (for example on the log of the mean for Poisson), which suits counts and
#' binary outcomes.
#'
#' The average treatment effects on the treated are then computed on the
#' response scale as marginal effects: for each treated observation, the
#' difference between the predicted mean with its cell effect and without
#' it, averaged over the treated observations (all of them, by cohort, by
#' calendar period or by time since treatment; see [summary.panglm_etwfe()]).
#' Standard errors use the delta method with analytic gradients.
#'
#' @section Specification:
#' Cohort and period effects enter as dummies (cohort effects rather than
#' individual effects: for linear models on balanced panels without
#' time-varying controls the two give the same cell effects, and for
#' nonlinear models the pooled
#' quasi-likelihood with cohort dummies is the estimator of Wooldridge
#' 2023). With `control = "notyet"` (the default) the cell dummies cover
#' the post-treatment cells only, so all untreated observations, of
#' never-treated and not-yet-treated cohorts, identify the counterfactual.
#' With `control = "never"` every cell of each treated cohort gets a
#' dummy except the period just before treatment, so only never-treated
#' individuals identify the period effects and the pre-treatment cells
#' give placebo estimates. The reference cohort is the never-treated
#' group; without never-treated individuals and with `control = "notyet"`
#' the last cohort to be treated is the reference and the periods from its
#' treatment date onwards are dropped. Individuals treated in the first
#' period are dropped.
#'
#' Control variables in `formula` enter with their main effect, their
#' interactions with the cohort and period dummies, and interactions of the
#' cell dummies with the control variables centred at their cohort means,
#' as in Wooldridge (2021) and the `etwfe` package, so that the cell
#' effects can vary with the controls.
#'
#' With `ivar = TRUE` (linear model only) individual fixed effects replace
#' the cohort effects, as `ivar` in `etwfe::etwfe()`: they are swept out by
#' within-individual demeaning (weighted when there are weights). On a
#' balanced panel without time-varying controls the cell effects and ATTs
#' are the same as with cohort effects. They differ with time-varying
#' controls (the individual effects then also absorb the part of the
#' individual heterogeneity that is correlated with the controls, as the
#' Mundlak means of Wooldridge 2021 would) and on unbalanced panels; the
#' individual effects also absorb time-invariant controls.
#' Nonlinear models keep the cohort effects (with few periods, individual
#' effects in a logit model suffer from the incidental parameters problem;
#' the pooled Poisson and logit quasi-likelihoods with cohort effects are
#' the estimators of Wooldridge 2023), and `ivar = TRUE` is an error there,
#' whereas `etwfe::etwfe()` ignores `ivar` for nonlinear families.
#'
#' @section Heterogeneity by a covariate:
#' `xvar` names a covariate by which the treatment effects may vary, as
#' `xvar` in `etwfe::etwfe()` (typically a time-invariant indicator, such
#' as sex, rural residence or a baseline risk group). The covariate is
#' centred at its cohort by period cell means and enters through its
#' interactions with the period dummies (except the first period) and with
#' every cell dummy, so each cohort by period effect becomes a linear
#' function of the covariate, with the cell coefficient as the effect at
#' the cell mean. This is the specification of the `etwfe` package; the
#' level of the covariate itself is not a regressor, so add it to
#' `formula` as a control to adjust for it as well. [summary.panglm_etwfe()]
#' then reports the ATTs separately for each value of the covariate
#' (`byxvar = TRUE`, the default for such a fit), as
#' `etwfe::emfx(by_xvar = TRUE)`; with a continuous covariate there is one
#' row per distinct value, so summarise a binary or discrete version of it.
#'
#' @param formula `outcome ~ 1`, or `outcome ~ controls`
#' @param data a data.frame in long format
#' @param index `c(individual, time)`; the time index must be numeric. The
#'   individual identifies the clusters for the standard errors.
#' @param treat name of the column holding each individual's first treated
#'   period; `NA`, `0` or `Inf` for units that are never treated
#' @param family `"gaussian"` (linear), `"poisson"` (counts, log link) or
#'   `"binomial"` (binary, logit link; `"logit"` is accepted)
#' @param control `"notyet"` or `"never"`; see Specification
#' @param vcov `"cluster"` (by individual, the default) or `"classical"`.
#'   The cluster covariance uses the factor
#'   `G / (G - 1) * (n - 1) / (n - k)` with `k` counting all coefficients,
#'   cohort and period dummies included (with `ivar = TRUE`, the individual
#'   effects, nested in the clusters, count as one), as `fixest` does.
#' @param level confidence level of the reported intervals
#' @param xvar optional name of one numeric or logical column of `data`: a
#'   covariate by which the treatment effects may vary (see Heterogeneity
#'   by a covariate)
#' @param ivar logical; `TRUE` uses individual (the first element of
#'   `index`) instead of cohort fixed effects, for `family = "gaussian"`
#'   only (see Specification)
#' @param weights optional sampling weights: the name of a column of `data`
#'   or a non-negative numeric vector with one value per row (rows with a
#'   missing weight are dropped). They weight the quasi-likelihood (as
#'   `weights = ~w` in `etwfe::etwfe()` and `fixest`) and the averages of
#'   the treatment effects over the treated observations in
#'   [summary.panglm_etwfe()] (as `wts` in `etwfe::emfx()`). The centring
#'   of the controls at their cohort means (and of `xvar` at its cell
#'   means) is unweighted, as in `etwfe`. Unlike in [attgt()], the weights
#'   need not be constant within individuals: the model is a weighted
#'   pooled (quasi-)likelihood over the rows, so observation-level weights
#'   are well defined, and the cluster-robust covariance by individual
#'   remains valid (as in `fixest`). Use weights constant within
#'   individuals when they are sampling weights of the individuals, so that
#'   the estimates are comparable with [attgt()].
#' @return an object of class `"panglm_etwfe"`: a list with
#'   `coefficients`, `vcov`, `cells` (a data frame of the cohort by period
#'   cells with their coefficient on the index scale), `family`, `control`,
#'   `reference` (the reference cohort), `xvar`, `ivar`, `nobs`,
#'   `n_groups`, the design
#'   (`X`, `y`, `rows`) and the call. Methods: `print`, `summary` (the
#'   ATTs), `coef`, `vcov`, `nobs`, `plot`.
#' @seealso [attgt()] for the Callaway and Sant'Anna estimator, [did()] and
#'   [eventstudy()] for the TWFE versions.
#' @references Wooldridge, J. M. (2021). Two-way fixed effects, the
#'   two-way Mundlak regression, and difference-in-differences estimators.
#'   Working paper. \doi{10.2139/ssrn.3906345}
#'
#'   Wooldridge, J. M. (2023). Simple approaches to nonlinear
#'   difference-in-differences with panel data. *The Econometrics Journal*,
#'   26, C31-C66. \doi{10.1093/ectj/utad016}
#' @family policy evaluation
#' @examples
#' # Staggered rollout of a health insurance scheme by region (for example
#' # RAMED, then AMO, in Morocco): regions join in 2016, 2018 or 2020, the
#' # last region not before 2022. Outcome: number of hospital admissions
#' # per district and year (a count).
#' set.seed(1)
#' d <- expand.grid(district = 1:160, year = 2012:2021)
#' d$adoption <- rep(c(2016, 2018, 2020, NA), each = 40)[d$district]
#' d$poverty <- rep(runif(160), 10)[d$district]
#' a <- rnorm(160, sd = 0.3)
#' treated <- !is.na(d$adoption) & d$year >= d$adoption
#' d$admissions <- rpois(nrow(d), exp(2 + a[d$district] + 0.03 * (d$year - 2012) +
#'                                     0.2 * treated))
#' fit <- etwfe(admissions ~ 1, data = d, index = c("district", "year"),
#'              treat = "adoption", family = "poisson")
#' fit                                   # overall ATT, admissions per district
#' summary(fit, type = "event")          # by time since adoption
#' summary(fit, type = "group")
#' plot(fit)
#' # linear model with a control variable
#' d$visits <- 50 + 5 * a[d$district] + 0.5 * (d$year - 2012) + 2 * treated +
#'   rnorm(nrow(d))
#' summary(etwfe(visits ~ poverty, data = d, index = c("district", "year"),
#'               treat = "adoption"), type = "calendar")
#' # individual instead of cohort fixed effects (linear model)
#' summary(etwfe(visits ~ poverty, data = d, index = c("district", "year"),
#'               treat = "adoption", ivar = TRUE), type = "group")
#' # effects by a covariate: poor versus other districts
#' d$poor <- as.numeric(d$poverty > 0.5)
#' fx <- etwfe(admissions ~ 1, data = d, index = c("district", "year"),
#'             treat = "adoption", family = "poisson", xvar = "poor")
#' summary(fx)                           # ATT for poor = 0 and poor = 1
#' summary(fx, type = "event", byxvar = FALSE)
#' @export
etwfe <- function(formula, data, index, treat,
                  family = c("gaussian", "poisson", "binomial", "logit"),
                  control = c("notyet", "never"), vcov = c("cluster", "classical"),
                  level = 0.95, weights = NULL, xvar = NULL, ivar = FALSE) {
  family <- match.arg(family)
  if (family == "logit") family <- "binomial"
  control <- match.arg(control)
  vcov <- match.arg(vcov)
  data <- as.data.frame(data)
  did_timing(data, index, treat)
  formula <- stats::as.formula(formula)
  if (length(formula) != 3L) stop("'formula' must have an outcome, e.g. y ~ 1", call. = FALSE)
  if (!is.logical(ivar) || length(ivar) != 1L || is.na(ivar)) {
    stop("'ivar' must be TRUE or FALSE", call. = FALSE)
  }
  if (ivar && family != "gaussian") {
    stop("ivar = TRUE (individual fixed effects) is available for family = 'gaussian' only; ",
         "nonlinear models use cohort effects (Wooldridge 2023)", call. = FALSE)
  }
  if (!is.null(xvar)) {
    if (!is.character(xvar) || length(xvar) != 1L) {
      stop("'xvar' must be the name of one column of data", call. = FALSE)
    }
    if (!xvar %in% names(data)) stop("xvar column '", xvar, "' not found in data", call. = FALSE)
    if (is.logical(data[[xvar]])) data[[xvar]] <- as.numeric(data[[xvar]])
    if (!is.numeric(data[[xvar]])) {
      stop("'xvar' must be numeric or logical (code a categorical moderator as 0/1 ",
           "indicators)", call. = FALSE)
    }
  }
  vars <- unique(c(all.vars(formula), index, xvar))
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars)) {
    stop("column(s) not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  wv <- did_weights(weights, data)
  ok <- stats::complete.cases(data[, vars, drop = FALSE])
  if (!is.null(wv)) ok <- ok & !is.na(wv)
  if (!all(ok)) warning("dropped ", sum(!ok), " rows with missing values", call. = FALSE)
  data <- data[ok, , drop = FALSE]
  wt <- if (is.null(wv)) rep(1, nrow(data)) else wv[ok]
  mf <- stats::model.frame(formula, data = data)
  y <- as.numeric(stats::model.response(mf))
  fam <- switch(family, gaussian = stats::gaussian(), poisson = stats::poisson(),
                binomial = stats::binomial())
  if (family == "poisson" && any(y < 0)) stop("a Poisson outcome must be non-negative", call. = FALSE)
  if (family == "binomial" && !all(y %in% c(0, 1))) {
    stop("a binary outcome must be coded 0/1", call. = FALSE)
  }
  ctrl <- stats::model.matrix(stats::delete.response(stats::terms(formula)), mf)
  ctrl <- ctrl[, colnames(ctrl) != "(Intercept)", drop = FALSE]
  id <- data[[index[1]]]
  xv <- if (is.null(xvar)) NULL else data[[xvar]]
  time <- data[[index[2]]]
  first <- data[[treat]]
  tlist <- sort(unique(time))
  first[is.na(first) | !is.finite(first) | first == 0 | first > max(tlist)] <- Inf
  always <- first <= min(tlist)
  if (any(always)) {
    warning("dropped ", length(unique(id[always])),
            " individuals already treated in the first period", call. = FALSE)
  }
  keep <- !always
  # controls are centred at their cohort means before any period is dropped
  ctrl_dm <- ctrl
  if (ncol(ctrl)) {
    cohort_key <- as.character(first)
    for (j in seq_len(ncol(ctrl))) {
      ctrl_dm[keep, j] <- ctrl[keep, j] - stats::ave(ctrl[keep, j], cohort_key[keep])
    }
  }
  # xvar is centred at its cohort by period cell means (as in etwfe)
  xv_dm <- if (is.null(xv)) NULL else xv - stats::ave(xv, paste(first, time))
  if (any(is.infinite(first[keep]))) {
    gref <- Inf
  } else if (control == "notyet") {
    gref <- max(first[keep])
    keep <- keep & time < gref
    warning("no never-treated individuals: cohort ", gref, " is the reference and periods from ",
            gref, " on are dropped", call. = FALSE)
  } else {
    stop("control = 'never' needs never-treated individuals", call. = FALSE)
  }
  id <- id[keep]; time <- time[keep]; first <- first[keep]; y <- y[keep]; wt <- wt[keep]
  if (!is.null(xv)) {
    xv <- xv[keep]
    xv_dm <- matrix(xv_dm[keep], ncol = 1L, dimnames = list(NULL, paste0(xvar, "_dm")))
  }
  if (!(sum(wt) > 0)) stop("the weights must have a positive sum", call. = FALSE)
  ctrl <- ctrl[keep, , drop = FALSE]; ctrl_dm <- ctrl_dm[keep, , drop = FALSE]
  tlist <- sort(unique(time))
  tref <- min(tlist)
  cohorts <- sort(unique(first[first != gref]))
  if (!length(cohorts)) stop("no treated cohort", call. = FALSE)

  # cells: cohort g in period t
  cells <- do.call(rbind, lapply(cohorts, function(g) {
    tt <- if (control == "notyet") {
      tlist[tlist >= g & tlist != tref]
    } else {
      prev <- max(tlist[tlist < g])
      tlist[tlist != prev]
    }
    if (length(tt)) data.frame(group = g, time = tt) else NULL
  }))
  cells$name <- paste0("cell::", cells$group, ":", cells$time)
  row_cell <- match(paste(first, time), paste(cells$group, cells$time))
  present <- seq_len(nrow(cells)) %in% row_cell
  cells <- cells[present, , drop = FALSE]
  row_cell <- match(paste(first, time), paste(cells$group, cells$time))

  n <- length(y)
  cn <- colnames(ctrl)
  dummies <- function(v, levels, prefix) {
    m <- outer(v, levels, `==`) * 1
    colnames(m) <- paste0(prefix, levels)
    m
  }
  Dg <- dummies(first, cohorts, "cohort::")
  Dt <- dummies(time, setdiff(tlist, tref), "period::")
  inter <- function(D, M) {
    if (!ncol(M) || !ncol(D)) return(NULL)
    do.call(cbind, lapply(seq_len(ncol(M)), function(j) {
      out <- D * M[, j]
      colnames(out) <- paste0(colnames(D), ":", colnames(M)[j])
      out
    }))
  }
  Dc <- matrix(0, n, nrow(cells), dimnames = list(NULL, cells$name))
  has <- !is.na(row_cell)
  Dc[cbind(which(has), row_cell[has])] <- 1
  if (is.null(xv)) {
    Dcx <- Dtx <- NULL
  } else {
    Dtx <- inter(Dt, xv_dm)
    Dcx <- inter(Dc, xv_dm)
  }
  X <- cbind(`(Intercept)` = 1, Dg, Dt, ctrl, inter(Dg, ctrl), inter(Dt, ctrl), Dtx,
             Dc, inter(Dc, ctrl_dm), Dcx)
  is_cell <- colnames(X) %in% c(cells$name, colnames(inter(Dc, ctrl_dm)), colnames(Dcx))

  # with individual fixed effects (linear model only), the individual
  # effects are swept out by (weighted) within-individual demeaning; the
  # cell effects are differences of linear predictions, so the ATTs use the
  # undemeaned design
  if (ivar) {
    wsum <- as.numeric(rowsum(wt, id))
    gi <- match(id, sort(unique(id)))
    within <- function(M) M - (rowsum(wt * M, id) / ifelse(wsum > 0, wsum, 1))[gi, , drop = FALSE]
    Xf <- within(X)
    # columns constant within individuals (intercept, cohort dummies) are
    # absorbed: set their rounding residue to zero so that they are dropped
    absorbed <- sqrt(colSums(Xf^2)) <= 1e-9 * pmax(sqrt(colSums(X^2)), 1e-300)
    Xf[, absorbed] <- 0
    yf <- as.numeric(within(matrix(y)))
  } else {
    Xf <- X
    yf <- y
  }
  fit <- etwfe_fit(Xf, yf, fam, wt)
  X <- X[, fit$keep, drop = FALSE]
  Xf <- Xf[, fit$keep, drop = FALSE]
  is_cell <- is_cell[fit$keep]
  b <- fit$coefficients
  eta <- as.numeric(Xf %*% b)
  mu <- fam$linkinv(eta)
  w <- wt * fam$mu.eta(eta)^2 / fam$variance(mu)
  XtWX <- crossprod(Xf, w * Xf)
  bread <- chol2inv(chol(XtWX))
  k <- ncol(Xf)
  G <- length(unique(id))
  if (vcov == "cluster") {
    # individual effects are nested in the clusters: as in fixest, they
    # count as one parameter (the intercept they absorb)
    kc <- if (ivar) k + 1 else k
    score <- Xf * (wt * (yf - mu) * fam$mu.eta(eta) / fam$variance(mu))
    meat <- crossprod(rowsum(score, id))
    V <- bread %*% meat %*% bread * (G / (G - 1)) * ((n - 1) / (n - kc))
  } else {
    kc <- if (ivar) k + G else k
    disp <- if (family == "gaussian") sum(wt * (yf - mu)^2) / (n - kc) else 1
    V <- bread * disp
  }
  dimnames(V) <- list(colnames(X), colnames(X))
  cells$estimate <- b[cells$name]
  cells$se <- sqrt(diag(V))[cells$name]
  cells$name <- NULL
  rownames(cells) <- NULL
  if (anyNA(cells$estimate)) cells <- cells[!is.na(cells$estimate), , drop = FALSE]

  structure(list(coefficients = b, vcov = V, vcov_type = vcov, cells = cells,
                 family = fam, control = control, reference = gref, X = X, y = y,
                 weights = if (is.null(wv)) NULL else wt,
                 is_cell = is_cell, xvar = xvar, ivar = ivar,
                 rows = data.frame(id = id, time = time, group = first, cell = has,
                                   xvar = if (is.null(xv)) NA_real_ else xv),
                 nobs = n, n_groups = G, level = level, outcome = all.vars(formula[[2]]),
                 controls = cn, converged = fit$converged, call = match.call()),
            class = "panglm_etwfe")
}

#' GLM fit for etwfe(), dropping aliased columns
#'
#' @keywords internal
#' @noRd
etwfe_fit <- function(X, y, fam, wt = NULL) {
  keep <- rep(TRUE, ncol(X))
  q <- qr(X)
  if (q$rank < ncol(X)) keep[q$pivot[seq.int(q$rank + 1L, ncol(X))]] <- FALSE
  Xk <- X[, keep, drop = FALSE]
  fit <- suppressWarnings(stats::glm.fit(Xk, y, weights = wt, family = fam,
                                         control = list(epsilon = 1e-12, maxit = 100)))
  if (!fit$converged) warning("the model did not converge", call. = FALSE)
  list(coefficients = stats::setNames(fit$coefficients, colnames(Xk)), keep = keep,
       converged = fit$converged)
}

#' Average treatment effects on the treated from an ETWFE fit
#'
#' Averages, over the treated observations, the difference between the
#' predicted mean with and without the observation's cohort by period
#' effect (on the response scale). The averages are taken over all treated
#' observations (`"simple"`), by cohort (`"group"`), by calendar period
#' (`"calendar"`) or by time since treatment (`"event"`), with delta-method
#' standard errors. With sampling weights in [etwfe()], the averages are
#' weighted. They correspond to `etwfe::emfx()`. With
#' `control = "never"`, the event study also reports the pre-treatment
#' placebo effects, with event time `-1` as the reference (zero).
#'
#' @param object an `"panglm_etwfe"` fit from [etwfe()]
#' @param type `"simple"`, `"group"`, `"calendar"` or `"event"`
#' @param level confidence level; defaults to that of the fit
#' @param byxvar logical; for a fit with `xvar`, report the ATTs separately
#'   for each value of the covariate (by aggregation level and value). The
#'   default is `TRUE` for such fits; `overall` is always the ATT over all
#'   treated observations.
#' @param ... unused
#' @return an object of class `"summary.panglm_etwfe"`: a list with
#'   `overall` (the simple ATT: `estimate`, `se`, `statistic`, `p.value`,
#'   `lower`, `upper`), `estimates` (the same columns by `level` of the
#'   aggregation, and by value of the `xvar` covariate in a column named
#'   after it when `byxvar = TRUE`; `NULL` for `"simple"` without
#'   `byxvar`), `cells` and the fit information.
#' @family etwfe methods
#' @examples
#' set.seed(2)
#' d <- expand.grid(id = 1:150, t = 1:8)
#' d$first <- rep(c(4, 6, NA), each = 50)[d$id]
#' d$y <- rbinom(nrow(d), 1, plogis(-1 + 0.1 * d$t +
#'   0.8 * (!is.na(d$first) & d$t >= d$first)))
#' fit <- etwfe(y ~ 1, d, index = c("id", "t"), treat = "first", family = "binomial")
#' summary(fit, type = "event")
#' @export
summary.panglm_etwfe <- function(object, type = c("simple", "group", "calendar", "event"),
                                 level = object$level, byxvar = !is.null(object$xvar), ...) {
  type <- match.arg(type)
  if (!is.logical(byxvar) || length(byxvar) != 1L || is.na(byxvar)) {
    stop("'byxvar' must be TRUE or FALSE", call. = FALSE)
  }
  if (byxvar && is.null(object$xvar)) {
    stop("byxvar = TRUE needs a fit with 'xvar'", call. = FALSE)
  }
  rows <- object$rows
  treated <- is.finite(rows$group) & rows$group != object$reference
  use_overall <- treated & rows$time >= rows$group & rows$cell
  overall <- etwfe_ate(object, use_overall, rep(1, sum(use_overall)), level, overall = TRUE)
  by <- NULL
  if (type != "simple" || byxvar) {
    use <- if (type == "event" && object$control == "never") treated else use_overall
    key <- switch(type, simple = rep(1, nrow(rows)), group = rows$group, calendar = rows$time,
                  event = rows$time - rows$group)[use]
    by <- etwfe_ate(object, use, key, level, overall = type == "simple",
                    key2 = if (byxvar) rows$xvar[use] else NULL)
    if (byxvar) names(by)[names(by) == "xvar"] <- object$xvar
  }
  out <- list(overall = overall, estimates = by, type = type, xvar = if (byxvar) object$xvar,
              cells = object$cells,
              family = object$family$family, control = object$control,
              vcov_type = object$vcov_type, ivar = isTRUE(object$ivar), nobs = object$nobs,
              n_groups = object$n_groups,
              level = level, call = object$call)
  class(out) <- "summary.panglm_etwfe"
  out
}

etwfe_ate <- function(object, use, key, level, overall = FALSE, key2 = NULL) {
  X1 <- object$X[use, , drop = FALSE]
  X0 <- X1
  X0[, object$is_cell] <- 0
  b <- object$coefficients
  fam <- object$family
  wt <- if (is.null(object$weights)) NULL else object$weights[use]
  eta1 <- as.numeric(X1 %*% b)
  eta0 <- as.numeric(X0 %*% b)
  effect <- fam$linkinv(eta1) - fam$linkinv(eta0)
  grad <- X1 * fam$mu.eta(eta1) - X0 * fam$mu.eta(eta0)
  if (is.null(key2)) {
    lev <- data.frame(level = sort(unique(key)))
    grp <- match(key, lev$level)
  } else {
    lev <- unique(data.frame(level = key, xvar = key2))
    lev <- lev[order(lev$level, lev$xvar), , drop = FALSE]
    grp <- match(paste(key, key2), paste(lev$level, lev$xvar))
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  rows <- lapply(seq_len(nrow(lev)), function(l) {
    s <- grp == l
    if (is.null(wt)) {
      est <- mean(effect[s])
      g <- colMeans(grad[s, , drop = FALSE])
    } else {
      ws <- wt[s] / sum(wt[s])
      est <- sum(ws * effect[s])
      g <- colSums(ws * grad[s, , drop = FALSE])
    }
    se <- sqrt(max(0, as.numeric(t(g) %*% object$vcov %*% g)))
    c(estimate = est, se = se)
  })
  out <- cbind(lev, as.data.frame(do.call(rbind, rows)))
  rownames(out) <- NULL
  out$statistic <- out$estimate / out$se
  out$p.value <- 2 * stats::pnorm(-abs(out$statistic))
  out$statistic[out$se == 0] <- NA
  out$p.value[out$se == 0] <- NA
  out$lower <- out$estimate - z * out$se
  out$upper <- out$estimate + z * out$se
  if (overall) out$level <- NULL
  out
}

#' @export
print.summary.panglm_etwfe <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  fam <- c(gaussian = "linear", poisson = "Poisson", binomial = "logit")[[x$family]]
  cat("\nExtended TWFE (Wooldridge), ", fam, " model, control group: ",
      if (x$control == "never") "never treated" else "not yet treated", "\n", sep = "")
  cat("Observations:", x$nobs, " Individuals:", x$n_groups, " vcov:", x$vcov_type, "\n")
  cat("Fixed effects: ", if (isTRUE(x$ivar)) "individual" else "cohort", " and period",
      if (!is.null(x$xvar)) paste0("; effects by ", x$xvar), "\n\n", sep = "")
  ov <- x$overall
  cat("ATT (response scale): ", format(ov$estimate, digits = digits), " (SE ",
      format(ov$se, digits = digits), ", ", format(100 * x$level), "% CI ",
      format(ov$lower, digits = digits), " to ", format(ov$upper, digits = digits),
      ", p = ", format.pval(ov$p.value, digits = digits), ")\n", sep = "")
  if (!is.null(x$estimates)) {
    tab <- x$estimates
    if (x$type != "simple") {
      names(tab)[1] <- c(group = "group", calendar = "time", event = "event")[[x$type]]
    }
    cat("\n")
    print(format(tab, digits = digits), row.names = FALSE)
  }
  cat("\n")
  invisible(x)
}

#' @export
print.panglm_etwfe <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  s <- summary.panglm_etwfe(x)
  print(s, digits = digits)
  cat("Use summary(x, type = ...) for effects by group, calendar period or event time.\n\n")
  invisible(x)
}

#' @export
coef.panglm_etwfe <- function(object, ...) object$coefficients

#' @export
vcov.panglm_etwfe <- function(object, ...) object$vcov

#' @export
nobs.panglm_etwfe <- function(object, ...) object$nobs

#' Plot ETWFE treatment effects
#'
#' Plots the event-study aggregation of an [etwfe()] fit (or the cohort or
#' calendar aggregation) with pointwise confidence intervals.
#'
#' @param x an `"panglm_etwfe"` fit from [etwfe()], or its summary
#' @param type aggregation to plot (for a fit): `"event"`, `"group"` or
#'   `"calendar"`
#' @param xlab,ylab axis labels
#' @param ... passed to [graphics::plot()]
#' @return the plotted summary, invisibly
#' @family etwfe methods
#' @examples
#' set.seed(2)
#' d <- expand.grid(id = 1:150, t = 1:8)
#' d$first <- rep(c(4, 6, NA), each = 50)[d$id]
#' d$y <- rnorm(150)[d$id] + 0.2 * d$t +
#'   ifelse(!is.na(d$first) & d$t >= d$first, 1, 0) + rnorm(nrow(d))
#' plot(etwfe(y ~ 1, d, index = c("id", "t"), treat = "first", control = "never"))
#' @export
plot.panglm_etwfe <- function(x, type = c("event", "group", "calendar"), xlab = NULL,
                              ylab = "ATT", ...) {
  type <- match.arg(type)
  plot.summary.panglm_etwfe(summary.panglm_etwfe(x, type = type, byxvar = FALSE), xlab = xlab,
                            ylab = ylab, ...)
}

#' @rdname plot.panglm_etwfe
#' @export
plot.summary.panglm_etwfe <- function(x, type = NULL, xlab = NULL, ylab = "ATT", ...) {
  tab <- x$estimates
  if (is.null(tab) || x$type == "simple") {
    stop("nothing to plot for the simple aggregation", call. = FALSE)
  }
  if (!is.null(x$xvar)) {
    stop("plot a summary with byxvar = FALSE, or subset the estimates by '", x$xvar, "'",
         call. = FALSE)
  }
  s <- list(estimates = tab[, c("level", "estimate", "se", "lower", "upper")],
            type = c(event = "dynamic", group = "group", calendar = "calendar")[[x$type]])
  plot.summary.panglm_attgt(s, xlab = xlab, ylab = ylab, ...)
}
