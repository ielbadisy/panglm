#' Group-time average treatment effects for staggered adoption
#'
#' The Callaway and Sant'Anna (2021) estimator of the group-time average
#' treatment effects on the treated `ATT(g, t)`: the effect at period `t`
#' for the cohort of individuals first treated at period `g`. Each
#' `ATT(g, t)` is a two-group, two-period comparison of the change in the
#' outcome between a base period and period `t` for cohort `g` against a
#' clean control group (never-treated individuals, or individuals not yet
#' treated by period `t`), so that already-treated individuals are never
#' used as controls. This avoids the negative weights of the two-way
#' fixed-effects regression of [did()] under staggered adoption with
#' heterogeneous effects. The `ATT(g, t)` are then aggregated into an
#' overall effect, an event study (effect by time since treatment), effects
#' by cohort or by calendar period with [summary.panglm_attgt()].
#'
#' @section Estimators:
#' With `y ~ 1` (no covariates) every `ATT(g, t)` is the difference between
#' the mean outcome change of cohort `g` and that of the control group,
#' whatever `method`. With covariates, parallel trends is assumed to hold
#' conditionally on them; their values in the base period (the earlier of
#' the two periods compared) are used, and the comparison is made by
#' * `method = "dr"`: the doubly robust estimator of Sant'Anna and Zhao
#'   (2020): a logit propensity score for belonging to cohort `g` and a
#'   linear regression of the outcome change among controls; consistent if
#'   either model is correct;
#' * `method = "ipw"`: normalized inverse probability weighting with the
#'   logit propensity score;
#' * `method = "reg"`: outcome regression (the linear regression of the
#'   outcome change fitted on the controls, averaged over the treated).
#'
#' Controls with a propensity score of 0.995 or more are trimmed, as in the
#' `DRDID` package. The panel is balanced first (individuals not observed in
#' every period are dropped), individuals already treated in the first
#' period are dropped, and cohorts first treated after the last period are
#' treated as never treated. When no individual is never treated, the last
#' cohort to be treated serves as the never-treated group and the periods
#' from its treatment date onwards are dropped (the convention of
#' `did::att_gt()`).
#'
#' With `base = "varying"` (the default) pre-treatment `ATT(g, t)` compare
#' period `t` with period `t - 1`, which gives placebo estimates of the
#' change in trends; with `base = "universal"` every `ATT(g, t)` is
#' measured against the last period before treatment, `g - 1 - anticipation`,
#' whose own `ATT` is zero. Post-treatment effects are the same either way.
#'
#' @section Inference:
#' Standard errors are analytic and come from the influence function of
#' each estimator, including the estimation of the propensity score and
#' outcome regression and, for the aggregations, of the cohort shares used
#' as weights; they are robust to serial correlation within individuals.
#' They match `did::att_gt(bstrap = FALSE)` and `did::aggte()`. Confidence
#' intervals are pointwise. The joint Wald test of all pre-treatment
#' `ATT(g, t)` is reported as a test of parallel pre-trends.
#'
#' @param formula `outcome ~ 1`, or `outcome ~ covariates` for conditional
#'   parallel trends (covariates are taken in the base period)
#' @param data a data.frame in long format
#' @param index `c(individual, time)`; the time index must be numeric
#' @param treat name of the column holding each individual's first treated
#'   period (in the units of the time index); `NA`, `0` or `Inf` for units
#'   that are never treated
#' @param control comparison group: `"never"` (never-treated individuals)
#'   or `"notyet"` (never treated or not yet treated, nor anticipating
#'   treatment, at both periods compared)
#' @param method `"dr"` (doubly robust), `"ipw"` or `"reg"`; only matters
#'   with covariates
#' @param base base period for pre-treatment `ATT(g, t)`: `"varying"` or
#'   `"universal"`
#' @param anticipation number of periods before `g` in which treatment can
#'   already affect the outcome (the base period moves back accordingly)
#' @param level confidence level of the reported intervals
#' @return an object of class `"panglm_attgt"`: a list with `estimates` (a
#'   data frame with `group`, `time`, `event` (`time - group`), `estimate`,
#'   `se`, `lower`, `upper` and `post`), `inf` (the matrix of influence
#'   functions, one row per individual and one column per `ATT(g, t)`),
#'   `units` (individual identifiers and cohorts), `pretrend` (the Wald test
#'   of the pre-treatment `ATT(g, t)`, or `NULL`), `n` (individuals),
#'   `groups`, `periods` and the settings. Methods: `print`, `summary`,
#'   `coef`, `vcov`, `nobs`, `plot`.
#' @seealso [summary.panglm_attgt()] for the aggregations, [etwfe()] for
#'   the extended two-way fixed-effects alternative (also for count and
#'   binary outcomes), [did()] and [eventstudy()] for the TWFE versions.
#' @references Callaway, B. and Sant'Anna, P. H. C. (2021).
#'   Difference-in-differences with multiple time periods. *Journal of
#'   Econometrics*, 225, 200-230. \doi{10.1016/j.jeconom.2020.12.001}
#'
#'   Sant'Anna, P. H. C. and Zhao, J. (2020). Doubly robust
#'   difference-in-differences estimators. *Journal of Econometrics*, 219,
#'   101-122. \doi{10.1016/j.jeconom.2020.06.003}
#' @examples
#' # A health insurance scheme (for example the AMO extension in Morocco)
#' # is rolled out region by region: provinces of region A join in 2016,
#' # region B in 2018, region C in 2020, and region D not before 2022.
#' # Outcome: outpatient visits per 1,000 inhabitants.
#' set.seed(1)
#' d <- expand.grid(province = 1:120, year = 2012:2021)
#' region <- rep(c(2016, 2018, 2020, NA), each = 30)
#' d$adoption <- region[d$province]
#' d$urban <- rep(runif(120), 10)[d$province]
#' a <- rnorm(120)
#' e <- d$year - d$adoption
#' d$visits <- 50 + 5 * a[d$province] + 0.5 * (d$year - 2012) * (1 + d$urban) +
#'   ifelse(!is.na(e) & e >= 0, 2 + 0.5 * e, 0) + rnorm(nrow(d))
#' fit <- attgt(visits ~ 1, data = d, index = c("province", "year"),
#'              treat = "adoption")
#' fit
#' summary(fit)                     # overall ATT
#' summary(fit, type = "dynamic")   # event study
#' summary(fit, type = "group")
#' plot(fit)
#' # conditional parallel trends (trends differ by urbanization)
#' fitx <- attgt(visits ~ urban, data = d, index = c("province", "year"),
#'               treat = "adoption", control = "notyet", method = "dr")
#' summary(fitx, type = "dynamic")
#' @export
attgt <- function(formula, data, index, treat, control = c("never", "notyet"),
                  method = c("dr", "ipw", "reg"), base = c("varying", "universal"),
                  anticipation = 0, level = 0.95) {
  control <- match.arg(control)
  method <- match.arg(method)
  base <- match.arg(base)
  if (!is.numeric(anticipation) || length(anticipation) != 1L || anticipation < 0) {
    stop("'anticipation' must be a non-negative number", call. = FALSE)
  }
  data <- as.data.frame(data)
  did_timing(data, index, treat)
  prep <- attgt_prepare(formula, data, index, treat, control, anticipation)
  n <- prep$n
  tlist <- prep$periods
  glist <- prep$groups
  G <- prep$G
  nT <- length(tlist)
  cells <- list()
  inf <- list()
  for (g in glist) {
    before <- which(tlist + anticipation < g)
    pret_g <- before[length(before)]
    tidx <- if (base == "varying") seq.int(2L, nT) else seq_len(nT)
    for (tc in tidx) {
      post <- g <= tlist[tc]
      b <- if (base == "universal" || post) pret_g else tc - 1L
      if (base == "universal" && b == tc) {
        cells[[length(cells) + 1L]] <- c(g, tlist[tc], 0, 0)
        inf[[length(inf) + 1L]] <- numeric(n)
        next
      }
      ctrl <- if (control == "never") {
        G == 0
      } else {
        G == 0 | (G > max(tlist[tc], tlist[b]) + anticipation & G != g)
      }
      keep <- which(G == g | ctrl)
      D <- as.numeric(G[keep] == g)
      res <- if (sum(D) == 0 || sum(1 - D) == 0) {
        warning("no control units for group ", g, " in period ", tlist[tc],
                "; ATT(g, t) set to NA", call. = FALSE)
        NULL
      } else {
        X <- prep$X[[min(tc, b)]][keep, , drop = FALSE]
        tryCatch(attgt_2x2(prep$Y[keep, tc] - prep$Y[keep, b], D, X, method),
                 error = function(e) {
                   warning("ATT(g, t) for group ", g, " in period ", tlist[tc],
                           " set to NA: ", conditionMessage(e), call. = FALSE)
                   NULL
                 })
      }
      phi <- rep(NA_real_, n)
      att <- NA_real_
      if (!is.null(res)) {
        phi <- numeric(n)
        phi[keep] <- (n / length(keep)) * res$inf
        att <- res$att
      }
      cells[[length(cells) + 1L]] <- c(g, tlist[tc], att, as.numeric(post))
      inf[[length(inf) + 1L]] <- phi
    }
  }
  cells <- do.call(rbind, cells)
  inf <- do.call(cbind, inf)
  se <- sqrt(colSums(inf^2)) / n
  se[!is.na(se) & se <= sqrt(.Machine$double.eps) * 10] <- NA
  z <- stats::qnorm(1 - (1 - level) / 2)
  est <- data.frame(group = cells[, 1], time = cells[, 2], event = cells[, 2] - cells[, 1],
                    estimate = cells[, 3], se = se)
  est$lower <- est$estimate - z * est$se
  est$upper <- est$estimate + z * est$se
  est$post <- cells[, 4] == 1
  colnames(inf) <- paste0("g", est$group, ".t", est$time)
  rownames(inf) <- as.character(prep$ids)

  pre <- which(!est$post & !is.na(est$se))
  pretrend <- NULL
  if (length(pre)) {
    V <- crossprod(inf[, pre, drop = FALSE]) / n^2
    if (rcond(V) > .Machine$double.eps) {
      stat <- as.numeric(t(est$estimate[pre]) %*% solve(V, est$estimate[pre]))
      pretrend <- list(statistic = stat, df = length(pre),
                       p.value = stats::pchisq(stat, length(pre), lower.tail = FALSE),
                       hypothesis = "all pre-treatment ATT(g, t) are zero")
    }
  }
  structure(list(estimates = est, inf = inf,
                 units = data.frame(id = prep$ids, group = G),
                 pretrend = pretrend, n = n, groups = glist, periods = tlist,
                 control = control, method = method, base = base,
                 anticipation = anticipation, level = level,
                 covariates = prep$covariates, outcome = prep$outcome,
                 dropped = prep$dropped, call = match.call()),
            class = "panglm_attgt")
}

#' Balanced wide data for attgt()
#'
#' @keywords internal
#' @noRd
attgt_prepare <- function(formula, data, index, treat, control, anticipation) {
  formula <- stats::as.formula(formula)
  if (length(formula) != 3L) stop("'formula' must have an outcome, e.g. y ~ 1", call. = FALSE)
  rhs <- stats::delete.response(stats::terms(formula))
  outcome <- all.vars(formula[[2]])
  vars <- unique(c(all.vars(formula), index))
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars)) {
    stop("column(s) not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  ok <- stats::complete.cases(data[, vars, drop = FALSE])
  mf <- stats::model.frame(formula, data = data[ok, , drop = FALSE], na.action = stats::na.pass)
  y <- stats::model.response(mf)
  if (!is.numeric(y) && !is.logical(y)) stop("the outcome must be numeric", call. = FALSE)
  X <- stats::model.matrix(rhs, mf)
  fin <- is.finite(as.numeric(y)) & rowSums(!is.finite(X)) == 0 &
    is.finite(data[[index[2]]][ok])
  ok[ok] <- fin
  dropped <- list(rows = sum(!ok))
  if (!all(ok)) {
    warning("dropped ", sum(!ok), " rows with missing or non-finite values", call. = FALSE)
  }
  id <- data[[index[1]]][ok]
  time <- data[[index[2]]][ok]
  first <- data[[treat]][ok]
  first[is.na(first) | !is.finite(first) | first == 0] <- 0
  y <- as.numeric(y[fin])
  X <- X[fin, , drop = FALSE]
  if (anyDuplicated(data.frame(id, time))) {
    stop("each individual must be observed at most once per period", call. = FALSE)
  }
  tlist <- sort(unique(time))
  first[first > max(tlist) + anticipation] <- 0

  keep <- rep(TRUE, length(id))
  glist <- sort(unique(first[first > 0]))
  if (!any(first == 0)) {
    latest <- max(glist)
    keep <- time < latest - anticipation
    if (control == "never") {
      warning("no never-treated individuals: the last cohort (", latest,
              ") is used as the never-treated group and periods from ", latest - anticipation,
              " on are dropped", call. = FALSE)
      first[first == latest] <- 0
    } else {
      warning("no never-treated individuals: periods from ", latest - anticipation,
              " on are dropped and the last cohort (", latest, ") serves only as a control",
              call. = FALSE)
    }
  }
  early <- first > 0 & first <= min(time[keep]) + anticipation
  if (any(early & keep)) {
    warning("dropped ", length(unique(id[early & keep])),
            " individuals already treated in the first period", call. = FALSE)
    dropped$early <- length(unique(id[early & keep]))
  }
  keep <- keep & !early
  id <- id[keep]; time <- time[keep]; first <- first[keep]; y <- y[keep]
  X <- X[keep, , drop = FALSE]
  tlist <- sort(unique(time))
  nT <- length(tlist)
  counts <- table(id)
  full <- names(counts)[counts == nT]
  bal <- as.character(id) %in% full
  if (!all(bal)) {
    warning(length(counts) - length(full), " individuals are not observed in every period ",
            "and are dropped to balance the panel", call. = FALSE)
    dropped$unbalanced <- length(counts) - length(full)
  }
  id <- id[bal]; time <- time[bal]; first <- first[bal]; y <- y[bal]
  X <- X[bal, , drop = FALSE]
  if (!length(id)) stop("no individuals left after balancing the panel", call. = FALSE)
  ids <- sort(unique(id))
  n <- length(ids)
  row <- match(id, ids)
  col <- match(time, tlist)
  Y <- matrix(NA_real_, n, nT)
  Y[cbind(row, col)] <- y
  Xs <- lapply(seq_len(nT), function(j) {
    sel <- which(col == j)
    out <- X[sel, , drop = FALSE][order(row[sel]), , drop = FALSE]
    dimnames(out) <- list(NULL, colnames(X))
    out
  })
  G <- numeric(n)
  G[row] <- first
  glist <- sort(unique(G[G > 0]))
  if (control == "notyet" && !any(G == 0) && length(glist)) glist <- glist[glist < max(G)]
  glist <- glist[glist > tlist[1] + anticipation]
  if (!length(glist)) stop("no treated cohort with a pre-treatment period", call. = FALSE)
  list(Y = Y, X = Xs, G = G, ids = ids, n = n, periods = tlist, groups = glist,
       outcome = outcome, covariates = setdiff(colnames(X), "(Intercept)"), dropped = dropped)
}

#' Two-group two-period ATT and its influence function (panel data)
#'
#' Outcome regression, normalized IPW and the doubly robust estimator of
#' Sant'Anna and Zhao (2020), as in DRDID::reg_did_panel(),
#' DRDID::std_ipw_did_panel() and DRDID::drdid_panel() with unit weights.
#'
#' @keywords internal
#' @noRd
attgt_2x2 <- function(dy, D, X, method, trim = 0.995) {
  n <- length(D)
  control <- D == 0
  need_or <- method %in% c("dr", "reg")
  need_ps <- method %in% c("dr", "ipw")
  if (need_or) {
    Xc <- X[control, , drop = FALSE]
    qx <- qr(Xc)
    if (qx$rank < ncol(X)) stop("outcome regression design is singular", call. = FALSE)
    beta <- qr.coef(qx, dy[control])
    fitted <- as.numeric(X %*% beta)
    XpX <- crossprod(X[control, , drop = FALSE]) / n
    lin_or <- ((1 - D) * (dy - fitted) * X) %*% solve(XpX)
  }
  if (need_ps) {
    ps_fit <- suppressWarnings(stats::glm.fit(X, D, family = stats::binomial(),
                                              control = list(epsilon = 1e-14, maxit = 100)))
    if (anyNA(ps_fit$coefficients)) stop("propensity score design is singular", call. = FALSE)
    if (!ps_fit$converged) warning("propensity score estimation did not converge", call. = FALSE)
    ps <- pmin(as.numeric(ps_fit$fitted.values), 1 - 1e-6)
    trimmed <- rep(TRUE, n)
    trimmed[control] <- ps[control] < trim
    H <- crossprod(X, ps * (1 - ps) * X)
    lin_ps <- ((D - ps) * X) %*% (chol2inv(chol(H)) * n)
  }
  if (method == "reg") {
    w1 <- D
    eta1 <- mean(w1 * dy) / mean(w1)
    eta0 <- mean(w1 * fitted) / mean(w1)
    M1 <- colSums(w1 * X) / n
    inf1 <- (w1 * dy - w1 * eta1) / mean(w1)
    inf0 <- (w1 * fitted - w1 * eta0 + as.numeric(lin_or %*% M1)) / mean(w1)
    return(list(att = eta1 - eta0, inf = inf1 - inf0))
  }
  w1 <- trimmed * D
  w0 <- trimmed * ps * (1 - D) / (1 - ps)
  if (method == "ipw") {
    eta1 <- mean(w1 * dy) / mean(w1)
    eta0 <- mean(w0 * dy) / mean(w0)
    M2 <- colSums(w0 * (dy - eta0) * X) / n
    inf1 <- (w1 * dy - w1 * eta1) / mean(w1)
    inf0 <- (w0 * dy - w0 * eta0 + as.numeric(lin_ps %*% M2)) / mean(w0)
    return(list(att = eta1 - eta0, inf = inf1 - inf0))
  }
  r <- dy - fitted
  eta1 <- mean(w1 * r) / mean(w1)
  eta0 <- mean(w0 * r) / mean(w0)
  M1 <- colSums(w1 * X) / n
  M2 <- colSums(w0 * (r - eta0) * X) / n
  M3 <- colSums(w0 * X) / n
  inf1 <- (w1 * r - w1 * eta1 - as.numeric(lin_or %*% M1)) / mean(w1)
  inf0 <- (w0 * r - w0 * eta0 + as.numeric(lin_ps %*% M2) - as.numeric(lin_or %*% M3)) /
    mean(w0)
  list(att = eta1 - eta0, inf = inf1 - inf0)
}

#' Aggregate group-time average treatment effects
#'
#' Combines the `ATT(g, t)` of an [attgt()] fit into summary effects, as in
#' Callaway and Sant'Anna (2021, Section 3) and `did::aggte()`:
#' * `"simple"`: the average of all post-treatment `ATT(g, t)` weighted by
#'   cohort size;
#' * `"dynamic"`: the event study, `ATT(e)` for each event time
#'   `e = t - g`, averaging across cohorts weighted by cohort size; the
#'   overall effect is the mean of the post-treatment `ATT(e)`, `e >= 0`;
#' * `"group"`: the average effect of each cohort over its post-treatment
#'   periods; the overall effect weights cohorts by size;
#' * `"calendar"`: the average effect in each calendar period across the
#'   cohorts treated by then (weighted by size); the overall effect is the
#'   mean over periods.
#'
#' Standard errors come from the influence functions and account for the
#' estimation of the cohort-size weights.
#'
#' @param object an `"panglm_attgt"` fit from [attgt()]
#' @param type `"simple"`, `"dynamic"`, `"group"` or `"calendar"`
#' @param min_e,max_e for `"dynamic"`: event-time window
#' @param balance_e for `"dynamic"`: if not `NULL`, only cohorts observed for
#'   at least `balance_e` periods after treatment are used, so that the
#'   composition of cohorts is the same at every event time up to
#'   `balance_e`
#' @param level confidence level; defaults to that of the fit
#' @param ... unused
#' @return an object of class `"summary.panglm_attgt"`: a list with
#'   `overall` (a one-row data frame with `estimate`, `se`, `lower`,
#'   `upper`), `estimates` (the effects by event time, cohort or period:
#'   column `level` plus the same columns; `NULL` for `"simple"`), `type`,
#'   `pretrend`, `inf` (the influence functions of the overall effect and,
#'   as attribute `"by"`, of the effects by level) and the fit settings.
#' @examples
#' set.seed(2)
#' d <- expand.grid(id = 1:150, t = 1:8)
#' d$first <- rep(c(4, 6, NA), each = 50)[d$id]
#' d$y <- rnorm(150)[d$id] + 0.2 * d$t +
#'   ifelse(!is.na(d$first) & d$t >= d$first, 1, 0) + rnorm(nrow(d))
#' fit <- attgt(y ~ 1, d, index = c("id", "t"), treat = "first")
#' summary(fit, type = "dynamic")
#' summary(fit, type = "calendar")
#' @export
summary.panglm_attgt <- function(object, type = c("simple", "dynamic", "group", "calendar"),
                                 min_e = -Inf, max_e = Inf, balance_e = NULL,
                                 level = object$level, ...) {
  type <- match.arg(type)
  est <- object$estimates
  ok <- !is.na(est$estimate) & colSums(is.na(object$inf)) == 0
  if (!all(ok)) {
    warning("dropping ", sum(!ok), " ATT(g, t) that could not be estimated", call. = FALSE)
  }
  est <- est[ok, , drop = FALSE]
  inf <- object$inf[, ok, drop = FALSE]
  n <- object$n
  Gu <- object$units$group
  att <- est$estimate
  grp <- est$group
  pg_all <- vapply(object$groups, function(g) mean(Gu == g), numeric(1))
  pg <- pg_all[match(grp, object$groups)]
  wif <- function(keep, p, gg) {
    S <- sum(p[keep])
    centered <- vapply(keep, function(k) as.numeric(Gu == gg[k]) - p[k], numeric(n))
    centered <- matrix(centered, nrow = n)
    centered / S - rowSums(centered) %*% t(p[keep] / S^2)
  }
  combine <- function(keep, w, with_wif = TRUE, a = att, phi = inf, p = pg, gg = grp) {
    estimate <- sum(w * a[keep])
    f <- as.numeric(phi[, keep, drop = FALSE] %*% w)
    if (with_wif) f <- f + as.numeric(wif(keep, p, gg) %*% a[keep])
    list(estimate = estimate, inf = f)
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  se_of <- function(f) sqrt(sum(f^2)) / n
  by <- NULL
  if (type == "simple") {
    keep <- which(est$post)
    overall <- combine(keep, pg[keep] / sum(pg[keep]))
  } else {
    if (type == "dynamic") {
      use <- rep(TRUE, nrow(est))
      levels_e <- sort(unique(est$event))
      if (!is.null(balance_e)) {
        use <- max(object$periods) - grp >= balance_e
        levels_e <- sort(unique(est$event[use]))
        levels_e <- levels_e[levels_e <= balance_e &
                               levels_e >= balance_e - max(object$periods) + min(object$periods)]
      }
      levels_e <- levels_e[levels_e >= min_e & levels_e <= max_e]
      if (!length(levels_e)) stop("no event time in the requested window", call. = FALSE)
      by <- lapply(levels_e, function(e) {
        keep <- which(est$event == e & use)
        combine(keep, pg[keep] / sum(pg[keep]))
      })
      lev <- levels_e
      final <- lev >= 0
      final_wif <- FALSE
    } else if (type == "group") {
      lev <- object$groups[object$groups %in% grp[est$post]]
      by <- lapply(lev, function(g) {
        keep <- which(grp == g & est$post)
        combine(keep, rep(1 / length(keep), length(keep)), with_wif = FALSE)
      })
      final <- rep(TRUE, length(lev))
      final_wif <- TRUE
    } else {
      lev <- sort(unique(est$time[est$post]))
      lev <- lev[lev >= min(grp)]
      by <- lapply(lev, function(t1) {
        keep <- which(est$time == t1 & est$post)
        combine(keep, pg[keep] / sum(pg[keep]))
      })
      final <- rep(TRUE, length(lev))
      final_wif <- FALSE
    }
    a_by <- vapply(by, `[[`, numeric(1), "estimate")
    inf_by <- matrix(vapply(by, `[[`, numeric(n), "inf"), nrow = n)
    k <- which(final)
    if (!length(k)) stop("no post-treatment effect to aggregate", call. = FALSE)
    if (final_wif) {
      p <- pg_all[match(lev, object$groups)]
      overall <- combine(k, p[k] / sum(p[k]), with_wif = TRUE, a = a_by, phi = inf_by,
                         p = p, gg = lev)
    } else {
      overall <- combine(k, rep(1 / length(k), length(k)), with_wif = FALSE,
                         a = a_by, phi = inf_by)
    }
    se_by <- apply(inf_by, 2, se_of)
    se_by[se_by <= sqrt(.Machine$double.eps) * 10] <- NA
    by <- data.frame(level = lev, estimate = a_by, se = se_by,
                     lower = a_by - z * se_by, upper = a_by + z * se_by)
    attr(by, "inf") <- inf_by
  }
  se <- se_of(overall$inf)
  out <- list(overall = data.frame(estimate = overall$estimate, se = se,
                                   lower = overall$estimate - z * se,
                                   upper = overall$estimate + z * se),
              estimates = by, type = type, inf = overall$inf, pretrend = object$pretrend,
              control = object$control, method = object$method,
              covariates = object$covariates, n = n, level = level, call = object$call)
  class(out) <- "summary.panglm_attgt"
  out
}

#' @export
print.summary.panglm_attgt <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  label <- c(simple = "overall (simple weighted average)",
             dynamic = "event study (by time since treatment)",
             group = "by treatment cohort", calendar = "by calendar period")
  cat("\nCallaway and Sant'Anna group-time ATT, aggregated ", label[[x$type]], "\n", sep = "")
  cat("Control group:", if (x$control == "never") "never treated" else "not yet treated",
      " Estimator:", if (length(x$covariates)) toupper(x$method) else "unconditional",
      " Individuals:", x$n, "\n\n")
  ov <- x$overall
  name <- c(simple = "ATT", dynamic = "ATT (mean of e >= 0)",
            group = "ATT (cohort-size weighted)", calendar = "ATT (mean over periods)")
  cat(name[[x$type]], ": ", format(ov$estimate, digits = digits), " (SE ",
      format(ov$se, digits = digits), ", ", format(100 * x$level), "% CI ",
      format(ov$lower, digits = digits), " to ", format(ov$upper, digits = digits), ")\n",
      sep = "")
  if (!is.null(x$estimates)) {
    tab <- x$estimates
    names(tab)[1] <- c(dynamic = "event", group = "group", calendar = "time")[[x$type]]
    cat("\n")
    print(format(tab, digits = digits), row.names = FALSE)
  }
  if (!is.null(x$pretrend)) {
    cat("\nPre-trend test (pre-treatment ATT(g, t) = 0): Chi-squared = ",
        format(x$pretrend$statistic, digits = digits), ", df = ", x$pretrend$df,
        ", p-value = ", format.pval(x$pretrend$p.value, digits = digits), "\n", sep = "")
  }
  cat("\n")
  invisible(x)
}

#' @export
print.panglm_attgt <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nGroup-time average treatment effects (Callaway and Sant'Anna 2021)\n")
  cat("Control group:", if (x$control == "never") "never treated" else "not yet treated",
      " Estimator:", if (length(x$covariates)) toupper(x$method) else "unconditional",
      " Base period:", x$base, "\n")
  cat("Individuals:", x$n, " Periods:", length(x$periods), " Cohorts:",
      paste(x$groups, collapse = ", "), "\n\n")
  tab <- x$estimates
  tab$post <- NULL
  print(format(tab, digits = digits), row.names = FALSE)
  if (!is.null(x$pretrend)) {
    cat("\nPre-trend test: Chi-squared = ", format(x$pretrend$statistic, digits = digits),
        ", df = ", x$pretrend$df, ", p-value = ",
        format.pval(x$pretrend$p.value, digits = digits), "\n", sep = "")
  }
  cat("Use summary(x, type = ...) for the simple, dynamic, group or calendar aggregation.\n\n")
  invisible(x)
}

#' @export
coef.panglm_attgt <- function(object, ...) {
  stats::setNames(object$estimates$estimate, colnames(object$inf))
}

#' @export
vcov.panglm_attgt <- function(object, ...) {
  crossprod(object$inf) / object$n^2
}

#' @export
nobs.panglm_attgt <- function(object, ...) object$n

#' Plot an event study from group-time effects
#'
#' Plots the dynamic aggregation of [attgt()] (`ATT(e)` by event time with
#' pointwise confidence intervals), or the cohort or calendar aggregation.
#'
#' @param x an `"panglm_attgt"` fit from [attgt()], or its summary
#' @param type aggregation to plot (for a fit): `"dynamic"`, `"group"` or
#'   `"calendar"`
#' @param xlab,ylab axis labels
#' @param ... passed to [summary.panglm_attgt()] (for example `min_e`,
#'   `max_e`) and then to [graphics::plot()]
#' @return the plotted summary, invisibly
#' @examples
#' set.seed(2)
#' d <- expand.grid(id = 1:150, t = 1:8)
#' d$first <- rep(c(4, 6, NA), each = 50)[d$id]
#' d$y <- rnorm(150)[d$id] + 0.2 * d$t +
#'   ifelse(!is.na(d$first) & d$t >= d$first, 1, 0) + rnorm(nrow(d))
#' plot(attgt(y ~ 1, d, index = c("id", "t"), treat = "first"))
#' @export
plot.panglm_attgt <- function(x, type = c("dynamic", "group", "calendar"), xlab = NULL,
                              ylab = "ATT", ...) {
  type <- match.arg(type)
  dots <- list(...)
  sargs <- intersect(names(dots), c("min_e", "max_e", "balance_e", "level"))
  s <- do.call(summary.panglm_attgt, c(list(x, type = type), dots[sargs]))
  do.call(plot.summary.panglm_attgt, c(list(s, xlab = xlab, ylab = ylab),
                                       dots[setdiff(names(dots), sargs)]))
}

#' @rdname plot.panglm_attgt
#' @export
plot.summary.panglm_attgt <- function(x, type = NULL, xlab = NULL, ylab = "ATT", ...) {
  tab <- x$estimates
  if (is.null(tab)) stop("nothing to plot for the simple aggregation", call. = FALSE)
  if (is.null(xlab)) {
    xlab <- c(dynamic = "Event time", group = "Treatment cohort",
              calendar = "Period")[[x$type]]
  }
  ok <- is.finite(tab$estimate)
  lo <- ifelse(is.finite(tab$lower), tab$lower, tab$estimate)
  hi <- ifelse(is.finite(tab$upper), tab$upper, tab$estimate)
  pre <- x$type == "dynamic" & tab$level < 0
  graphics::plot(tab$level[ok], tab$estimate[ok], pch = ifelse(pre[ok], 1, 19),
                 xlab = xlab, ylab = ylab, ylim = range(lo[ok], hi[ok], 0), xaxt = "n", ...)
  graphics::axis(1, at = tab$level)
  graphics::segments(tab$level[ok], lo[ok], tab$level[ok], hi[ok])
  graphics::abline(h = 0, lty = 3)
  if (x$type == "dynamic") graphics::abline(v = -0.5, lty = 2)
  invisible(x)
}
