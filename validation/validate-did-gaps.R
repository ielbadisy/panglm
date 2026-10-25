# Difference-in-differences: binary outcomes in did() and eventstudy().
#
# did(family = "binomial") and eventstudy(family = "binomial") fit the
# two-way conditional (Chamberlain) logit: individual effects conditioned
# out, time effects as dummies inside the conditional likelihood.
# References:
#   - survival::clogit(y ~ treated + controls + factor(time) + strata(id),
#     method = "exact"), the same estimator: coefficients (log odds ratios)
#     and classical standard errors; for eventstudy() the event-time
#     dummies are coded by hand with -1 (and never-treated units) as the
#     reference.
#   - Cluster-robust standard errors: clogit has no robust covariance with
#     method = "exact", so the stratum scores are computed by brute force
#     (enumerating every 0/1 sequence with the observed number of
#     successes, combn) at clogit's estimate and combined with clogit's
#     bread over the full parameter vector (treatment, controls, time
#     dummies), with the factor G / (G - 1) (m - 1) / (m - K) (m informative
#     strata, K parameters), as in validate-twowaylogit.R. The pre-trend
#     Wald statistic is recomputed from that covariance.
#   - fixest::feglm(family = binomial) with individual and time fixed
#     effects is the unconditional logit, a different estimator (incidental
#     parameters bias). With two periods, however, its estimate is exactly
#     twice the conditional logit estimate (Abrevaya 1997), which gives an
#     exact check of the 2 x 2 design; the cluster standard error of that
#     design is checked against the logit of the second-period outcome on
#     the change in treatment among switchers (glm + sandwich::vcovCL).
#
# etwfe(): heterogeneity by a covariate (xvar) and individual fixed effects
# (ivar), against etwfe::etwfe(xvar = , ivar = ) and etwfe::emfx() (with
# by_xvar = TRUE by default for an xvar fit): simple, group, calendar and
# event ATTs and their standard errors, by value of the covariate, for the
# linear, Poisson and logit models, not-yet and never-treated controls,
# with and without a control, without and with weights (weights = ~w and
# emfx(wts = "w")). ivar is for the linear model only (etwfe ignores it for
# nonlinear families): cluster and classical (iid) covariances, with and
# without weights and xvar, an unbalanced panel and a sample without
# never-treated units. Weights that vary within individuals (row-level
# weights, which attgt() refuses) for every family, and with ivar and xvar
# for the linear model. For the linear model emfx() is called with
# numderiv = "richardson": its default forward differences move some
# weighted standard errors in the sixth decimal, while Richardson
# extrapolation agrees with panglm's analytic gradients (it is too slow for
# the nonlinear models, compared with relative tolerance 1e-4 as in
# validate-did.R). Tolerances and the nonlinear glm.tol = 1e-12 are
# otherwise those of validate-did.R and validate-did-weights.R.
# Data: did::mpdta (500 counties, 2003-2007; binary outcome a noisy
# indicator of high teen employment, as in validate-did-weights.R) and a
# simulated staggered panel with a time-varying control (300 individuals,
# 8 periods, cohorts 2014 and 2016 and never treated). For etwfe, mpdta
# with counts round(exp(lemp)), county weights uniform on [0.2, 3] and the
# covariate big = 1(log population above its median).
source("validation/common.R")
suppressPackageStartupMessages({
  library(survival)
  library(sandwich)
})
data(mpdta, package = "did")
set.seed(20261025)
mpdta$hiw <- stats::rbinom(nrow(mpdta), 1,
                           stats::plogis(2 * (mpdta$lemp - stats::median(mpdta$lemp))))
idx <- c("countyreal", "year")
se <- function(v) sqrt(diag(v))
out <- list()
add <- function(...) out[[length(out) + 1L]] <<- vcompare(...)

sim_did <- function(N = 300, seed = 6) {
  set.seed(seed)
  d <- expand.grid(id = seq_len(N), year = 2011:2018)
  d$first <- sample(c(2014, 2016, NA), N, replace = TRUE)[d$id]
  a <- stats::rnorm(N)
  d$x <- stats::rnorm(nrow(d))
  post <- !is.na(d$first) & d$year >= d$first
  d$yb <- stats::rbinom(nrow(d), 1, stats::plogis(a[d$id] - 0.2 + 0.1 * (d$year - 2011) +
                                                    0.4 * d$x + 0.7 * post))
  d
}
sim <- sim_did()

## brute-force stratum scores (independent of panglm's recursion) ------------
brute_scores <- function(Z, y, id, theta) {
  eta <- as.numeric(Z %*% theta)
  groups <- split(seq_along(y), id)
  t(vapply(groups, function(rows) {
    S <- sum(y[rows]); Tn <- length(rows)
    if (S == 0 || S == Tn) return(rep(NA_real_, ncol(Z)))
    sets <- utils::combn(Tn, S, simplify = FALSE)
    lw <- vapply(sets, function(s) sum(eta[rows][s]), 0)
    w <- exp(lw - max(lw)); w <- w / sum(w)
    ex <- Reduce(`+`, Map(function(s, wi) wi * colSums(Z[rows[s], , drop = FALSE]), sets, w))
    colSums(Z[rows[y[rows] == 1], , drop = FALSE]) - ex
  }, numeric(ncol(Z))))
}
brute_cluster <- function(cl, Z, y, id) {
  U <- brute_scores(Z, y, id, coef(cl))
  U <- U[stats::complete.cases(U), , drop = FALSE]
  m <- nrow(U); K <- ncol(U)
  B <- vcov(cl)
  V <- B %*% crossprod(U) %*% B * (m / (m - 1)) * ((m - 1) / (m - K))
  dimnames(V) <- dimnames(B)
  V
}
clogit_fit <- function(f, data) {
  suppressWarnings(clogit(f, data, method = "exact",
                          control = coxph.control(eps = 1e-12, toler.chol = 1e-14,
                                                  iter.max = 100)))
}

## did() -----------------------------------------------------------------------
did_check <- function(label, data, index, treat, y, xnames) {
  tv <- data[[treat]]
  data$tr <- as.numeric(!is.na(tv) & tv > 0 & data[[index[2]]] >= tv)
  rhs <- paste(c("tr", xnames), collapse = " + ")
  cl <- clogit_fit(stats::as.formula(sprintf("%s ~ %s + factor(%s) + strata(%s)", y, rhs,
                                             index[2], index[1])), data)
  Z <- stats::model.matrix(stats::as.formula(sprintf("~ %s + factor(%s)", rhs, index[2])),
                           data)[, -1, drop = FALSE]
  Z <- Z[, names(coef(cl)), drop = FALSE]
  Vb <- brute_cluster(cl, Z, data[[y]], data[[index[1]]])
  b <- c("tr", xnames)
  f <- stats::as.formula(paste(y, "~", if (length(xnames)) paste(xnames, collapse = " + ") else "1"))
  fc <- suppressMessages(did(f, data, index, treat, family = "binomial", vcov = "classical"))
  fr <- suppressMessages(did(f, data, index, treat, family = "binomial"))
  mod <- "did binomial (conditional logit)"
  add(mod, label, "coef", coef(fr), list(survival = coef(cl)[b]), 1e-6)
  add(mod, label, "SE classical", se(vcov(fc)), list(survival = se(vcov(cl))[b]), 1e-6)
  add(mod, label, "SE cluster", se(vcov(fr)), list(brute_force = se(Vb)[b]), 1e-6)
}
did_check("mpdta", mpdta, idx, "first.treat", "hiw", character(0))
did_check("simulated, control x", sim, c("id", "year"), "first", "yb", "x")

## eventstudy() ----------------------------------------------------------------
es_check <- function(label, data, index, treat, y, xnames, leads, lags) {
  tv <- data[[treat]]
  tv[!is.na(tv) & tv == 0] <- NA
  rel <- pmin(pmax(data[[index[2]]] - tv, -leads), lags)
  rel[is.na(rel)] <- -1
  data$rel <- stats::relevel(factor(rel), "-1")
  rhs <- paste(c("rel", xnames), collapse = " + ")
  cl <- clogit_fit(stats::as.formula(sprintf("%s ~ %s + factor(%s) + strata(%s)", y, rhs,
                                             index[2], index[1])), data)
  Z <- stats::model.matrix(stats::as.formula(sprintf("~ %s + factor(%s)", rhs, index[2])),
                           data)[, -1, drop = FALSE]
  Z <- Z[, names(coef(cl)), drop = FALSE]
  Vb <- brute_cluster(cl, Z, data[[y]], data[[index[1]]])
  ev <- setdiff(seq(-leads, lags), -1)
  b <- paste0("rel", ev)
  f <- stats::as.formula(paste(y, "~", if (length(xnames)) paste(xnames, collapse = " + ") else "1"))
  ec <- suppressMessages(eventstudy(f, data, index, treat, leads = leads, lags = lags,
                                    family = "binomial", vcov = "classical"))
  er <- suppressMessages(eventstudy(f, data, index, treat, leads = leads, lags = lags,
                                    family = "binomial"))
  tc <- ec$estimates[ec$estimates$event != -1, ]
  tr <- er$estimates[er$estimates$event != -1, ]
  mod <- "eventstudy binomial (conditional logit)"
  add(mod, label, "event coef", tr$estimate, list(survival = coef(cl)[b]), 1e-6)
  add(mod, label, "event SE classical", tc$se, list(survival = se(vcov(cl))[b]), 1e-6)
  add(mod, label, "event SE cluster", tr$se, list(brute_force = se(Vb)[b]), 1e-6)
  lead <- paste0("rel", ev[ev < 0])
  W <- as.numeric(t(coef(cl)[lead]) %*% solve(Vb[lead, lead], coef(cl)[lead]))
  add(mod, label, "pre-trend Wald (cluster)", er$pretrend$statistic, list(brute_force = W), 1e-6)
}
es_check("mpdta", mpdta, idx, "first.treat", "hiw", character(0), leads = 3, lags = 3)
es_check("simulated, control x", sim, c("id", "year"), "first", "yb", "x", leads = 3, lags = 2)

## 2 x 2 design: fixest (Abrevaya 1997) and the differenced logit -------------
two <- subset(sim, year %in% c(2013, 2014) & (is.na(first) | first == 2014))
two$tr <- as.numeric(!is.na(two$first) & two$year == 2014)
f2 <- suppressMessages(did(yb ~ x, two, c("id", "year"), "first", family = "binomial"))
fx <- fixest::feglm(yb ~ tr + x | id + year, two, family = binomial, glm.tol = 1e-12,
                    notes = FALSE)
y1 <- two[two$year == 2013, ]; y2 <- two[two$year == 2014, ]
y1 <- y1[order(y1$id), ]; y2 <- y2[order(y2$id), ]
stopifnot(identical(y1$id, y2$id))
fd <- data.frame(y = y2$yb, dtr = y2$tr - y1$tr, dx = y2$x - y1$x)[y1$yb != y2$yb, ]
gd <- stats::glm(y ~ dtr + dx, binomial, fd,
                 control = stats::glm.control(epsilon = 1e-14, maxit = 100))
mod <- "did binomial (conditional logit), 2 x 2"
lab <- "simulated, 2013-2014, cohort 2014 vs never, control x"
add(mod, lab, "coef", coef(f2), list(fixest_half = coef(fx)[c("tr", "x")] / 2,
                                     glm_fd = coef(gd)[c("dtr", "dx")]), 1e-6)
add(mod, lab, "SE cluster", se(vcov(f2)),
    list(sandwich = se(vcovCL(gd, cluster = seq_len(nrow(fd)), type = "HC1"))[c("dtr", "dx")]),
    1e-6)

## etwfe(): xvar and ivar ------------------------------------------------------
mpdta$emp <- round(exp(mpdta$lemp))
mpdta$big <- as.integer(mpdta$lpop > stats::median(mpdta$lpop))
wc <- stats::setNames(stats::runif(500, 0.2, 3), unique(mpdta$countyreal))
mpdta$w <- unname(wc[as.character(mpdta$countyreal)])
et_check <- function(mod, label, ref, fit, weighted, tol_b, tol_se) {
  for (type in c("simple", "group", "calendar", "event")) {
    args <- list(ref, type = type)
    if (weighted) args$wts <- "w"
    if (is.null(ref$family)) args$numderiv <- "richardson"
    e <- do.call(etwfe::emfx, args)
    s <- summary(fit, type = type)
    tab <- if (is.null(s$estimates)) s$overall else s$estimates
    add(mod, label, paste(type, "ATT"), tab$estimate, list(etwfe = e$estimate), tol_b)
    ok <- !is.na(e$std.error)
    add(mod, label, paste(type, "SE"), tab$se[ok], list(etwfe = e$std.error[ok]), tol_se)
  }
}
et_ref <- function(f, data, fam, cg, weighted, xvar = NULL, ivar = NULL, vcov = ~countyreal) {
  args <- list(f, tvar = quote(year), gvar = quote(first.treat), data = data, cgroup = cg,
               vcov = vcov)
  if (!is.null(xvar)) args$xvar <- as.name(xvar)
  if (!is.null(ivar)) args$ivar <- as.name(ivar)
  if (weighted) args$weights <- ~w
  if (fam != "gaussian") args <- c(args, list(family = fam, glm.tol = 1e-12))
  suppressMessages(do.call(etwfe::etwfe, args))
}
for (fam in c("gaussian", "poisson", "logit")) for (cg in c("notyet", "never")) {
  for (ctrl in c(FALSE, TRUE)) for (wtd in c(FALSE, TRUE)) {
    yv <- c(gaussian = "lemp", poisson = "emp", logit = "hiw")[[fam]]
    f <- stats::as.formula(paste(yv, if (ctrl) "~ lpop" else "~ 1"))
    ref <- et_ref(f, mpdta, fam, cg, wtd, xvar = "big")
    fit <- etwfe(f, mpdta, idx, "first.treat", family = fam, control = cg, xvar = "big",
                 weights = if (wtd) "w")
    mod <- sprintf("etwfe xvar %s %s%s%s", fam, cg, if (ctrl) " + lpop" else "",
                   if (wtd) ", weighted" else "")
    et_check(mod, "mpdta, by big", ref, fit, wtd, if (fam == "gaussian") 1e-8 else 1e-6,
             if (fam == "gaussian") 1e-6 else 1e-4)
  }
}
nonever <- subset(mpdta, first.treat != 0)
unb <- mpdta[-sample(nrow(mpdta), 50), ]
for (cg in c("notyet", "never")) for (ctrl in c(FALSE, TRUE)) for (wtd in c(FALSE, TRUE)) {
  for (vc in c("cluster", "classical")) {
    f <- if (ctrl) lemp ~ lpop else lemp ~ 1
    ref <- et_ref(f, mpdta, "gaussian", cg, wtd, ivar = "countyreal",
                  vcov = if (vc == "cluster") ~countyreal else "iid")
    fit <- etwfe(f, mpdta, idx, "first.treat", control = cg, ivar = TRUE, vcov = vc,
                 weights = if (wtd) "w")
    mod <- sprintf("etwfe ivar gaussian %s%s%s", cg, if (ctrl) " + lpop" else "",
                   if (wtd) ", weighted" else "")
    et_check(mod, paste0("mpdta, ", vc, " vcov"), ref, fit, wtd, 1e-8, 1e-6)
  }
}
for (wtd in c(FALSE, TRUE)) {
  for (cg in c("notyet", "never")) {
    ref <- et_ref(lemp ~ lpop, mpdta, "gaussian", cg, wtd, xvar = "big", ivar = "countyreal")
    fit <- etwfe(lemp ~ lpop, mpdta, idx, "first.treat", control = cg, ivar = TRUE,
                 xvar = "big", weights = if (wtd) "w")
    et_check(sprintf("etwfe ivar + xvar gaussian %s + lpop%s", cg, if (wtd) ", weighted" else ""),
             "mpdta, by big", ref, fit, wtd, 1e-8, 1e-6)
  }
  ref <- et_ref(lemp ~ lpop, unb, "gaussian", "notyet", wtd, ivar = "countyreal")
  fit <- etwfe(lemp ~ lpop, unb, idx, "first.treat", ivar = TRUE, weights = if (wtd) "w")
  et_check(sprintf("etwfe ivar gaussian notyet + lpop%s", if (wtd) ", weighted" else ""),
           "mpdta, unbalanced", ref, fit, wtd, 1e-8, 1e-6)
  ref <- et_ref(lemp ~ lpop, nonever, "gaussian", "notyet", wtd, ivar = "countyreal")
  fit <- suppressWarnings(etwfe(lemp ~ lpop, nonever, idx, "first.treat", ivar = TRUE,
                                weights = if (wtd) "w"))
  et_check(sprintf("etwfe ivar gaussian notyet + lpop%s", if (wtd) ", weighted" else ""),
           "mpdta, no never-treated", ref, fit, wtd, 1e-8, 1e-6)
}

## etwfe(): weights varying within individuals --------------------------------
mpdta$wr <- stats::runif(nrow(mpdta), 0.2, 3)
mpdta_r <- transform(mpdta, w = wr)
for (fam in c("gaussian", "poisson", "logit")) for (cg in c("notyet", "never")) {
  yv <- c(gaussian = "lemp", poisson = "emp", logit = "hiw")[[fam]]
  f <- stats::as.formula(paste(yv, "~ lpop"))
  ref <- et_ref(f, mpdta_r, fam, cg, TRUE)
  fit <- etwfe(f, mpdta_r, idx, "first.treat", family = fam, control = cg, weights = "w")
  et_check(sprintf("etwfe row weights %s %s + lpop", fam, cg), "mpdta, row-level weights",
           ref, fit, TRUE, if (fam == "gaussian") 1e-8 else 1e-6,
           if (fam == "gaussian") 1e-6 else 1e-4)
}
ref <- et_ref(lemp ~ lpop, mpdta_r, "gaussian", "notyet", TRUE, xvar = "big",
              ivar = "countyreal")
fit <- etwfe(lemp ~ lpop, mpdta_r, idx, "first.treat", ivar = TRUE, xvar = "big",
             weights = "w")
et_check("etwfe row weights ivar + xvar gaussian notyet + lpop", "mpdta, row-level weights",
         ref, fit, TRUE, 1e-8, 1e-6)

do.call(vsave, c(out, name = "did-gaps"))
