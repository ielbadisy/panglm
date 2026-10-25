# Heterogeneity-robust difference-in-differences with sampling weights.
# attgt(weights = ) against did::att_gt(weightsname = ) and did::aggte():
# ATT(g,t), analytic standard errors, the pre-trend Wald statistic and the
# simple, dynamic, group and calendar aggregations, for never- and
# not-yet-treated controls, the dr, ipw and reg estimators, varying and
# universal base periods, with and without a covariate; plus a sample
# without never-treated units and an unbalanced panel.
# etwfe(weights = ) against etwfe::etwfe(weights = ~w) (fixest weights) and
# etwfe::emfx(wts = "w"): simple, group, calendar and event ATTs for linear,
# Poisson and logit models, not-yet and never-treated controls, with and
# without a control. The nonlinear references are fitted with
# glm.tol = 1e-12 (fixest stops at a deviance tolerance of 1e-8 by default
# and evaluates the Hessian at its last iterate, which moves the weighted
# standard errors in the sixth decimal). emfx() differentiates numerically
# (marginaleffects), so the nonlinear ATTs are compared with relative
# tolerance 1e-6 and their standard errors with 1e-4, as the Poisson ones in
# validate-did.R. The logit outcome is a noisy indicator of high employment
# (the deterministic indicator of validate-did.R is quasi-separated by lpop
# within the 2004 cohort, so its logit MLE does not exist and the
# comparison would depend on where each optimizer stops).
# Data: did::mpdta (500 US counties, 2003-2007). County-level weights are
# drawn once (uniform on [0.2, 3]) and are constant within county.
source("validation/common.R")
suppressPackageStartupMessages(library(panglm))
data(mpdta, package = "did")
set.seed(20261025)
wc <- stats::setNames(stats::runif(500, 0.2, 3), unique(mpdta$countyreal))
mpdta$w <- unname(wc[as.character(mpdta$countyreal)])
mpdta$emp <- round(exp(mpdta$lemp))
mpdta$hiw <- stats::rbinom(nrow(mpdta), 1,
                           stats::plogis(2 * (mpdta$lemp - stats::median(mpdta$lemp))))
idx <- c("countyreal", "year")
out <- list()
add <- function(...) out[[length(out) + 1L]] <<- vcompare(...)

## attgt ---------------------------------------------------------------------
att_check <- function(label, data, xf, control, method, base = "varying") {
  ref <- did::att_gt("lemp", tname = "year", idname = "countyreal", gname = "first.treat",
                     xformla = xf, data = data, est_method = method, weightsname = "w",
                     control_group = if (control == "never") "nevertreated" else "notyettreated",
                     base_period = base, bstrap = FALSE, cband = FALSE,
                     allow_unbalanced_panel = FALSE)
  f <- if (is.null(xf)) lemp ~ 1 else lemp ~ lpop
  fit <- suppressWarnings(attgt(f, data, idx, "first.treat", control = control, method = method,
                                base = base, weights = "w"))
  mod <- sprintf("attgt weighted %s %s %s%s", control, method, base,
                 if (is.null(xf)) "" else " + lpop")
  add(mod, label, "ATT(g,t)", fit$estimates$estimate, list(did = ref$att), 1e-8)
  ok <- !is.na(ref$se)
  add(mod, label, "SE ATT(g,t)", fit$estimates$se[ok], list(did = ref$se[ok]), 1e-6)
  if (!is.null(ref$W)) add(mod, label, "pre-trend Wald", fit$pretrend$statistic, list(did = ref$W), 1e-6)
  for (type in c("simple", "dynamic", "group", "calendar")) {
    ag <- did::aggte(ref, type = type, bstrap = FALSE, cband = FALSE)
    s <- summary(fit, type = type)
    add(mod, label, paste(type, "overall"), s$overall$estimate, list(did = ag$overall.att), 1e-8)
    add(mod, label, paste(type, "overall SE"), s$overall$se, list(did = ag$overall.se), 1e-6)
    if (type != "simple") {
      add(mod, label, paste(type, "ATT"), s$estimates$estimate, list(did = ag$att.egt), 1e-8)
      ok <- !is.na(ag$se.egt)
      add(mod, label, paste(type, "SE"), s$estimates$se[ok], list(did = ag$se.egt[ok]), 1e-6)
    }
  }
}
for (cg in c("never", "notyet")) for (m in c("dr", "ipw", "reg")) for (b in c("varying", "universal")) {
  for (xf in list(NULL, ~lpop)) att_check("mpdta", mpdta, xf, cg, m, b)
}
nonever <- subset(mpdta, first.treat != 0)
for (cg in c("never", "notyet")) att_check("mpdta, no never-treated", nonever, ~lpop, cg, "dr")
unb <- mpdta[-sample(nrow(mpdta), 50), ]
for (cg in c("never", "notyet")) att_check("mpdta, unbalanced", unb, ~lpop, cg, "ipw")

## etwfe ---------------------------------------------------------------------
for (fam in c("gaussian", "poisson", "logit")) for (cg in c("notyet", "never")) {
  for (ctrl in c(FALSE, TRUE)) {
    yv <- c(gaussian = "lemp", poisson = "emp", logit = "hiw")[[fam]]
    f <- stats::as.formula(paste(yv, if (ctrl) "~ lpop" else "~ 1"))
    args <- list(f, tvar = quote(year), gvar = quote(first.treat), data = mpdta, cgroup = cg,
                 vcov = ~countyreal, weights = ~w)
    if (fam != "gaussian") args <- c(args, list(family = fam, glm.tol = 1e-12))
    ref <- suppressMessages(do.call(etwfe::etwfe, args))
    fit <- etwfe(f, mpdta, idx, "first.treat", family = fam, control = cg, weights = "w")
    mod <- sprintf("etwfe weighted %s %s%s", fam, cg, if (ctrl) " + lpop" else "")
    tol_b <- if (fam == "gaussian") 1e-8 else 1e-6
    tol_se <- if (fam == "gaussian") 1e-6 else 1e-4
    for (type in c("simple", "group", "calendar", "event")) {
      e <- etwfe::emfx(ref, type = type, wts = "w")
      s <- summary(fit, type = type)
      tab <- if (type == "simple") s$overall else s$estimates
      add(mod, "mpdta", paste(type, "ATT"), tab$estimate, list(etwfe = e$estimate), tol_b)
      ok <- !is.na(e$std.error)
      add(mod, "mpdta", paste(type, "SE"), tab$se[ok], list(etwfe = e$std.error[ok]), tol_se)
    }
  }
}
ref <- suppressMessages(etwfe::etwfe(lemp ~ lpop, tvar = year, gvar = first.treat, data = mpdta,
                                     vcov = "iid", weights = ~w))
fit <- etwfe(lemp ~ lpop, mpdta, idx, "first.treat", vcov = "classical", weights = "w")
e <- etwfe::emfx(ref, type = "event", wts = "w")
s <- summary(fit, type = "event")$estimates
add("etwfe weighted gaussian notyet + lpop", "mpdta, classical vcov", "event ATT", s$estimate,
    list(etwfe = e$estimate), 1e-8)
add("etwfe weighted gaussian notyet + lpop", "mpdta, classical vcov", "event SE", s$se,
    list(etwfe = e$std.error), 1e-6)

do.call(vsave, c(out, name = "did-weights"))
