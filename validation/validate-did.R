# Heterogeneity-robust difference-in-differences.
# attgt() against did::att_gt() and did::aggte() (Callaway and Sant'Anna
# 2021): ATT(g,t), analytic standard errors, the pre-trend Wald statistic
# and the simple, dynamic, group and calendar aggregations, for never- and
# not-yet-treated controls, the dr, ipw and reg estimators, varying and
# universal base periods, with and without a covariate; plus a sample
# without never-treated units, anticipation, an unbalanced panel and a
# balanced event window, and a subset whose factor individual index keeps
# unused levels.
# etwfe() against etwfe::etwfe() and etwfe::emfx() (Wooldridge 2021, 2023):
# simple, group, calendar and event ATTs for linear, Poisson and logit
# models, not-yet and never-treated controls, with and without a control.
# Standard errors are clustered by county; panglm's small-sample factor
# G/(G-1)*(n-1)/(n-k) matches fixest's default. emfx() differentiates
# numerically (marginaleffects), so Poisson ATTs are compared with relative
# tolerance 1e-6 and their standard errors with 1e-4; the event time -1
# reference of the never-treated event study (zero, no standard error in
# emfx) is left out of the standard-error comparison.
# Data: did::mpdta (500 US counties, 2003-2007, minimum wage, log teen
# employment lemp, log population lpop). Counts: round(exp(lemp)); binary:
# lemp above its median.
source("validation/common.R")
suppressPackageStartupMessages(library(panglm))
data(mpdta, package = "did")
mpdta$emp <- round(exp(mpdta$lemp))
mpdta$hi <- as.integer(mpdta$lemp > stats::median(mpdta$lemp))
idx <- c("countyreal", "year")
out <- list()
add <- function(...) out[[length(out) + 1L]] <<- vcompare(...)

## attgt ---------------------------------------------------------------------
att_check <- function(label, data, xf, control, method, base = "varying", anticipation = 0) {
  ref <- did::att_gt("lemp", tname = "year", idname = "countyreal", gname = "first.treat",
                     xformla = xf, data = data, est_method = method,
                     control_group = if (control == "never") "nevertreated" else "notyettreated",
                     base_period = base, anticipation = anticipation, bstrap = FALSE,
                     cband = FALSE, allow_unbalanced_panel = FALSE)
  f <- if (is.null(xf)) lemp ~ 1 else lemp ~ lpop
  fit <- suppressWarnings(attgt(f, data, idx, "first.treat", control = control, method = method,
                                base = base, anticipation = anticipation))
  mod <- sprintf("attgt %s %s %s%s", control, method, base,
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
att_check("mpdta, anticipation 1", mpdta, ~lpop, "notyet", "dr", anticipation = 1)
set.seed(20261025)
unb <- mpdta[-sample(nrow(mpdta), 50), ]
for (cg in c("never", "notyet")) att_check("mpdta, unbalanced", unb, ~lpop, cg, "dr")
ref <- did::att_gt("lemp", tname = "year", idname = "countyreal", gname = "first.treat",
                   xformla = ~lpop, data = mpdta, bstrap = FALSE, cband = FALSE)
ag <- did::aggte(ref, type = "dynamic", balance_e = 1, min_e = -2, max_e = 2,
                 bstrap = FALSE, cband = FALSE)
s <- summary(attgt(lemp ~ lpop, mpdta, idx, "first.treat"), type = "dynamic",
             balance_e = 1, min_e = -2, max_e = 2)
add("attgt never dr varying + lpop", "mpdta, balance_e 1, e in [-2, 2]", "dynamic ATT",
    c(s$overall$estimate, s$estimates$estimate), list(did = c(ag$overall.att, ag$att.egt)), 1e-8)
add("attgt never dr varying + lpop", "mpdta, balance_e 1, e in [-2, 2]", "dynamic SE",
    c(s$overall$se, s$estimates$se), list(did = c(ag$overall.se, ag$se.egt)), 1e-6)

## etwfe ---------------------------------------------------------------------
for (fam in c("gaussian", "poisson", "logit")) for (cg in c("notyet", "never")) {
  for (ctrl in c(FALSE, TRUE)) {
    yv <- c(gaussian = "lemp", poisson = "emp", logit = "hi")[[fam]]
    f <- stats::as.formula(paste(yv, if (ctrl) "~ lpop" else "~ 1"))
    ref <- suppressMessages(etwfe::etwfe(f, tvar = year, gvar = first.treat, data = mpdta,
                                         cgroup = cg, vcov = ~countyreal,
                                         family = if (fam == "gaussian") NULL else fam))
    fit <- etwfe(f, mpdta, idx, "first.treat", family = fam, control = cg)
    mod <- sprintf("etwfe %s %s%s", fam, cg, if (ctrl) " + lpop" else "")
    tol_b <- if (fam == "poisson") 1e-6 else 1e-8
    tol_se <- if (fam == "poisson") 1e-4 else 1e-6
    for (type in c("simple", "group", "calendar", "event")) {
      e <- etwfe::emfx(ref, type = type)
      s <- summary(fit, type = type)
      tab <- if (type == "simple") s$overall else s$estimates
      add(mod, "mpdta", paste(type, "ATT"), tab$estimate, list(etwfe = e$estimate), tol_b)
      ok <- !is.na(e$std.error)
      add(mod, "mpdta", paste(type, "SE"), tab$se[ok], list(etwfe = e$std.error[ok]), tol_se)
    }
  }
}
for (vc in c("cluster", "classical")) {
  d <- if (vc == "cluster") nonever else mpdta
  ref <- suppressMessages(etwfe::etwfe(lemp ~ lpop, tvar = year, gvar = first.treat, data = d,
                                       vcov = if (vc == "cluster") ~countyreal else "iid"))
  fit <- suppressWarnings(etwfe(lemp ~ lpop, d, idx, "first.treat", vcov = vc))
  lab <- if (vc == "cluster") "mpdta, no never-treated" else "mpdta, classical vcov"
  e <- etwfe::emfx(ref, type = "event")
  s <- summary(fit, type = "event")$estimates
  add("etwfe gaussian notyet + lpop", lab, "event ATT", s$estimate, list(etwfe = e$estimate), 1e-8)
  add("etwfe gaussian notyet + lpop", lab, "event SE", s$se, list(etwfe = e$std.error), 1e-6)
}
e <- etwfe::emfx(suppressMessages(etwfe::etwfe(lemp ~ lpop, tvar = year, gvar = first.treat,
                                               data = unb, vcov = ~countyreal)), type = "group")
s <- summary(etwfe(lemp ~ lpop, unb, idx, "first.treat"), type = "group")$estimates
add("etwfe gaussian notyet + lpop", "mpdta, unbalanced", "group ATT", s$estimate,
    list(etwfe = e$estimate), 1e-8)
add("etwfe gaussian notyet + lpop", "mpdta, unbalanced", "group SE", s$se,
    list(etwfe = e$std.error), 1e-6)

## factor index with unused levels (a subset of a factor-indexed panel) ----
sub <- mpdta[mpdta$countyreal %in% sort(unique(mpdta$countyreal))[-(1:100)], ]
subf <- sub
subf$countyreal <- factor(subf$countyreal, levels = sort(unique(mpdta$countyreal)))
ref <- did::att_gt("lemp", tname = "year", idname = "countyreal", gname = "first.treat",
                   data = sub, bstrap = FALSE, cband = FALSE)
fit <- attgt(lemp ~ 1, subf, idx, "first.treat")
add("attgt never dr varying", "mpdta subset, factor id with unused levels", "ATT(g,t)",
    fit$estimates$estimate, list(did = ref$att), 1e-8)
add("attgt never dr varying", "mpdta subset, factor id with unused levels", "SE ATT(g,t)",
    fit$estimates$se, list(did = ref$se), 1e-6)

do.call(vsave, c(out, name = "did"))
