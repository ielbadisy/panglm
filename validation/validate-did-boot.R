# Clustered standard errors, multiplier bootstrap and uniform confidence
# bands of attgt() (Callaway and Sant'Anna 2021).
# 1. Clustered analytic standard errors (influence functions summed by
#    state, 29 states) against did::att_gt(clustervars = ) and
#    did::aggte(), bstrap = FALSE: ATT(g,t) and all aggregations, 1e-6.
# 2. The bootstrap draws are random, so the bootstrap results cannot match
#    did exactly. Instead:
#    (a) with many draws (20000) the bootstrap standard errors (IQR-based
#        scale, as in did) are close to the analytic ones (relative 2%) on a
#        large simulated panel (4000 individuals), for Rademacher and
#        Mammen multipliers, unclustered and clustered. On mpdta the
#        interquartile scale of the Rademacher bootstrap distribution is up
#        to 8% above the analytic standard error (a finite-sample feature
#        of the IQR scale with skewed influence functions, shared by did),
#        so on mpdta the bootstrap standard errors are compared with did's
#        own bootstrap standard errors with the same number of draws
#        (relative 3%; 6% when clustering by the 29 states, where the
#        interquartile scale of either package moves by 2 to 4% across
#        seeds);
#    (b) the simultaneous critical value of the ATT(g,t) and of the
#        dynamic, group and calendar aggregations is close to did's with
#        the same number of draws (20000) on mpdta (relative 5%), and at
#        least the normal 1.96;
#    (c) a small seeded Monte Carlo (300 samples, 499 draws each) of the
#        simultaneous coverage of the uniform band of the event study:
#        within 0.04 of 0.95;
#    (d) a fit with a seed is reproduced exactly and leaves the global
#        random number generator state unchanged.
# Data: did::mpdta (500 US counties, 2003-2007; state = county FIPS
# code %/% 1000).
source("validation/common.R")
suppressPackageStartupMessages(library(panglm))
data(mpdta, package = "did")
mpdta$state <- mpdta$countyreal %/% 1000
idx <- c("countyreal", "year")
out <- list()
add <- function(...) out[[length(out) + 1L]] <<- vcompare(...)
types <- c("simple", "dynamic", "group", "calendar")

## 1. clustered analytic standard errors ------------------------------------
for (m in c("dr", "reg")) for (cg in c("never", "notyet")) {
  ref <- did::att_gt("lemp", tname = "year", idname = "countyreal", gname = "first.treat",
                     xformla = ~lpop, data = mpdta, est_method = m,
                     control_group = if (cg == "never") "nevertreated" else "notyettreated",
                     clustervars = c("countyreal", "state"), bstrap = FALSE, cband = FALSE)
  fit <- attgt(lemp ~ lpop, mpdta, idx, "first.treat", control = cg, method = m,
               cluster = "state")
  mod <- sprintf("attgt %s %s + lpop, clustered by state", cg, m)
  add(mod, "mpdta", "ATT(g,t)", fit$estimates$estimate, list(did = ref$att), 1e-8)
  add(mod, "mpdta", "SE ATT(g,t)", fit$estimates$se, list(did = ref$se), 1e-6)
  add(mod, "mpdta", "pre-trend Wald", fit$pretrend$statistic, list(did = ref$W), 1e-6)
  for (type in types) {
    ag <- did::aggte(ref, type = type, bstrap = FALSE, cband = FALSE)
    s <- summary(fit, type = type)
    add(mod, "mpdta", paste(type, "overall SE"), s$overall$se, list(did = ag$overall.se), 1e-6)
    if (type != "simple") {
      add(mod, "mpdta", paste(type, "SE"), s$estimates$se, list(did = ag$se.egt), 1e-6)
    }
  }
}

## 2a. bootstrap standard errors --------------------------------------------
B <- 20000
sim_panel <- function(seed, N, years = 2011:2016) {
  set.seed(seed)
  d <- expand.grid(id = seq_len(N), year = years)
  coh <- sample(c(2013, 2015, NA), N, replace = TRUE)
  d$first <- coh[d$id]
  d$region <- ceiling(d$id / 20)
  a <- rnorm(N)
  u <- rnorm(max(d$region), sd = 0.3)
  e <- d$year - d$first
  d$y <- a[d$id] + 0.2 * (d$year - min(years)) + u[d$region] * (d$year - min(years)) / 3 +
    ifelse(!is.na(e) & e >= 0, 1 + 0.5 * e, 0) + rnorm(nrow(d))
  d
}
big <- sim_panel(1, 4000)
for (mult in c("rademacher", "mammen")) for (clus in c(FALSE, TRUE)) {
  cl <- if (clus) "region" else NULL
  fa <- attgt(y ~ 1, big, c("id", "year"), "first", cluster = cl)
  fb <- attgt(y ~ 1, big, c("id", "year"), "first", cluster = cl, bstrap = TRUE, biters = B,
              multiplier = mult, seed = 11)
  mod <- sprintf("attgt bootstrap %s%s", mult, if (clus) ", clustered by region" else "")
  add(mod, "simulated, 4000 individuals", "bootstrap SE / analytic SE ATT(g,t)",
      fb$estimates$se / fa$estimates$se, list(analytic = rep(1, nrow(fa$estimates))), 0.02)
  for (type in c("dynamic", "group")) {
    sa <- summary(fa, type = type)
    sb <- summary(fb, type = type)
    add(mod, "simulated, 4000 individuals", paste("bootstrap SE / analytic SE", type),
        c(sb$overall$se, sb$estimates$se) / c(sa$overall$se, sa$estimates$se),
        list(analytic = rep(1, nrow(sa$estimates) + 1L)), 0.02)
  }
}

## 2a, 2b. against did's bootstrap on mpdta ---------------------------------
for (clus in c(FALSE, TRUE)) {
  cv <- if (clus) c("countyreal", "state") else NULL
  set.seed(101)
  ref <- did::att_gt("lemp", tname = "year", idname = "countyreal", gname = "first.treat",
                     xformla = ~lpop, data = mpdta, bstrap = TRUE, cband = TRUE, biters = B,
                     clustervars = cv)
  fit <- attgt(lemp ~ lpop, mpdta, idx, "first.treat", cband = TRUE, biters = B, seed = 7,
               cluster = if (clus) "state" else NULL)
  mod <- sprintf("attgt never dr + lpop, bootstrap%s", if (clus) ", clustered by state" else "")
  tol_se <- if (clus) 0.06 else 0.03
  add(mod, "mpdta", "bootstrap SE ATT(g,t) / did", fit$estimates$se / ref$se,
      list(did = rep(1, length(ref$se))), tol_se)
  add(mod, "mpdta", "uniform critical value ATT(g,t) / did", fit$boot$crit / ref$c,
      list(did = 1), 0.05)
  add(mod, "mpdta", "uniform critical value ATT(g,t) >= 1.96",
      as.numeric(fit$boot$crit >= stats::qnorm(0.975)), list(rule = 1), 0)
  for (type in c("dynamic", "group", "calendar")) {
    set.seed(102)
    ag <- did::aggte(ref, type = type, biters = B)
    s <- summary(fit, type = type)
    add(mod, "mpdta", paste("bootstrap SE", type, "/ did"),
        c(s$overall$se, s$estimates$se) / c(ag$overall.se, ag$se.egt),
        list(did = rep(1, length(ag$se.egt) + 1L)), tol_se)
    add(mod, "mpdta", paste("uniform critical value", type, "/ did"), s$crit / ag$crit.val.egt,
        list(did = 1), 0.05)
    add(mod, "mpdta", paste("uniform critical value", type, ">= 1.96"),
        as.numeric(s$crit >= stats::qnorm(0.975)), list(rule = 1), 0)
  }
}

## 2c. coverage of the uniform band -----------------------------------------
R <- 300
cover <- vapply(seq_len(R), function(r) {
  d <- sim_panel(1000 + r, 300)
  fit <- attgt(y ~ 1, d, c("id", "year"), "first", cband = TRUE, biters = 499, seed = r)
  s <- summary(fit, type = "dynamic")$estimates
  truth <- ifelse(s$level >= 0, 1 + 0.5 * s$level, 0)
  z <- stats::qnorm(0.975)
  c(band = all(s$band_lower <= truth & truth <= s$band_upper),
    pointwise = all(abs(s$estimate - truth) <= z * s$se))
}, numeric(2))
cat(sprintf("Simultaneous coverage over %d samples: uniform band %.3f, pointwise intervals %.3f\n",
            R, mean(cover["band", ]), mean(cover["pointwise", ])))
add("attgt dynamic, uniform band", "simulated, 300 individuals, 300 samples",
    "simultaneous coverage", mean(cover["band", ]), list(nominal = 0.95), 0.04)

## 2d. reproducibility ------------------------------------------------------
set.seed(2026)
before <- .Random.seed
f1 <- attgt(lemp ~ lpop, mpdta, idx, "first.treat", cband = TRUE, biters = 999, seed = 3)
after <- .Random.seed
f2 <- attgt(lemp ~ lpop, mpdta, idx, "first.treat", cband = TRUE, biters = 999, seed = 3)
s1 <- summary(f1, type = "dynamic")
s2 <- summary(f2, type = "dynamic")
add("attgt bootstrap", "mpdta", "same seed, same fit and summary",
    as.numeric(identical(f1$estimates, f2$estimates) && identical(s1$estimates, s2$estimates)),
    list(rule = 1), 0)
add("attgt bootstrap", "mpdta", "global RNG state restored", as.numeric(identical(before, after)),
    list(rule = 1), 0)

do.call(vsave, c(out, name = "did-boot"))
