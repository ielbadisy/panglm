# Generic methods: predict() and tidy() against established packages.
#
# predict():
# * paniv() pooled, fixed-effects (one-way and two-way) and random-effects
#   2SLS against predict() of plm::plm() IV fits (in sample, fixed effects
#   included for the within models) and AER::ivreg() (pooled, on newdata);
#   panglm's newdata predictions on the fitting data must reproduce its
#   in-sample predictions. Data: plm::Crime, and an unbalanced version.
# * pangmm() difference and system GMM against plm::pgmm(): pgmm stores, per
#   firm, the transformed equations (first differences, then levels for
#   system GMM) as [y, X]; X b is compared with predict(type = "difference")
#   and predict(type = "level") at the same firm-years. Rows that plm pads
#   for missing years are left out (panglm predicts NA there). Data:
#   plm::EmplUK (unbalanced).
# * panselect() on one period (pooled, no Mundlak means, one lambda) is
#   Heckman's two-step estimator: unconditional, conditional and selection
#   predictions against predict() of sampleSelection::selection(method =
#   "2step"). Data: sampleSelection::Mroz87. The first-step probits are
#   maximised by different routines (glm.fit, maxLik), so the comparison
#   uses a tolerance of 1e-6.
#
# tidy():
# * estimates and standard errors of tidy() against the reference fits'
#   coefficients and covariances: paniv vs plm (classical), pangmm vs pgmm
#   (vcovHC, Windmeijer), pooled panqr vs quantreg::rq() estimates, attgt vs
#   did::att_gt() (analytic standard errors), etwfe vs etwfe::emfx() (simple
#   ATT). tidy() intervals must equal confint() of the same fit.
source("validation/common.R")
suppressPackageStartupMessages({ library(plm); library(sampleSelection) })
out <- list()
add <- function(...) out[[length(out) + 1L]] <<- vcompare(...)
se <- function(v) sqrt(diag(v))

## paniv ---------------------------------------------------------------------
data("Crime", package = "plm")
set.seed(20261025)
CrimeU <- Crime[-sample(nrow(Crime), 60), ]
idx <- c("county", "year")
exog <- "lprbconv + lprbpris + lavgsen + ldensity + lwcon + lpctymle"
f_iv <- stats::as.formula(paste("lcrmrte ~ lprbarr + lpolpc +", exog, "| ltaxpc + lmix +", exog))
for (dname in c("Crime", "Crime unbalanced")) {
  d <- if (dname == "Crime") Crime else CrimeU
  pd <- plm::pdata.frame(d, index = idx)
  key <- match(paste(d$county, d$year), paste(pd$county, pd$year))
  specs <- list(pooling = list("pooling", "individual"), within = list("within", "individual"),
                `within twoways` = list("within", "twoways"), random = list("random", "individual"))
  for (nm in names(specs)) {
    a <- paniv(f_iv, d, index = idx, model = specs[[nm]][[1]], effect = specs[[nm]][[2]])
    b <- plm::plm(f_iv, pd, model = specs[[nm]][[1]], effect = specs[[nm]][[2]])
    m <- paste("paniv", nm)
    add(m, dname, "predict (in sample)", predict(a), list(plm = as.numeric(predict(b))[key]), 1e-8)
    add(m, dname, "predict (newdata = data)", predict(a, newdata = d),
        list(`panglm in sample` = predict(a)), 1e-10)
    td <- generics::tidy(a, conf.int = TRUE)
    add(m, dname, "tidy estimate", td$estimate, list(plm = coef(b)[td$term]), 1e-8)
    add(m, dname, "tidy std.error", td$std.error, list(plm = se(vcov(b))[td$term]), 1e-8)
    add(m, dname, "tidy conf.low", td$conf.low, list(confint = confint(a)[td$term, 1]), 1e-12)
  }
  a <- paniv(f_iv, d, index = idx, model = "pooling")
  iv <- AER::ivreg(f_iv, data = d)
  add("paniv pooling", dname, "predict (newdata)", predict(a, newdata = d),
      list(AER = predict(iv, newdata = d)), 1e-8)
}

## pangmm --------------------------------------------------------------------
data("EmplUK", package = "plm")
gmm_specs <- list(
  `difference GMM` = list(
    f = log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) + log(capital) +
      lag(log(output), 0:1) | lag(log(emp), 2:99),
    transformation = "d"),
  `system GMM` = list(
    f = log(emp) ~ lag(log(emp), 1) + lag(log(wage), 0:1) + lag(log(capital), 0:1) |
      lag(log(emp), 2:99) + lag(log(wage), 2:99) + lag(log(capital), 2:99),
    transformation = "ld"))
for (nm in names(gmm_specs)) {
  s <- gmm_specs[[nm]]
  a <- pangmm(s$f, EmplUK, index = c("firm", "year"), model = "twostep",
              transformation = s$transformation)
  b <- plm::pgmm(s$f, EmplUK, index = c("firm", "year"), effect = "twoways",
                 model = "twosteps", transformation = s$transformation)
  dif <- predict(a, newdata = EmplUK, type = "difference")
  lev <- predict(a, newdata = EmplUK, type = "level")
  mine <- ref <- numeric(0)
  for (firm in names(b$model)) {
    M <- b$model[[firm]]
    yr <- as.numeric(rownames(M))
    block_start <- c(1L, which(diff(yr) < 0) + 1L)
    block <- findInterval(seq_along(yr), block_start)
    rows <- which(EmplUK$firm == as.numeric(firm))
    p <- ifelse(block == 1L, dif[rows][match(yr, EmplUK$year[rows])],
                lev[rows][match(yr, EmplUK$year[rows])])
    keep <- !is.na(p)
    mine <- c(mine, p[keep])
    ref <- c(ref, as.numeric(M[keep, -1, drop = FALSE] %*% coef(b)))
  }
  add(paste("pangmm", nm), "EmplUK", "predict (transformed equations)", mine, list(plm = ref), 1e-8)
  td <- generics::tidy(a)
  # plm names lag terms differently; the coefficients are in the same order
  add(paste("pangmm", nm), "EmplUK", "tidy estimate", td$estimate,
      list(plm = unname(coef(b))), 1e-6)
  add(paste("pangmm", nm), "EmplUK", "tidy std.error", td$std.error,
      list(plm = unname(se(plm::vcovHC(b)))), 1e-6)
}

## panselect -----------------------------------------------------------------
data("Mroz87", package = "sampleSelection")
Mroz87$id <- seq_len(nrow(Mroz87))
Mroz87$t <- 1L
Mroz87$wage[Mroz87$lfp == 0] <- NA
fs <- lfp ~ age + I(age^2) + faminc + kids5 + educ
fo <- wage ~ exper + I(exper^2) + educ + city
a <- panselect(fs, fo, data = Mroz87, index = c("id", "t"), model = "pooling",
               mundlak = FALSE, lambda = "common", boot = 0)
h <- sampleSelection::selection(fs, fo, data = Mroz87, method = "2step")
ho <- coef(h, part = "outcome")
add("panselect (one period)", "Mroz87", "coefficients", coef(a),
    list(sampleSelection = c(ho[names(ho) != "invMillsRatio"], ho[["invMillsRatio"]])), 1e-6)
add("panselect (one period)", "Mroz87", "predict unconditional", predict(a, newdata = Mroz87),
    list(sampleSelection = predict(h, newdata = Mroz87, part = "outcome", type = "unconditional")),
    1e-6)
add("panselect (one period)", "Mroz87", "predict conditional",
    predict(a, newdata = Mroz87, type = "conditional"),
    list(sampleSelection = predict(h, newdata = Mroz87, part = "outcome",
                                   type = "conditional")[, "E[yo|ys=1]"]), 1e-6)
add("panselect (one period)", "Mroz87", "predict selection",
    predict(a, newdata = Mroz87, type = "selection"),
    list(sampleSelection = predict(h, newdata = Mroz87, part = "selection", type = "response")),
    1e-6)
add("panselect (one period)", "Mroz87", "tidy estimate", generics::tidy(a)$estimate,
    list(coef = coef(a)), 1e-12)

## panqr (pooled) ------------------------------------------------------------
data(copd)
taus <- c(0.25, 0.5, 0.75)
a <- panqr(fev1 ~ crp + age + treatment, data = copd, index = c("id", "visit"), tau = taus,
           model = "pooling")
td <- generics::tidy(a, conf.int = TRUE)
r <- quantreg::rq(fev1 ~ crp + age + treatment, tau = taus, data = copd)
add("panqr pooling", "copd", "tidy estimate", td$estimate,
    list(quantreg = as.numeric(stats::coef(r))), 1e-8)
ci <- do.call(rbind, lapply(taus, function(tau) confint(a, tau = tau)))
add("panqr pooling", "copd", "tidy conf.high", td$conf.high, list(confint = ci[, 2]), 1e-12)

## attgt and etwfe -----------------------------------------------------------
data(mpdta, package = "did")
ref <- did::att_gt("lemp", tname = "year", idname = "countyreal", gname = "first.treat",
                   xformla = ~lpop, data = mpdta, bstrap = FALSE, cband = FALSE)
fit <- attgt(lemp ~ lpop, mpdta, c("countyreal", "year"), "first.treat")
td <- generics::tidy(fit)
add("attgt", "mpdta", "tidy estimate", td$estimate, list(did = ref$att), 1e-8)
add("attgt", "mpdta", "tidy std.error", td$std.error, list(did = ref$se), 1e-6)
ag <- did::aggte(ref, type = "dynamic", bstrap = FALSE, cband = FALSE)
tdd <- generics::tidy(fit, type = "dynamic")
add("attgt", "mpdta", "tidy dynamic estimate", tdd$estimate, list(did = ag$att.egt), 1e-8)
add("attgt", "mpdta", "tidy dynamic std.error", tdd$std.error, list(did = ag$se.egt), 1e-6)

et <- etwfe(lemp ~ lpop, mpdta, c("countyreal", "year"), "first.treat")
e <- etwfe::emfx(suppressMessages(etwfe::etwfe(lemp ~ lpop, tvar = year, gvar = first.treat,
                                               data = mpdta, vcov = ~countyreal)),
                 type = "simple")
tde <- generics::tidy(et)
add("etwfe", "mpdta", "tidy simple ATT", tde$estimate, list(etwfe = e$estimate), 1e-8)
add("etwfe", "mpdta", "tidy simple SE", tde$std.error, list(etwfe = e$std.error), 1e-6)

vsave(do.call(rbind, out), name = "methods")
