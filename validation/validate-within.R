# Validation of the fixed-effects ("within") panel GLMs.
#
# Same model, same data, different packages:
#   gaussian one-way / two-way  plm within, fixest::feols, lm with dummies
#   poisson one-way / two-way   fixest::fepois, glm with dummies, pglm
#                               (conditional log-likelihood, one-way)
#   negbin one-way / two-way    MASS::glm.nb with dummies, fixest::fenegbin
#   gamma one-way               glm(Gamma("log")) with dummies, fixest::feglm
#   logit (conditional)         survival::clogit (exact), and for T = 2 the
#                               logit on first differences (glm + sandwich)
#
# Estimator notes:
# - panglm's fixed-effects negative binomial is the unconditional
#   dummy-variable NB2 MLE (Allison and Waterman 2002), not the conditional
#   model of Hausman, Hall and Griliches (pglm model = "within"), so it is
#   compared with glm.nb and fenegbin, which fit the same dummy model.
#   Classical covariances hold theta fixed (as glm.nb does); fenegbin's
#   covariance treats theta as a free parameter in the information, so its
#   standard errors are a different quantity and are not compared.
# - The fixed-effects logit is the conditional (Chamberlain) logit.
#   fixest::feglm(binomial) | id is the unconditional dummy-variable logit,
#   a different estimator (incidental-parameter bias), so it is not used.
#
# Small-sample conventions (k = number of slopes, G = number of clusters):
#   classical  residual variance with df n - k - (fixed effects)  (plm, lm)
#   HC1        n / (n - k)                 (plm white1 HC1; fixest hetero
#                                            with ssc(K.fixef = "none"))
#   cluster    G / (G - 1) (n - 1) / (n - k) (plm arellano "sss"; fixest
#                                            with ssc(K.fixef = "none"))
#   The fixest default with clustering counts one parameter for the nested
#   individual effects (K = k + 1, as in Stata's xtreg, fe), hence the
#   explicit ssc(). The one-way negative binomial counts the k + G
#   dummy-model parameters, n / (n - k - G) and G / (G - 1) (n - 1) /
#   (n - k - G), which is what sandwich::vcovHC / vcovCL (type "HC1") give
#   on the glm.nb dummy fit. The two-way negative binomial uses
#   n / (n - k): the reference is the glm.nb dummy-fit HC0 sandwich times
#   that same factor.
#
# The simulated count panels have no individual with an all-zero outcome:
# fixest drops such individuals from n and G, while panglm's Poisson
# keeps them in n and G (their score is zero), a convention difference
# that changes the HC1 / cluster factors slightly.
#
# Run from the package root: Rscript validation/validate-within.R
source("validation/common.R")
suppressPackageStartupMessages({
  library(plm)
  library(fixest)
  library(sandwich)
  library(MASS)
  library(survival)
})
setFixest_nthreads(1)
se <- function(v) sqrt(diag(v))
nofe <- ssc(K.fixef = "none")
tol <- 1e-6
out <- list()

# ---- data -----------------------------------------------------------------
data("Grunfeld", package = "plm")
grun <- Grunfeld[-c(3, 17, 45, 46, 120, 151:155), ] # unbalanced
gname <- "Grunfeld unbalanced (190 obs)"

# simulated health-like panel: 150 individuals, 4 to 6 visits, individual
# effects correlated with x1, a binary exposure x2 and a visit trend
sim_panel <- function(G = 150, T = 6, seed = 1) {
  set.seed(seed)
  d <- expand.grid(time = seq_len(T), id = seq_len(G))[, c("id", "time")]
  a <- rnorm(G, 0, 0.7)
  d$x1 <- rnorm(nrow(d)) + 0.5 * a[d$id]
  d$x2 <- rbinom(nrow(d), 1, 0.4)
  d <- d[!(d$time > 4 & d$id %% 5 == 0), ]
  eta <- a[d$id] + 0.4 * d$x1 - 0.3 * d$x2 + 0.1 * d$time
  d$count <- rpois(nrow(d), rgamma(nrow(d), 2, 2) * exp(1 + eta))
  d$cost <- rgamma(nrow(d), shape = 2, rate = 2 / exp(5 + eta))
  d$y <- rbinom(nrow(d), 1, plogis(eta))
  rownames(d) <- NULL
  d
}
sim <- sim_panel()
stopifnot(all(tapply(sim$count, sim$id, sum) > 0))
sname <- "simulated (150 ids, 840 obs)"
idx <- c("id", "time")

# ---- gaussian one-way and two-way -----------------------------------------
pg <- pdata.frame(grun, c("firm", "year"))
for (eff in c("individual", "twoways")) {
  label <- paste("FE gaussian", eff)
  p <- panglm(inv ~ value + capital, grun, c("firm", "year"), model = "within", effect = eff)
  pl <- plm(inv ~ value + capital, pg, model = "within", effect = eff)
  if (eff == "individual") {
    fx <- feols(inv ~ value + capital | firm, grun)
    lsdv <- lm(inv ~ value + capital + factor(firm), grun)
  } else {
    fx <- feols(inv ~ value + capital | firm + year, grun)
    lsdv <- lm(inv ~ value + capital + factor(firm) + factor(year), grun)
  }
  b <- c("value", "capital")
  out[[label]] <- rbind(
    vcompare(label, gname, "coef", coef(p),
             list(plm = coef(pl), fixest = coef(fx), lm = coef(lsdv)[b]), tol),
    vcompare(label, gname, "se_classical", se(vcov(p)),
             list(plm = se(vcov(pl)), fixest = se(vcov(fx, vcov = "iid")),
                  lm = se(vcov(lsdv))[b]), tol),
    vcompare(label, gname, "se_HC1", se(vcov(p, type = "HC1")),
             list(plm = se(plm::vcovHC(pl, method = "white1", type = "HC1")),
                  fixest = se(vcov(fx, vcov = "hetero", ssc = nofe))), tol),
    vcompare(label, gname, "se_cluster", se(vcov(p, type = "cluster")),
             list(plm = se(plm::vcovHC(pl, method = "arellano", type = "sss")),
                  fixest = se(vcov(fx, vcov = ~firm, ssc = nofe))), tol),
    vcompare(label, gname, "logLik", as.numeric(logLik(p)),
             list(fixest = as.numeric(logLik(fx)), lm = as.numeric(logLik(lsdv))), tol)
  )
}

# ---- poisson one-way and two-way --------------------------------------------
for (eff in c("individual", "twoways")) {
  label <- paste("FE poisson", eff)
  p <- panglm(count ~ x1 + x2, sim, idx, model = "within", family = "poisson", effect = eff)
  if (eff == "individual") {
    fx <- fepois(count ~ x1 + x2 | id, sim)
    g <- glm(count ~ x1 + x2 + factor(id), poisson, sim, control = glm.control(epsilon = 1e-12))
  } else {
    fx <- fepois(count ~ x1 + x2 | id + time, sim)
    g <- glm(count ~ x1 + x2 + factor(id) + factor(time), poisson, sim,
             control = glm.control(epsilon = 1e-12))
  }
  b <- c("x1", "x2")
  out[[label]] <- rbind(
    vcompare(label, sname, "coef", coef(p), list(fixest = coef(fx), glm = coef(g)[b]), tol),
    vcompare(label, sname, "se_classical", se(vcov(p)),
             list(fixest = se(vcov(fx, vcov = "iid", ssc = ssc(K.adj = FALSE))),
                  glm = se(vcov(g))[b]), tol),
    vcompare(label, sname, "se_HC1", se(vcov(p, type = "HC1")),
             list(fixest = se(vcov(fx, vcov = "hetero", ssc = nofe))), tol),
    vcompare(label, sname, "se_cluster", se(vcov(p, type = "cluster")),
             list(fixest = se(vcov(fx, vcov = ~id, ssc = nofe))), tol),
    if (eff == "individual") {
      # one-way: panglm reports the conditional (multinomial) log-likelihood,
      # as pglm does; the dummy-model log-likelihood of glm / fixest exceeds
      # it by the Poisson log-likelihood of the individual totals at their
      # MLE, sum_i log dpois(Y_i, Y_i)
      Yi <- tapply(sim$count, sim$id, sum)
      vcompare(label, sname, "logLik_conditional", as.numeric(logLik(p)),
               list(pglm = as.numeric(logLik(pglm::pglm(count ~ x1 + x2, sim, index = idx,
                                                        model = "within", family = poisson))),
                    glm = as.numeric(logLik(g)) - sum(dpois(Yi, Yi, log = TRUE))), tol)
    } else {
      # two-way: panglm reports the dummy-model (unconditional) log-likelihood
      vcompare(label, sname, "logLik", as.numeric(logLik(p)),
               list(fixest = as.numeric(logLik(fx)), glm = as.numeric(logLik(g))), tol)
    }
  )
}
# a cluster that cuts across individuals (by visit) uses observation scores
p <- panglm(count ~ x1 + x2, sim, idx, model = "within", family = "poisson")
fx <- fepois(count ~ x1 + x2 | id, sim)
out$pois_time_cluster <- vcompare(
  "FE poisson individual", sname, "se_cluster_by_time",
  se(vcov(p, type = "cluster", cluster = sim$time)),
  list(fixest = se(vcov(fx, vcov = ~time, ssc = nofe))), tol
)

# ---- negative binomial (dummy-variable NB2) ---------------------------------
# glm.nb(), fenegbin() and panglm agree to about 1e-8 on the one-way model.
# The two-way fit profiles theta with stats::optimize() and stops on a
# relative log-likelihood change below tol, so its coefficients and theta are
# accurate to about 1e-6 only: compared at 1e-5.
n <- nrow(sim)
for (eff in c("individual", "twoways")) {
  label <- paste("FE negbin", eff)
  tl <- if (eff == "individual") tol else 1e-5
  p <- panglm(count ~ x1 + x2, sim, idx, model = "within", family = "negbin", effect = eff)
  if (eff == "individual") {
    fx <- fenegbin(count ~ x1 + x2 | id, sim)
    nb <- glm.nb(count ~ x1 + x2 + factor(id), sim,
                 control = glm.control(epsilon = 1e-13, maxit = 200))
    ref_hc1 <- vcovHC(nb, type = "HC1")
    ref_cl <- vcovCL(nb, cluster = ~id, type = "HC1")
  } else {
    fx <- fenegbin(count ~ x1 + x2 | id + time, sim)
    nb <- glm.nb(count ~ x1 + x2 + factor(id) + factor(time), sim,
                 control = glm.control(epsilon = 1e-13, maxit = 200))
    ref_hc1 <- vcovHC(nb, type = "HC0") * n / (n - 2)
    ref_cl <- vcovCL(nb, cluster = ~id, type = "HC0", cadjust = TRUE) * (n - 1) / (n - 2)
  }
  b <- c("x1", "x2")
  out[[label]] <- rbind(
    vcompare(label, sname, "coef", coef(p), list(MASS = coef(nb)[b], fixest = coef(fx)[b]), tl),
    vcompare(label, sname, "theta", p$theta, list(MASS = nb$theta, fixest = fx$theta), tl),
    vcompare(label, sname, "se_classical", se(vcov(p)), list(MASS = se(vcov(nb))[b]), tl),
    vcompare(label, sname, "se_HC1", se(vcov(p, type = "HC1")), list(sandwich = se(ref_hc1)[b]), tl),
    vcompare(label, sname, "se_cluster", se(vcov(p, type = "cluster")),
             list(sandwich = se(ref_cl)[b]), tl),
    vcompare(label, sname, "logLik", as.numeric(logLik(p)),
             list(MASS = as.numeric(logLik(nb)), fixest = as.numeric(logLik(fx))), tl)
  )
}

# ---- gamma (log link) one-way -----------------------------------------------
# The IRLS stops on a relative deviance change below tol, so with the default
# tol = 1e-10 the coefficients are accurate to about 1e-6 only; the estimator
# is compared at tol = 1e-14 (1e-6) and the default fit at 1e-5.
p <- panglm(cost ~ x1 + x2, sim, idx, model = "within", family = "gamma", tol = 1e-14)
p_default <- panglm(cost ~ x1 + x2, sim, idx, model = "within", family = "gamma")
g <- glm(cost ~ x1 + x2 + factor(id), Gamma("log"), sim, control = glm.control(epsilon = 1e-13, maxit = 100))
fx <- feglm(cost ~ x1 + x2 | id, sim, family = Gamma("log"), glm.tol = 1e-12)
b <- c("x1", "x2")
out$gamma <- rbind(
  vcompare("FE gamma individual", sname, "coef", coef(p),
           list(glm = coef(g)[b], fixest = coef(fx)), tol),
  vcompare("FE gamma individual", sname, "coef_default_tol", coef(p_default),
           list(glm = coef(g)[b]), 1e-5),
  vcompare("FE gamma individual", sname, "se_classical", se(vcov(p)),
           list(glm = se(vcov(g))[b]), tol),
  vcompare("FE gamma individual", sname, "se_HC1", se(vcov(p, type = "HC1")),
           list(fixest = se(vcov(fx, vcov = "hetero", ssc = nofe))), tol),
  vcompare("FE gamma individual", sname, "se_cluster", se(vcov(p, type = "cluster")),
           list(fixest = se(vcov(fx, vcov = ~id, ssc = nofe))), tol),
  vcompare("FE gamma individual", sname, "logLik", as.numeric(logLik(p)),
           list(glm = as.numeric(logLik(g))), tol)
)

# ---- conditional logit ------------------------------------------------------
# panglm maximizes the conditional likelihood with BFGS (reltol = tol on the
# objective), so with the default tol = 1e-10 the coefficients are accurate
# to about 1e-5 only; tol = 1e-14 is used here to compare the estimator
# itself at 1e-6, and the default fit is compared at 1e-4.
pb <- suppressMessages(panglm(y ~ x1 + x2, sim, idx, model = "within", family = "binomial",
                              tol = 1e-14))
pb_default <- suppressMessages(panglm(y ~ x1 + x2, sim, idx, model = "within",
                                      family = "binomial"))
cl <- suppressWarnings(clogit(y ~ x1 + x2 + strata(id), sim, method = "exact",
                              control = coxph.control(eps = 1e-12, iter.max = 100)))
out$clogit <- rbind(
  vcompare("FE logit (conditional)", sname, "coef", coef(pb), list(survival = coef(cl)), tol),
  vcompare("FE logit (conditional)", sname, "coef_default_tol", coef(pb_default),
           list(survival = coef(cl)), 1e-4),
  vcompare("FE logit (conditional)", sname, "se_classical", se(vcov(pb)),
           list(survival = se(vcov(cl))), tol),
  vcompare("FE logit (conditional)", sname, "logLik", as.numeric(logLik(pb)),
           list(survival = cl$loglik[2]), tol)
)

# Robust covariances: with T = 2 the conditional logit is exactly the logit
# without intercept of y_i2 on x_i2 - x_i1 among individuals whose outcome
# changes, so glm + sandwich give independent stratum-level scores. panglm's
# HC1 and cluster factors are then both m / (m - k) for m switchers, which is
# sandwich's HC1 (and vcovCL type HC1 with one cluster per switcher).
sim2 <- sim[sim$time <= 2, ]
p2 <- suppressMessages(panglm(y ~ x1 + x2, sim2, idx, model = "within", family = "binomial",
                              tol = 1e-14))
w1 <- sim2[sim2$time == 1, ]; w2 <- sim2[sim2$time == 2, ]
stopifnot(identical(w1$id, w2$id))
fd <- data.frame(y = w2$y, dx1 = w2$x1 - w1$x1, dx2 = w2$x2 - w1$x2)[w1$y != w2$y, ]
gd <- glm(y ~ 0 + dx1 + dx2, binomial, fd, control = glm.control(epsilon = 1e-14, maxit = 100))
cl2 <- suppressWarnings(clogit(y ~ x1 + x2 + strata(id), sim2, method = "exact"))
lab2 <- "FE logit (conditional), T = 2"
s2 <- "simulated, waves 1-2"
out$clogit2 <- rbind(
  vcompare(lab2, s2, "coef", coef(p2), list(glm_fd = coef(gd), survival = coef(cl2)), tol),
  vcompare(lab2, s2, "se_classical", se(vcov(p2)),
           list(glm_fd = se(vcov(gd)), survival = se(vcov(cl2))), tol),
  vcompare(lab2, s2, "se_HC1", se(vcov(p2, type = "HC1")),
           list(sandwich = se(vcovHC(gd, type = "HC1"))), tol),
  vcompare(lab2, s2, "se_cluster", se(vcov(p2, type = "cluster")),
           list(sandwich = se(vcovCL(gd, cluster = seq_len(nrow(fd)), type = "HC1"))), tol)
)

do.call(vsave, c(unname(out), name = "within"))
