# Random-intercept models.
#
# Gaussian: panglm(model = "random") is the Swamy-Arora FGLS estimator and
# is compared with plm::plm(model = "random", random.method = "swar") on a
# balanced panel (Grunfeld), an unbalanced one (Grunfeld with 50 of 200
# rows dropped at random) and the unbalanced RAND Health Insurance
# Experiment panel (pglm::HealthIns): coefficients, standard errors, the
# variance components and the quasi-demeaning weight theta_i of every
# individual. The Hausman test on the unbalanced panel is compared with
# plm::phtest().
#
# Binomial (logit, probit) and ordered (logit, probit): random-intercept
# likelihoods integrated by adaptive Gauss-Hermite quadrature. At the
# default R = 21, binomial fits are compared with lme4::glmer(nAGQ = 21)
# on pglm::UnionWage and on the RAND HIE indicator of any medical spending
# (med > 0), and ordered fits with ordinal::clmm(nAGQ = 21) on
# lme4::schizophrenia. glmer's own vcov() is an approximation (it differs
# from the curvature of its own nAGQ = 21 deviance by 1 to 2 percent), so
# standard errors are compared with the inverse of a numerical Hessian
# (numDeriv) of the deviance function glmer optimizes. With R = 1 the fit
# is the Laplace approximation and is compared with glmmTMB (binomial) and
# clmm(nAGQ = 1) (ordered). UnionWage probit, whose random-intercept spread
# is large, is compared at 21 nodes with GLMMadaptive (same node scaling as
# panglm) and at 61 nodes with GLMMadaptive and glmer(nAGQ = 25, its
# maximum). glmer(nAGQ = 1) is not used as the Laplace
# reference: its reported log-likelihood is not the exact Laplace value
# (the exact Laplace value at glmer's own estimates is higher), and for the
# probit link it scales by the expected rather than the observed curvature.
source("validation/common.R")
suppressPackageStartupMessages(library(plm))
data("Grunfeld", package = "plm")
data("HealthIns", package = "pglm")

gaussian_re <- function(formula, data, index, label) {
  f <- panglm(formula, data = data, index = index, model = "random")
  p <- plm(formula, data = data, index = index, model = "random", random.method = "swar")
  ec <- ercomp(p)
  ids <- index(p)[[1L]]
  th <- if (length(ec$theta) == 1L) rep(ec$theta, length(ids)) else ec$theta
  th <- as.numeric(tapply(th, ids, `[`, 1L))
  rbind(
    vcompare("gaussian random", label, "coef", coef(f), list(plm = coef(p)), 1e-8),
    vcompare("gaussian random", label, "se", sqrt(diag(vcov(f))),
             list(plm = sqrt(diag(vcov(p)))), 1e-8),
    vcompare("gaussian random", label, "sigma2 (idios, id)",
             c(idios = f$sigma_v2, id = f$sigma_mu2),
             list(plm = ec$sigma2[c("idios", "id")]), 1e-8),
    vcompare("gaussian random", label, "theta_i", f$theta, list(plm = th), 1e-8)
  )
}

set.seed(1)
grun_u <- Grunfeld[sort(sample(nrow(Grunfeld), 150)), ]
hie <- transform(HealthIns, lcoins = log(1 + coins))
hie_fm <- med ~ lcoins + disease + sex + age + size + child

hausman_u <- local({
  fe <- panglm(inv ~ value + capital, data = grun_u, index = c("firm", "year"), model = "within")
  re <- panglm(inv ~ value + capital, data = grun_u, index = c("firm", "year"), model = "random")
  h <- panhausman(fe, re)
  ph <- phtest(plm(inv ~ value + capital, data = grun_u, index = c("firm", "year"), model = "within"),
               plm(inv ~ value + capital, data = grun_u, index = c("firm", "year"), model = "random"))
  vcompare("gaussian random", "Grunfeld unbalanced", "Hausman statistic",
           c(statistic = h$statistic), list(plm = ph$statistic), 1e-8)
})

suppressPackageStartupMessages({ library(lme4); library(glmmTMB); library(ordinal) })
options(panglm.threads = 6)

binomial_re <- function(formula, data, label, links = c("logit", "probit"),
                        agq21 = links) {
  out <- list()
  for (lk in links) {
    model <- paste("binomial random", lk)
    fr <- update(formula, . ~ . + (1 | id))
    if (lk %in% agq21) {
      f <- panglm(formula, data = data, index = c("id", "year"), model = "random",
                  family = binomialfamily(lk))
      g <- glmer(fr, data = data, family = binomial(lk), nAGQ = 21,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
      dev <- glmer(fr, data = data, family = binomial(lk), nAGQ = 21, devFunOnly = TRUE)
      V <- solve(numDeriv::hessian(dev, c(getME(g, "theta"), fixef(g))) / 2)
      b <- names(fixef(g))
      # probit: glmer scales the nodes by the expected curvature, panglm by
      # the observed one, so the two 21-node rules differ slightly
      tq <- if (lk == "logit") 1e-5 else 5e-5
      out <- c(out, list(
        vcompare(model, label, "coef (R = 21)", coef(f)[b], list(`lme4::glmer nAGQ=21` = fixef(g)), tq),
        vcompare(model, label, "sigma_u (R = 21)", f$sigma,
                 list(`lme4::glmer nAGQ=21` = sqrt(VarCorr(g)$id[1])), tq),
        vcompare(model, label, "logLik (R = 21)", f$loglik,
                 list(`lme4::glmer nAGQ=21` = as.numeric(logLik(g))), 1e-7),
        vcompare(model, label, "se (R = 21)", sqrt(diag(vcov(f)))[b],
                 list(`Hessian of glmer nAGQ=21 deviance` = sqrt(diag(V))[-1]), 1e-5)
      ))
    }
    b <- colnames(model.matrix(formula, data))
    f1 <- panglm(formula, data = data, index = c("id", "year"), model = "random",
                 family = binomialfamily(lk), R = 1)
    t <- glmmTMB(fr, data = data, family = binomial(lk))
    out <- c(out, list(
      vcompare(model, label, "coef (R = 1)", coef(f1)[b], list(glmmTMB = fixef(t)$cond), 1e-5),
      vcompare(model, label, "sigma_u (R = 1)", f1$sigma,
               list(glmmTMB = attr(VarCorr(t)$cond$id, "stddev")), 1e-5),
      vcompare(model, label, "logLik (R = 1)", f1$loglik, list(glmmTMB = as.numeric(logLik(t))), 1e-7),
      vcompare(model, label, "se (R = 1)", sqrt(diag(vcov(f1)))[b],
               list(glmmTMB = sqrt(diag(vcov(t)$cond))), 1e-5)
    ))
  }
  do.call(rbind, out)
}

data("UnionWage", package = "pglm")

# UnionWage probit has a large random-intercept spread (sigma_u = 1.7 on the
# probit scale). 21 nodes then leave a quadrature error of about 1e-3 in the
# log-likelihood, of different size for panglm and glmer because they scale
# the nodes differently (observed versus expected curvature at the mode).
# GLMMadaptive uses the same rule as panglm and agrees at 21 nodes; with
# enough nodes all three converge to the same integrated likelihood.
probit_converged <- local({
  fm <- union ~ wage + exper + rural
  model <- "binomial random probit"
  f21 <- panglm(fm, data = UnionWage, index = c("id", "year"), model = "random",
                family = binomialfamily("probit"), R = 21)
  a21 <- GLMMadaptive::mixed_model(fm, random = ~ 1 | id, data = UnionWage,
                                   family = binomial("probit"), nAGQ = 21)
  f61 <- panglm(fm, data = UnionWage, index = c("id", "year"), model = "random",
                family = binomialfamily("probit"), R = 61)
  a61 <- GLMMadaptive::mixed_model(fm, random = ~ 1 | id, data = UnionWage,
                                   family = binomial("probit"), nAGQ = 61,
                                   control = list(tol3 = 1e-12, iter_EM = 0))
  g25 <- glmer(update(fm, . ~ . + (1 | id)), data = UnionWage, family = binomial("probit"),
               nAGQ = 25, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  b <- names(fixef(g25))
  rbind(
    vcompare(model, "UnionWage", "logLik (R = 21)", f21$loglik,
             list(`GLMMadaptive nAGQ=21` = as.numeric(logLik(a21))), 1e-7),
    vcompare(model, "UnionWage", "coef (R = 61)", coef(f61)[b],
             list(`GLMMadaptive nAGQ=61` = GLMMadaptive::fixef(a61), `lme4::glmer nAGQ=25` = fixef(g25)), 1e-5),
    vcompare(model, "UnionWage", "sigma_u (R = 61)", f61$sigma,
             list(`lme4::glmer nAGQ=25` = sqrt(VarCorr(g25)$id[1])), 1e-5),
    vcompare(model, "UnionWage", "logLik (R = 61)", f61$loglik,
             list(`GLMMadaptive nAGQ=61` = as.numeric(logLik(a61)),
                  `lme4::glmer nAGQ=25` = as.numeric(logLik(g25))), 1e-7)
  )
})
hie$anymed <- as.numeric(hie$med > 0)
hie$age10 <- hie$age / 10

data("schizophrenia", package = "lme4")
ordered_re <- function() {
  out <- list()
  b <- c("TxDrug", "SqrtWeek", "TxSWeek")
  for (lk in c("logit", "probit")) for (R in c(21L, 1L)) {
    model <- paste("ordered random", lk)
    tag <- sprintf("(R = %d)", R)
    f <- panglm(imps79o ~ TxDrug + SqrtWeek + TxSWeek, data = schizophrenia,
                index = c("id", "Week"), model = "random", family = orderedfamily(lk), R = R)
    # clmm's inner mode search warns on a few individuals at this tolerance
    m <- suppressWarnings(clmm(factor(imps79o) ~ TxDrug + SqrtWeek + TxSWeek + (1 | id),
              data = schizophrenia, link = lk, nAGQ = R,
              control = clmm.control(gradTol = 1e-10, optCtrl = list(rel.tol = 1e-14))))
    ref <- list(c(m$beta, m$alpha)); names(ref) <- sprintf("ordinal::clmm nAGQ=%d", R)
    sref <- function(x) { r <- list(x); names(r) <- names(ref); r }
    out <- c(out, list(
      vcompare(model, "schizophrenia", paste("coef", tag), coef(f), ref, 1e-5),
      vcompare(model, "schizophrenia", paste("sigma_u", tag), f$sigma, sref(m$ST$id[1]), 1e-5),
      vcompare(model, "schizophrenia", paste("logLik", tag), f$loglik, sref(as.numeric(logLik(m))), 1e-7),
      vcompare(model, "schizophrenia", paste("se", tag), sqrt(diag(vcov(f)))[b],
               sref(sqrt(diag(vcov(m)))[b]), 1e-3)
    ))
  }
  do.call(rbind, out)
}

vsave(gaussian_re(inv ~ value + capital, Grunfeld, c("firm", "year"), "Grunfeld"),
      gaussian_re(inv ~ value + capital, grun_u, c("firm", "year"), "Grunfeld unbalanced"),
      gaussian_re(hie_fm, hie, c("id", "year"), "HealthIns (unbalanced)"),
      hausman_u,
      binomial_re(union ~ wage + exper + rural, UnionWage, "UnionWage", agq21 = "logit"),
      probit_converged,
      binomial_re(anymed ~ lcoins + disease + sex + age10 + size + child, hie,
                  "HealthIns (med > 0)"),
      ordered_re(),
      name = "random")
