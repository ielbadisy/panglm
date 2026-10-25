# Validation of the pooled panel GLMs (model = "pooling").
#
# Every family is fit on the same data with panglm and with stats::lm /
# stats::glm, MASS::glm.nb, plm (linear pooling) and fixest (feols /
# feglm / fenegbin without fixed effects). Compared: coefficients, the
# classical, HC1 and cluster-by-individual standard errors, and the
# log-likelihood. All these packages maximize the same likelihood, so the
# match is demanded to 1e-6 (relative to max(1, |reference|)).
#
# Note on convergence: panglm's IRLS stops on a relative deviance change
# below tol (1e-10). The deviance is quadratic in the coefficient error, so
# with the default tol the probit coefficients are accurate to about 1e-6
# only (glm with epsilon = 1e-12 is closer to the score root). The default
# is validated here; tol = 1e-14 brings the probit gap to about 1e-9.
#
# Small-sample conventions matched by the references:
#   HC1:     n / (n - k)                          (sandwich::vcovHC type "HC1",
#                                                  fixest vcov = "hetero")
#   cluster: G / (G - 1) * (n - 1) / (n - k)      (sandwich::vcovCL type "HC1",
#                                                  plm::vcovHC arellano "sss",
#                                                  fixest vcov = ~id)
#
# Run from the package root: Rscript validation/validate-pooled.R
source("validation/common.R")
suppressPackageStartupMessages({
  library(sandwich)
  library(plm)
  library(fixest)
  library(MASS)
})
setFixest_nthreads(1)

data(copd, package = "panglm")
copd$id <- as.integer(copd$id)
idx <- c("id", "visit")
dname <- "copd (60 x 5)"
tol <- 1e-6
se <- function(v) sqrt(diag(v))
out <- list()

# compare coefficients and the three standard errors of a panglm fit with
# per-reference lists of coefficients and standard errors
compare_fit <- function(model, fit, coefs, ses, loglik = NULL) {
  r <- list(vcompare(model, dname, "coef", coef(fit), coefs, tol))
  for (type in names(ses)) {
    r[[length(r) + 1L]] <- vcompare(model, dname, paste0("se_", type),
                                    se(vcov(fit, type = type)), ses[[type]], tol)
  }
  if (!is.null(loglik)) {
    r[[length(r) + 1L]] <- vcompare(model, dname, "logLik", as.numeric(logLik(fit)),
                                    loglik, tol)
  }
  do.call(rbind, r)
}

# ---- gaussian: pooled OLS -------------------------------------------------
f_g <- fev1 ~ treatment + age + smoker + crp
p <- panglm(f_g, copd, idx, model = "pooling", family = "gaussian")
m <- lm(f_g, copd)
pl <- plm(f_g, pdata.frame(copd, idx), model = "pooling")
fx <- feols(f_g, copd)
out$gaussian <- compare_fit(
  "pooled gaussian", p,
  coefs = list(lm = coef(m), plm = coef(pl), fixest = coef(fx)),
  ses = list(
    classical = list(lm = se(vcov(m)), plm = se(vcov(pl)), fixest = se(vcov(fx, vcov = "iid"))),
    HC1 = list(sandwich = se(vcovHC(m, type = "HC1")),
               plm = se(plm::vcovHC(pl, method = "white1", type = "HC1")),
               fixest = se(vcov(fx, vcov = "hetero"))),
    cluster = list(sandwich = se(vcovCL(m, cluster = ~id, type = "HC1")),
                   plm = se(plm::vcovHC(pl, method = "arellano", type = "sss")),
                   fixest = se(vcov(fx, vcov = ~id)))
  ),
  loglik = list(lm = as.numeric(logLik(m)), fixest = as.numeric(logLik(fx)))
)

# ---- GLM families vs glm, fixest::feglm and sandwich ------------------------
glm_case <- function(label, formula, family, glm_family) {
  p <- panglm(formula, copd, idx, model = "pooling", family = family)
  g <- glm(formula, family = glm_family, data = copd, control = glm.control(epsilon = 1e-12, maxit = 100))
  fx <- feglm(formula, data = copd, family = glm_family, glm.tol = 1e-12, glm.iter = 100)
  compare_fit(
    label, p,
    coefs = list(glm = coef(g), fixest = coef(fx)),
    ses = list(
      classical = list(glm = se(vcov(g)), fixest = se(vcov(fx, vcov = "iid", ssc = ssc(K.adj = FALSE)))),
      HC1 = list(sandwich = se(vcovHC(g, type = "HC1")), fixest = se(vcov(fx, vcov = "hetero"))),
      cluster = list(sandwich = se(vcovCL(g, cluster = ~id, type = "HC1")),
                     fixest = se(vcov(fx, vcov = ~id)))
    ),
    loglik = list(glm = as.numeric(logLik(g)))
  )
}
f_c <- exacerbations ~ treatment + age + smoker + crp
f_b <- hospitalized ~ treatment + age + smoker + crp
out$poisson <- glm_case("pooled poisson", f_c, "poisson", poisson())
out$logit <- glm_case("pooled logit", f_b, "binomial", binomial("logit"))
out$probit <- glm_case("pooled probit", f_b, "probit", binomial("probit"))
out$gamma <- glm_case("pooled gamma (log)", f_g, "gamma", Gamma("log"))

# ---- negative binomial (NB2) vs MASS::glm.nb and fixest::fenegbin -----------
# Overdispersed counts: a gamma individual frailty on the copd design.
set.seed(20261025)
copd_nb <- copd
frailty <- rgamma(60, shape = 2, rate = 2)[copd$id]
copd_nb$count <- rpois(nrow(copd), frailty * exp(-0.5 - 0.6 * copd$treatment +
                                                  0.01 * copd$age + 0.15 * copd$crp))
p <- panglm(count ~ treatment + age + smoker + crp, copd_nb, idx, model = "pooling", family = "negbin")
nb <- glm.nb(count ~ treatment + age + smoker + crp, data = copd_nb,
             control = glm.control(epsilon = 1e-12, maxit = 100))
fx <- fenegbin(count ~ treatment + age + smoker + crp, data = copd_nb)
out$negbin <- rbind(
  compare_fit(
    "pooled negbin", p,
    coefs = list(MASS = coef(nb), fixest = coef(fx)),
    ses = list(
      classical = list(MASS = se(vcov(nb))),
      HC1 = list(sandwich = se(vcovHC(nb, type = "HC1"))),
      cluster = list(sandwich = se(vcovCL(nb, cluster = ~id, type = "HC1")))
    ),
    loglik = list(MASS = as.numeric(logLik(nb)), fixest = as.numeric(logLik(fx)))
  ),
  vcompare("pooled negbin", dname, "theta", p$theta,
           list(MASS = nb$theta, fixest = fx$theta), tol)
)

do.call(vsave, c(unname(out), name = "pooled"))
