# Latent class (finite mixture) panel models: panfmm() against flexmix.
#
# flexmix(y ~ x | id, ...) fits the same model: the "| id" grouping keeps
# class membership constant within individual. flexmix is an EM
# implementation; it is started from the panfmm() posterior and run with a
# very small tolerance, so both packages report the same optimum and the
# log-likelihood, class coefficients, class probabilities, concomitant
# coefficients, BIC and posteriors can be compared directly. A separate
# check runs flexmix from its own random starts (stepFlexmix) and confirms
# that it finds no better optimum.
#
# flexmix has no NB2 driver, and its Gaussian driver uses the
# degrees-of-freedom corrected sigma sqrt(RSS_w / (n - p)) in the M step,
# which is not the maximum likelihood estimate. Both are handled with small
# custom flexmix drivers: NB2 M steps by MASS::glm.nb() with prior weights,
# Gaussian M steps by weighted least squares with the ML sigma.
#
# Standard errors: flexmix::refit() inverts the Hessian computed by
# stats::optim(hessian = TRUE) (finite differences of its analytic
# gradient), with the same parameterization as panfmm() for Poisson and
# logit components (class coefficients, multinomial logit coefficients of
# classes 2..k against class 1). The two numerical Hessians agree to about
# 1e-4 relative, so the comparison uses a tolerance of 1e-3.
#
# k = 1 reduces to the pooled GLM: compared with glm(), MASS::glm.nb() and
# panglm(model = "pooling").
#
# Data: a simulated panel of doctor visits (500 patients, 4 years, low and
# high users), and two cross-sections with one row per individual (T = 1):
# AER::DoctorVisits (5190 Australians, 1977-78) and AER::NMES1988 (4406 US
# elderly, the data of Deb and Trivedi 1997).
source("validation/common.R")
suppressPackageStartupMessages({
  library(flexmix)
  library(MASS)
})
out <- list()
add <- function(...) out[[length(out) + 1L]] <<- vcompare(...)
se <- function(v) sqrt(diag(v))
ctrl <- list(tolerance = 1e-14, iter.max = 20000, minprior = 0)

# custom flexmix drivers -------------------------------------------------
FLXMRmlnorm <- function(formula = . ~ .) {
  z <- new("FLXMRglm", weighted = TRUE, formula = formula,
           name = "FLXMRglm:gaussian (ML sigma)", family = "gaussian")
  z@defineComponent <- function(para) {
    predict <- function(x, ...) x %*% para$coef
    logLik <- function(x, y, ...) dnorm(y, predict(x), para$sigma, log = TRUE)
    new("FLXcomponent", parameters = list(coef = para$coef, sigma = para$sigma),
        logLik = logLik, predict = predict, df = para$df)
  }
  z@fit <- function(x, y, w, component) {
    fit <- lm.wfit(x, y, w = w)
    z@defineComponent(list(coef = coef(fit), df = ncol(x) + 1,
                           sigma = sqrt(sum(w * fit$residuals^2) / sum(w))))
  }
  z
}

FLXMRnb2 <- function(formula = . ~ .) {
  z <- new("FLXMR", weighted = TRUE, formula = formula, name = "FLXMR:negbin (MASS::glm.nb)")
  z@defineComponent <- function(para) {
    predict <- function(x, ...) exp(x %*% para$coef)
    logLik <- function(x, y, ...) dnbinom(y, size = para$theta, mu = predict(x), log = TRUE)
    new("FLXcomponent", parameters = list(coef = para$coef, theta = para$theta),
        logLik = logLik, predict = predict, df = para$df)
  }
  z@fit <- function(x, y, w, component) {
    yy <- as.numeric(y)
    fit <- suppressWarnings(MASS::glm.nb(yy ~ x - 1, weights = w,
                                         control = glm.control(epsilon = 1e-13, maxit = 200)))
    cf <- coef(fit); names(cf) <- colnames(x)
    z@defineComponent(list(coef = cf, theta = fit$theta, df = ncol(x) + 1))
  }
  z
}

# compare a panfmm fit with flexmix started from its posterior ----------
compare_flexmix <- function(label, dname, fit, ffit, aux = NULL, conc = FALSE) {
  k <- fit$k
  add(label, dname, "logLik", fit$loglik, list(flexmix = as.numeric(logLik(ffit))), 1e-7)
  add(label, dname, "BIC", BIC(fit), list(flexmix = BIC(ffit)), 1e-7)
  pars <- parameters(ffit)
  cf <- pars[grepl("^coef\\.", rownames(pars)), , drop = FALSE]
  add(label, dname, "class coefficients", as.numeric(fit$classcoef),
      list(flexmix = as.numeric(cf)), 1e-5)
  if (!is.null(aux)) {
    add(label, dname, aux, fit[[aux]], list(flexmix = pars[aux, ]), 1e-5)
  }
  if (conc) {
    pc <- parameters(ffit, which = "concomitant")
    add(label, dname, "concomitant coefficients", as.numeric(fit$prior_coef[, -1]),
        list(flexmix = as.numeric(pc[, -1])), 1e-5)
  } else {
    add(label, dname, "class probabilities", fit$prior_mean, list(flexmix = prior(ffit)), 1e-6)
  }
  pf <- posterior(ffit)
  if (length(ffit@group)) pf <- pf[!duplicated(ffit@group), , drop = FALSE]
  add(label, dname, "posterior", as.numeric(fit$posterior), list(flexmix = as.numeric(pf)), 1e-5)
}

# simulated panel of doctor visits -----------------------------------------
set.seed(20261025)
N <- 500; TT <- 4
d <- data.frame(id = rep(seq_len(N), each = TT), year = rep(seq_len(TT), N))
female <- rbinom(N, 1, 0.5)
d$female <- rep(female, each = TT)
high <- rep(rbinom(N, 1, plogis(-1.2 + 0.8 * female)), each = TT)
d$chronic <- rbinom(N * TT, 1, 0.4)
d$age <- rep(round(runif(N, 40, 80)), each = TT) / 10
d$income <- rnorm(N * TT)
eta <- ifelse(high == 1, 1.2 + 0.6 * d$chronic - 0.15 * d$income,
              -0.6 + 0.3 * d$chronic) + 0.1 * (d$age - 6)
d$visits <- rpois(N * TT, exp(eta))
d$visits_nb <- rnbinom(N * TT, size = ifelse(high == 1, 3, 1.5), mu = exp(eta))
d$anyvisit <- as.numeric(d$visits > 0)
d$logcost <- eta + rnorm(N * TT, sd = ifelse(high == 1, 0.5, 0.8))
idx <- c("id", "year")
dname <- "simulated panel (N = 500, T = 4)"
post_obs <- function(fit, ids) fit$posterior[as.character(ids), , drop = FALSE]

# Poisson, k = 2 and k = 3
for (k in 2:3) {
  a <- panfmm(visits ~ chronic + age + income, d, idx, k = k, family = "poisson",
              nstart = 10, seed = 1)
  b <- flexmix(visits ~ chronic + age + income | id, data = d, k = k,
               cluster = post_obs(a, d$id), model = FLXMRglm(family = "poisson"),
               control = ctrl)
  lab <- sprintf("panfmm poisson k = %d", k)
  compare_flexmix(lab, dname, a, b)
  rf <- refit(b)
  add(lab, dname, "classical SE (vs refit optim Hessian)", se(vcov(a)),
      list(flexmix = sqrt(diag(rf@vcov))), 1e-3)
}

# flexmix from its own random starts reaches the same (not a better) optimum
a2 <- panfmm(visits ~ chronic + age + income, d, idx, k = 2, family = "poisson",
             nstart = 10, seed = 1)
set.seed(7)
sf <- stepFlexmix(visits ~ chronic + age + income | id, data = d, k = 2, nrep = 5,
                  model = FLXMRglm(family = "poisson"), verbose = FALSE,
                  control = list(tolerance = 1e-10, iter.max = 5000))
add("panfmm poisson k = 2", dname, "logLik (flexmix own random starts)", a2$loglik,
    list(flexmix = as.numeric(logLik(sf))), 1e-6)

# Poisson with a concomitant variable (female) for the class probabilities
a <- panfmm(visits ~ chronic + age + income, d, idx, k = 2, family = "poisson",
            concomitant = ~ female, nstart = 10, seed = 1)
b <- flexmix(visits ~ chronic + age + income | id, data = d, k = 2,
             cluster = post_obs(a, d$id), model = FLXMRglm(family = "poisson"),
             concomitant = FLXPmultinom(~ female), control = ctrl)
compare_flexmix("panfmm poisson k = 2, concomitant", dname, a, b, conc = TRUE)
rf <- refit(b)
add("panfmm poisson k = 2, concomitant", dname, "classical SE (vs refit optim Hessian)",
    se(vcov(a)), list(flexmix = sqrt(diag(rf@vcov))), 1e-3)

# NB2, k = 2 (custom driver with MASS::glm.nb M steps)
a <- panfmm(visits_nb ~ chronic + age + income, d, idx, k = 2, family = "negbin",
            nstart = 10, seed = 1)
b <- flexmix(visits_nb ~ chronic + age + income | id, data = d, k = 2,
             cluster = post_obs(a, d$id), model = FLXMRnb2(), control = ctrl)
compare_flexmix("panfmm negbin k = 2", dname, a, b, aux = "theta")

# logit, k = 2
a <- panfmm(anyvisit ~ chronic + age + income, d, idx, k = 2, family = "binomial",
            nstart = 10, seed = 1)
b <- flexmix(cbind(anyvisit, 1 - anyvisit) ~ chronic + age + income | id, data = d, k = 2,
             cluster = post_obs(a, d$id), model = FLXMRglm(family = "binomial"),
             control = ctrl)
compare_flexmix("panfmm binomial k = 2", dname, a, b)
rf <- refit(b)
add("panfmm binomial k = 2", dname, "classical SE (vs refit optim Hessian)", se(vcov(a)),
    list(flexmix = sqrt(diag(rf@vcov))), 1e-3)

# Gaussian, k = 2 (custom driver with the ML sigma)
a <- panfmm(logcost ~ chronic + age + income, d, idx, k = 2, family = "gaussian",
            nstart = 10, seed = 1)
b <- flexmix(logcost ~ chronic + age + income | id, data = d, k = 2,
             cluster = post_obs(a, d$id), model = FLXMRmlnorm(), control = ctrl)
compare_flexmix("panfmm gaussian k = 2", dname, a, b, aux = "sigma")

# k = 1 is the pooled GLM -------------------------------------------------------
lab <- "panfmm k = 1"
a <- panfmm(visits ~ chronic + age + income, d, idx, k = 1, family = "poisson")
g <- glm(visits ~ chronic + age + income, poisson, d)
p <- panglm(visits ~ chronic + age + income, d, idx, model = "pooling", family = "poisson")
add(paste(lab, "poisson"), dname, "coefficients", coef(a),
    list(glm = coef(g), panglm = coef(p)), 1e-8)
add(paste(lab, "poisson"), dname, "logLik", a$loglik, list(glm = as.numeric(logLik(g))), 1e-8)
add(paste(lab, "poisson"), dname, "classical SE", se(vcov(a)), list(glm = se(vcov(g))), 1e-5)

a <- panfmm(visits_nb ~ chronic + age + income, d, idx, k = 1, family = "negbin")
g <- glm.nb(visits_nb ~ chronic + age + income, d, control = glm.control(epsilon = 1e-13))
p <- panglm(visits_nb ~ chronic + age + income, d, idx, model = "pooling", family = "negbin")
add(paste(lab, "negbin"), dname, "coefficients", coef(a)[1:4],
    list(MASS = coef(g), panglm = coef(p)), 1e-6)
add(paste(lab, "negbin"), dname, "theta", a$theta, list(MASS = g$theta, panglm = p$theta), 1e-5)
add(paste(lab, "negbin"), dname, "logLik", a$loglik,
    list(MASS = as.numeric(logLik(g)), panglm = as.numeric(logLik(p))), 1e-8)
add(paste(lab, "negbin"), dname, "theta SE (log scale, delta method)",
    se(vcov(a))[5], list(MASS = g$SE.theta / g$theta), 1e-4)

a <- panfmm(anyvisit ~ chronic + age + income, d, idx, k = 1, family = "binomial")
g <- glm(anyvisit ~ chronic + age + income, binomial, d)
add(paste(lab, "binomial"), dname, "coefficients", coef(a), list(glm = coef(g)), 1e-8)
add(paste(lab, "binomial"), dname, "logLik", a$loglik, list(glm = as.numeric(logLik(g))), 1e-8)
add(paste(lab, "binomial"), dname, "classical SE", se(vcov(a)), list(glm = se(vcov(g))), 1e-5)

a <- panfmm(logcost ~ chronic + age + income, d, idx, k = 1, family = "gaussian")
g <- glm(logcost ~ chronic + age + income, gaussian, d)
add(paste(lab, "gaussian"), dname, "coefficients", coef(a)[1:4], list(glm = coef(g)), 1e-8)
add(paste(lab, "gaussian"), dname, "logLik", a$loglik, list(glm = as.numeric(logLik(g))), 1e-8)
add(paste(lab, "gaussian"), dname, "sigma (ML)", a$sigma,
    list(glm = sqrt(mean(residuals(g)^2))), 1e-8)

# cross-section, T = 1: AER::DoctorVisits -----------------------------------------
data("DoctorVisits", package = "AER")
dv <- DoctorVisits
dv$row <- seq_len(nrow(dv))
dname <- "AER::DoctorVisits (cross-section, id = row)"
f <- visits ~ gender + age + income + illness + reduced + health
a <- panfmm(f, dv, "row", k = 2, family = "poisson", nstart = 5, seed = 1)
b <- flexmix(f, data = dv, k = 2, cluster = a$posterior, model = FLXMRglm(family = "poisson"),
             control = ctrl)
compare_flexmix("panfmm poisson k = 2", dname, a, b)
rf <- refit(b)
add("panfmm poisson k = 2", dname, "classical SE (vs refit optim Hessian)", se(vcov(a)),
    list(flexmix = sqrt(diag(rf@vcov))), 1e-3)

# NB2 on DoctorVisits drives one class to the Poisson boundary (theta -> Inf),
# so the NB2 cross-section check uses AER::NMES1988, the data of Deb and
# Trivedi (1997): physician office visits of 4406 US elderly
data("NMES1988", package = "AER")
nm <- NMES1988
nm$row <- seq_len(nrow(nm))
dname <- "AER::NMES1988 (cross-section, id = row)"
f <- visits ~ hospital + health + chronic + gender + school + insurance
a <- panfmm(f, nm, "row", k = 2, family = "negbin", nstart = 5, seed = 1)
b <- flexmix(f, data = nm, k = 2, cluster = a$posterior, model = FLXMRnb2(), control = ctrl)
compare_flexmix("panfmm negbin k = 2", dname, a, b, aux = "theta")

vsave(do.call(rbind, out), name = "fmm")
