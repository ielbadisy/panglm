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
# Gaussian M steps by weighted least squares with the ML sigma. flexmix's
# binomial driver has the logit link only; probit M steps use a custom driver
# with glm.fit(family = binomial(link = "probit")).
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
# Individual weights: flexmix(weights = ) with integer weights identical
# within individuals, the unweighted fit on replicated individuals, and at
# k = 1 the weighted glm(), MASS::glm.nb() and lm() with sandwich::vcovCL()
# for the cluster covariance.
#
# Offsets: flexmix with FLXMRglm(offset = ) (and the NB2 driver with an
# offset), and at k = 1 glm() / MASS::glm.nb() with the same offset.
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

FLXMRnb2 <- function(formula = . ~ ., offset = NULL) {
  z <- new("FLXMR", weighted = TRUE, formula = formula, name = "FLXMR:negbin (MASS::glm.nb)")
  off <- if (is.null(offset)) 0 else offset
  z@defineComponent <- function(para) {
    predict <- function(x, ...) exp(x %*% para$coef + off)
    logLik <- function(x, y, ...) dnbinom(y, size = para$theta, mu = predict(x), log = TRUE)
    new("FLXcomponent", parameters = list(coef = para$coef, theta = para$theta),
        logLik = logLik, predict = predict, df = para$df)
  }
  z@fit <- function(x, y, w, component) {
    yy <- as.numeric(y)
    oo <- rep_len(off, length(yy))
    fit <- suppressWarnings(MASS::glm.nb(yy ~ x - 1 + offset(oo), weights = w,
                                         control = glm.control(epsilon = 1e-13, maxit = 200)))
    cf <- coef(fit); names(cf) <- colnames(x)
    z@defineComponent(list(coef = cf, theta = fit$theta, df = ncol(x) + 1))
  }
  z
}

FLXMRprobit <- function(formula = . ~ ., offset = NULL) {
  z <- new("FLXMR", weighted = TRUE, formula = formula, name = "FLXMR:binomial (probit)")
  off <- if (is.null(offset)) 0 else offset
  z@defineComponent <- function(para) {
    predict <- function(x, ...) pnorm(x %*% para$coef + off)
    logLik <- function(x, y, ...) dbinom(y, 1, predict(x), log = TRUE)
    new("FLXcomponent", parameters = list(coef = para$coef),
        logLik = logLik, predict = predict, df = para$df)
  }
  z@fit <- function(x, y, w, component) {
    fit <- suppressWarnings(glm.fit(x, as.numeric(y), weights = w,
                                    offset = rep_len(off, nrow(x)),
                                    family = binomial(link = "probit"),
                                    control = glm.control(epsilon = 1e-13, maxit = 200)))
    z@defineComponent(list(coef = fit$coefficients, df = ncol(x)))
  }
  z
}

# compare a panfmm fit with flexmix started from its posterior ----------
compare_flexmix <- function(label, dname, fit, ffit, aux = NULL, conc = FALSE, bic = TRUE,
                            ctol = 1e-5) {
  k <- fit$k
  add(label, dname, "logLik", fit$loglik, list(flexmix = as.numeric(logLik(ffit))), 1e-7)
  if (bic) add(label, dname, "BIC", BIC(fit), list(flexmix = BIC(ffit)), 1e-7)
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
        list(flexmix = as.numeric(pc[, -1])), ctol)
  } else {
    add(label, dname, "class probabilities", fit$prior_mean, list(flexmix = prior(ffit)), 1e-6)
  }
  pf <- posterior(ffit)
  if (length(ffit@group)) pf <- pf[!duplicated(ffit@group), , drop = FALSE]
  add(label, dname, "posterior", as.numeric(fit$posterior), list(flexmix = as.numeric(pf)), ctol)
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

# individual (sampling) weights ----------------------------------------------
# flexmix(weights = ) takes integer weights identical within groups and
# multiplies each group's log-likelihood by its weight, the panfmm()
# convention. flexmix computes BIC with n = sum of the weights, panfmm() with
# n = number of observations, so BIC is not compared here. With integer
# weights the fit equals the unweighted fit on the data with each
# individual replicated w_i times (checked below), and at k = 1 the
# weighted fit is the weighted pooled GLM, with non-integer weights.
wid <- sample(1:3, N, replace = TRUE)
d$w <- rep(wid, each = TT)
lab <- "panfmm weighted"
a <- panfmm(visits ~ chronic + age + income, d, idx, k = 2, family = "poisson",
            weights = d$w, nstart = 10, seed = 1)
b <- flexmix(visits ~ chronic + age + income | id, data = d, k = 2, weights = d$w,
             cluster = post_obs(a, d$id), model = FLXMRglm(family = "poisson"),
             control = ctrl)
compare_flexmix(paste(lab, "poisson k = 2"), dname, a, b, bic = FALSE)
rf <- refit(b)
add(paste(lab, "poisson k = 2"), dname, "classical SE (vs refit optim Hessian)", se(vcov(a)),
    list(flexmix = sqrt(diag(rf@vcov))), 1e-3)
rep_rows <- rep(seq_len(nrow(d)), d$w)
dr <- d[rep_rows, ]
dr$rid <- paste(dr$id, ave(dr$year, dr$id, dr$year, FUN = seq_along), sep = ".")
ar <- panfmm(visits ~ chronic + age + income, dr, c("rid", "year"), k = 2,
             family = "poisson", nstart = 10, seed = 1)
add(paste(lab, "poisson k = 2"), dname, "coefficients (vs replicated individuals)",
    coef(a), list(panfmm_replicated = coef(ar)), 1e-6)
add(paste(lab, "poisson k = 2"), dname, "classical SE (vs replicated individuals)",
    se(vcov(a)), list(panfmm_replicated = se(vcov(ar))), 1e-5)

a <- panfmm(visits ~ chronic + age + income, d, idx, k = 2, family = "poisson",
            concomitant = ~ female, weights = d$w, nstart = 10, seed = 1)
b <- flexmix(visits ~ chronic + age + income | id, data = d, k = 2, weights = d$w,
             cluster = post_obs(a, d$id), model = FLXMRglm(family = "poisson"),
             concomitant = FLXPmultinom(~ female), control = ctrl)
# flexmix stops here after 5 EM iterations with a log-likelihood 4e-10
# below the optimum (its nnet::multinom concomitant M step has relative
# tolerance 1e-8), so the concomitant coefficients and posteriors agree to
# 1e-3 only; the replicated-individuals fit, at the optimum, agrees exactly
compare_flexmix(paste(lab, "poisson k = 2, concomitant"), dname, a, b, conc = TRUE,
                bic = FALSE, ctol = 1e-3)
ar <- panfmm(visits ~ chronic + age + income, dr, c("rid", "year"), k = 2,
             family = "poisson", concomitant = ~ female, nstart = 10, seed = 1)
add(paste(lab, "poisson k = 2, concomitant"), dname,
    "coefficients (vs replicated individuals)", coef(a),
    list(panfmm_replicated = coef(ar)), 1e-6)
add(paste(lab, "poisson k = 2, concomitant"), dname,
    "classical SE (vs replicated individuals)", se(vcov(a)),
    list(panfmm_replicated = se(vcov(ar))), 1e-5)

a <- panfmm(visits_nb ~ chronic + age + income, d, idx, k = 2, family = "negbin",
            weights = d$w, nstart = 10, seed = 1)
b <- flexmix(visits_nb ~ chronic + age + income | id, data = d, k = 2, weights = d$w,
             cluster = post_obs(a, d$id), model = FLXMRnb2(), control = ctrl)
compare_flexmix(paste(lab, "negbin k = 2"), dname, a, b, aux = "theta", bic = FALSE)

a <- panfmm(anyvisit ~ chronic + age + income, d, idx, k = 2, family = "binomial",
            weights = d$w, nstart = 10, seed = 1)
b <- flexmix(cbind(anyvisit, 1 - anyvisit) ~ chronic + age + income | id, data = d, k = 2,
             weights = d$w, cluster = post_obs(a, d$id),
             model = FLXMRglm(family = "binomial"), control = ctrl)
compare_flexmix(paste(lab, "binomial k = 2"), dname, a, b, bic = FALSE)

a <- panfmm(logcost ~ chronic + age + income, d, idx, k = 2, family = "gaussian",
            weights = d$w, nstart = 10, seed = 1)
b <- flexmix(logcost ~ chronic + age + income | id, data = d, k = 2, weights = d$w,
             cluster = post_obs(a, d$id), model = FLXMRmlnorm(), control = ctrl)
compare_flexmix(paste(lab, "gaussian k = 2"), dname, a, b, aux = "sigma", bic = FALSE)

# k = 1 with non-integer weights: the weighted pooled GLM; the cluster
# sandwich equals sandwich::vcovCL(type = "HC0") clustered by individual
# (factor G / (G - 1))
d$sw <- rep(runif(N, 0.2, 3), each = TT)
lab <- "panfmm weighted k = 1"
a <- panfmm(visits ~ chronic + age + income, d, idx, k = 1, family = "poisson",
            weights = d$sw)
g <- glm(visits ~ chronic + age + income, poisson, d, weights = sw,
         control = glm.control(epsilon = 1e-14))
add(paste(lab, "poisson"), dname, "coefficients", coef(a), list(glm = coef(g)), 1e-8)
add(paste(lab, "poisson"), dname, "logLik", a$loglik, list(glm = as.numeric(logLik(g))), 1e-8)
add(paste(lab, "poisson"), dname, "classical SE", se(vcov(a)), list(glm = se(vcov(g))), 1e-5)
add(paste(lab, "poisson"), dname, "cluster SE", se(vcov(a, type = "cluster")),
    list(sandwich = se(sandwich::vcovCL(g, cluster = ~ id, type = "HC0"))), 1e-5)

a <- panfmm(visits_nb ~ chronic + age + income, d, idx, k = 1, family = "negbin",
            weights = d$sw)
g <- suppressWarnings(glm.nb(visits_nb ~ chronic + age + income, d, weights = sw,
                             control = glm.control(epsilon = 1e-13)))
add(paste(lab, "negbin"), dname, "coefficients", coef(a)[1:4], list(MASS = coef(g)), 1e-6)
add(paste(lab, "negbin"), dname, "theta", a$theta, list(MASS = g$theta), 1e-5)
add(paste(lab, "negbin"), dname, "logLik", a$loglik, list(MASS = as.numeric(logLik(g))), 1e-8)

a <- panfmm(anyvisit ~ chronic + age + income, d, idx, k = 1, family = "binomial",
            weights = d$sw)
g <- suppressWarnings(glm(anyvisit ~ chronic + age + income, binomial, d, weights = sw,
                          control = glm.control(epsilon = 1e-14)))
add(paste(lab, "binomial"), dname, "coefficients", coef(a), list(glm = coef(g)), 1e-8)
add(paste(lab, "binomial"), dname, "classical SE", se(vcov(a)), list(glm = se(vcov(g))), 1e-5)
add(paste(lab, "binomial"), dname, "cluster SE", se(vcov(a, type = "cluster")),
    list(sandwich = se(sandwich::vcovCL(g, cluster = ~ id, type = "HC0"))), 1e-5)

a <- panfmm(logcost ~ chronic + age + income, d, idx, k = 1, family = "gaussian",
            weights = d$sw)
g <- lm(logcost ~ chronic + age + income, d, weights = sw)
add(paste(lab, "gaussian"), dname, "coefficients", coef(a)[1:4], list(lm = coef(g)), 1e-8)
add(paste(lab, "gaussian"), dname, "sigma (weighted ML)", a$sigma,
    list(lm = sqrt(sum(d$sw * residuals(g)^2) / sum(d$sw))), 1e-8)

# offsets ----------------------------------------------------------------------
# Visits over observation periods of 3 to 12 months, offset(log(months / 12)).
# flexmix silently ignores offset() terms in its formula; the offset is
# given to FLXMRglm(offset = ) and to the custom NB2 driver instead.
d$months <- sample(3:12, N * TT, replace = TRUE)
d$lexp <- log(d$months / 12)
d$visits_e <- rpois(N * TT, exp(eta + d$lexp))
d$visits_enb <- rnbinom(N * TT, size = ifelse(high == 1, 3, 1.5), mu = exp(eta + d$lexp))
lab <- "panfmm offset"
fo <- visits_e ~ chronic + age + income + offset(log(months / 12))
a <- panfmm(fo, d, idx, k = 2, family = "poisson", nstart = 10, seed = 1)
a2 <- panfmm(visits_e ~ chronic + age + income, d, idx, k = 2, family = "poisson",
             offset = d$lexp, nstart = 10, seed = 1)
add(paste(lab, "poisson k = 2"), dname, "coefficients (offset argument vs formula)",
    coef(a2), list(panfmm_formula = coef(a)), 1e-10)
b <- flexmix(visits_e ~ chronic + age + income | id, data = d, k = 2,
             cluster = post_obs(a, d$id),
             model = FLXMRglm(family = "poisson", offset = d$lexp), control = ctrl)
compare_flexmix(paste(lab, "poisson k = 2"), dname, a, b)
rf <- refit(b)
add(paste(lab, "poisson k = 2"), dname, "classical SE (vs refit optim Hessian)", se(vcov(a)),
    list(flexmix = sqrt(diag(rf@vcov))), 1e-3)

a <- panfmm(visits_e ~ chronic + age + income + offset(log(months / 12)), d, idx, k = 2,
            family = "poisson", weights = d$w, nstart = 10, seed = 1)
b <- flexmix(visits_e ~ chronic + age + income | id, data = d, k = 2, weights = d$w,
             cluster = post_obs(a, d$id),
             model = FLXMRglm(family = "poisson", offset = d$lexp), control = ctrl)
compare_flexmix(paste(lab, "poisson k = 2, weighted"), dname, a, b, bic = FALSE)

a <- panfmm(visits_enb ~ chronic + age + income + offset(log(months / 12)), d, idx, k = 2,
            family = "negbin", nstart = 10, seed = 1)
b <- flexmix(visits_enb ~ chronic + age + income | id, data = d, k = 2,
             cluster = post_obs(a, d$id), model = FLXMRnb2(offset = d$lexp), control = ctrl)
compare_flexmix(paste(lab, "negbin k = 2"), dname, a, b, aux = "theta")

lab <- "panfmm offset k = 1"
a <- panfmm(fo, d, idx, k = 1, family = "poisson")
g <- glm(fo, poisson, d, control = glm.control(epsilon = 1e-14))
add(paste(lab, "poisson"), dname, "coefficients", coef(a), list(glm = coef(g)), 1e-8)
add(paste(lab, "poisson"), dname, "logLik", a$loglik, list(glm = as.numeric(logLik(g))), 1e-8)
add(paste(lab, "poisson"), dname, "classical SE", se(vcov(a)), list(glm = se(vcov(g))), 1e-5)
add(paste(lab, "poisson"), dname, "fitted values", predict(a), list(glm = fitted(g)), 1e-8)
add(paste(lab, "poisson"), dname, "predict(newdata)", predict(a, newdata = d[1:20, ]),
    list(glm = predict(g, newdata = d[1:20, ], type = "response")), 1e-8)

a <- panfmm(visits_enb ~ chronic + age + income + offset(log(months / 12)), d, idx, k = 1,
            family = "negbin")
g <- glm.nb(visits_enb ~ chronic + age + income + offset(log(months / 12)), d,
            control = glm.control(epsilon = 1e-13))
add(paste(lab, "negbin"), dname, "coefficients", coef(a)[1:4], list(MASS = coef(g)), 1e-6)
add(paste(lab, "negbin"), dname, "theta", a$theta, list(MASS = g$theta), 1e-5)
add(paste(lab, "negbin"), dname, "logLik", a$loglik, list(MASS = as.numeric(logLik(g))), 1e-8)

a <- panfmm(anyvisit ~ chronic + age + income, d, idx, k = 1, family = "binomial",
            offset = 0.5 * d$income)
g <- glm(anyvisit ~ chronic + age + income, binomial, d, offset = 0.5 * income,
         control = glm.control(epsilon = 1e-14))
add(paste(lab, "binomial"), dname, "coefficients", coef(a), list(glm = coef(g)), 1e-8)
add(paste(lab, "binomial"), dname, "logLik", a$loglik, list(glm = as.numeric(logLik(g))), 1e-8)
add(paste(lab, "binomial"), dname, "classical SE", se(vcov(a)), list(glm = se(vcov(g))), 1e-5)

a <- panfmm(logcost ~ chronic + age + income + offset(lexp), d, idx, k = 1,
            family = "gaussian")
g <- glm(logcost ~ chronic + age + income + offset(lexp), gaussian, d)
add(paste(lab, "gaussian"), dname, "coefficients", coef(a)[1:4], list(glm = coef(g)), 1e-8)
add(paste(lab, "gaussian"), dname, "logLik", a$loglik, list(glm = as.numeric(logLik(g))), 1e-8)

# probit components -----------------------------------------------------------
# k = 2 against flexmix with the probit driver (unweighted, and weighted with
# an offset); k = 1 against glm(binomial(link = "probit")) and
# panglm(model = "pooling", family = "probit"). For the non-canonical
# probit link glm() reports standard errors from the expected information
# (Fisher scoring), which differ from the observed information by about
# 1e-4 here; panfmm() inverts the observed Hessian, as does
# sampleSelection::probit() (maxLik, analytic Hessian), the reference for
# the classical SE and, through sandwich::vcovCL(type = "HC0"), for the
# cluster SE.
lab <- "panfmm probit"
a <- panfmm(anyvisit ~ chronic + age + income, d, idx, k = 2, family = "probit",
            nstart = 10, seed = 1)
b <- flexmix(anyvisit ~ chronic + age + income | id, data = d, k = 2,
             cluster = post_obs(a, d$id), model = FLXMRprobit(), control = ctrl)
compare_flexmix(paste(lab, "k = 2"), dname, a, b)
a2 <- panfmm(anyvisit ~ chronic + age + income, d, idx, k = 2,
             family = binomial_family("probit"), nstart = 10, seed = 1)
add(paste(lab, "k = 2"), dname, "coefficients (binomial_family(\"probit\") vs \"probit\")",
    coef(a2), list(panfmm_probit = coef(a)), 1e-12)

a <- panfmm(anyvisit ~ chronic + age + income, d, idx, k = 2, family = "probit",
            weights = d$w, offset = 0.3 * d$female, nstart = 10, seed = 1)
b <- flexmix(anyvisit ~ chronic + age + income | id, data = d, k = 2, weights = d$w,
             cluster = post_obs(a, d$id), model = FLXMRprobit(offset = 0.3 * d$female),
             control = ctrl)
compare_flexmix(paste(lab, "k = 2, weighted, offset"), dname, a, b, bic = FALSE)

lab <- "panfmm probit k = 1"
a <- panfmm(anyvisit ~ chronic + age + income, d, idx, k = 1, family = "probit")
g <- glm(anyvisit ~ chronic + age + income, binomial(link = "probit"), d,
         control = glm.control(epsilon = 1e-14))
p <- panglm(anyvisit ~ chronic + age + income, d, idx, model = "pooling", family = "probit")
add(lab, dname, "coefficients", coef(a), list(glm = coef(g), panglm = coef(p)), 1e-8)
add(lab, dname, "logLik", a$loglik,
    list(glm = as.numeric(logLik(g)), panglm = as.numeric(logLik(p))), 1e-8)
sp <- sampleSelection::probit(anyvisit ~ chronic + age + income, data = d)
add(lab, dname, "coefficients", coef(a), list(sampleSelection = coef(sp)), 1e-6)
add(lab, dname, "classical SE (observed information)", se(vcov(a)),
    list(sampleSelection = se(vcov(sp))), 1e-6)
add(lab, dname, "cluster SE", se(vcov(a, type = "cluster")),
    list(sandwich = se(sandwich::vcovCL(sp, cluster = d$id, type = "HC0"))), 1e-6)
add(lab, dname, "fitted values", predict(a), list(glm = fitted(g)), 1e-8)

a <- panfmm(anyvisit ~ chronic + age + income, d, idx, k = 1, family = "probit",
            weights = d$sw)
g <- suppressWarnings(glm(anyvisit ~ chronic + age + income, binomial(link = "probit"), d,
                          weights = sw, control = glm.control(epsilon = 1e-14)))
sp <- sampleSelection::probit(anyvisit ~ chronic + age + income, data = d, weights = d$sw)
add(paste(lab, "weighted"), dname, "coefficients", coef(a), list(glm = coef(g)), 1e-8)
add(paste(lab, "weighted"), dname, "classical SE (observed information)", se(vcov(a)),
    list(sampleSelection = se(vcov(sp))), 1e-6)
add(paste(lab, "weighted"), dname, "cluster SE", se(vcov(a, type = "cluster")),
    list(sandwich = se(sandwich::vcovCL(sp, cluster = d$id, type = "HC0"))), 1e-6)

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
