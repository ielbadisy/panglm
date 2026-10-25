# Random slopes: panglm(model = "random", random = ~ time) vs lme4, nlme +
# clubSandwich, glmmTMB and GLMMadaptive.
# Gaussian: exact ML vs lme4::lmer(REML = FALSE) with one and two slopes,
# with cre = TRUE, and the cluster-robust covariance vs clubSandwich CR0 on
# nlme::lme. Poisson, logit and probit: the Laplace approximation (R = 1) vs
# glmmTMB (and lme4::glmer for logit and Poisson), adaptive Gauss-Hermite
# quadrature (R = 7) vs GLMMadaptive, the quadrature likelihood vs a
# brute-force dense grid, and population-averaged partial effects by hand.
# Two slopes (a three-dimensional random effect): Laplace vs glmmTMB and
# adaptive quadrature (R = 7, 343 nodes) vs GLMMadaptive. A singular
# (boundary) GLM fit: finite standard errors vs glmmTMB. Negative binomial
# (NB2) and gamma (log link) with lognormal random effects vs glmmTMB
# (Laplace) and, for NB2, GLMMadaptive (R = 7).
source("validation/common.R")
suppressPackageStartupMessages({
  library(lme4); library(glmmTMB); library(GLMMadaptive); library(nlme); library(clubSandwich)
})

vech <- function(m) m[lower.tri(m, diag = TRUE)]
idx <- c("id", "time")

# ---- Gaussian: patient FEV1 trajectories, random intercept and slope
set.seed(20)
G <- 150; T <- 6
pat <- data.frame(id = rep(seq_len(G), each = T), time = rep(0:(T - 1), G))
pat$age <- rep(round(rnorm(G, 65, 8)), each = T)
pat$smoker <- rep(rbinom(G, 1, 0.3), each = T)
pat$crp <- rgamma(G * T, 2, 0.5)
Lg <- matrix(c(0.45, 0, -0.02, 0.06), 2, 2, byrow = TRUE)
bg <- matrix(rnorm(G * 2), G) %*% t(Lg)
pat$fev1 <- 3.5 - 0.08 * pat$time - 0.015 * (pat$age - 65) - 0.25 * pat$smoker -
  0.02 * pat$crp + bg[pat$id, 1] + bg[pat$id, 2] * pat$time + rnorm(G * T, 0, 0.15)
# second slope variable (time-varying CRP) for a three-dimensional effect
pat$fev1b <- pat$fev1 + rep(rnorm(G, 0, 0.03), each = T) * pat$crp

gaussian_block <- function(fm, random, lmer_fm, label, data = pat) {
  f <- panglm(fm, data = data, index = idx, model = "random", random = random)
  m <- lmer(lmer_fm, data = data, REML = FALSE,
            control = lmerControl(optimizer = "bobyqa",
                                  optCtrl = list(maxfun = 1e5, rhobeg = 1e-2, rhoend = 1e-12)))
  vc <- VarCorr(m)$id
  rbind(
    vcompare("slopes gaussian", label, "coef", coef(f), list(lme4 = fixef(m)), 1e-5),
    vcompare("slopes gaussian", label, "se", sqrt(diag(vcov(f))),
             list(lme4 = sqrt(diag(as.matrix(vcov(m))))), 1e-5),
    vcompare("slopes gaussian", label, "varcomp (vech Sigma)", vech(f$varcomp),
             list(lme4 = vech(vc)), 1e-5),
    vcompare("slopes gaussian", label, "correlation", vech(cov2cor(f$varcomp))[-1],
             list(lme4 = vech(attr(vc, "correlation"))[-1]), 1e-4),
    vcompare("slopes gaussian", label, "residual variance", c(sigma2 = f$sigma_v2),
             list(lme4 = sigma(m)^2), 1e-5),
    vcompare("slopes gaussian", label, "loglik", c(loglik = f$loglik),
             list(lme4 = as.numeric(logLik(m))), 1e-6),
    vcompare("slopes gaussian", label, "ranef (BLUP)", c(f$ranef),
             list(lme4 = c(as.matrix(ranef(m)$id))), 1e-5)
  )
}

g1 <- gaussian_block(fev1 ~ time + age + smoker + crp, ~ time,
                     fev1 ~ time + age + smoker + crp + (1 + time | id), "sim, 1 slope")
g2 <- gaussian_block(fev1b ~ time + age + smoker + crp, ~ time + crp,
                     fev1b ~ time + age + smoker + crp + (1 + time + crp | id), "sim, 2 slopes")

# cre = TRUE: Mundlak mean of crp added as a regressor
pat$mean_crp <- ave(pat$crp, pat$id)
fc <- suppressWarnings(panglm(fev1 ~ time + age + smoker + crp, data = pat, index = idx, model = "random",
             random = ~ time, cre = TRUE)) # mean_time is aliased (balanced panel)
mc <- lmer(fev1 ~ time + age + smoker + crp + mean_crp + (1 + time | id), data = pat,
           REML = FALSE, control = lmerControl(optimizer = "bobyqa",
                                               optCtrl = list(maxfun = 1e5, rhoend = 1e-12)))
g3 <- rbind(
  vcompare("slopes gaussian", "sim, cre = TRUE", "coef", coef(fc), list(lme4 = fixef(mc)), 1e-5),
  vcompare("slopes gaussian", "sim, cre = TRUE", "se", sqrt(diag(vcov(fc))),
           list(lme4 = sqrt(diag(as.matrix(vcov(mc))))), 1e-5),
  vcompare("slopes gaussian", "sim, cre = TRUE", "loglik", c(loglik = fc$loglik),
           list(lme4 = as.numeric(logLik(mc))), 1e-6)
)

# cluster-robust covariance: CR0 sandwich of nlme::lme (ML) times G / (G - 1)
f1 <- panglm(fev1 ~ time + age + smoker + crp, data = pat, index = idx, model = "random",
             random = ~ time, vcov = "cluster")
ml <- lme(fev1 ~ time + age + smoker + crp, random = ~ time | id, data = pat, method = "ML",
          control = lmeControl(tolerance = 1e-12, msTol = 1e-12, maxIter = 500, msMaxIter = 500))
cr0 <- as.matrix(vcovCR(ml, cluster = pat$id, type = "CR0")) * G / (G - 1)
g4 <- rbind(
  vcompare("slopes gaussian", "sim, 1 slope", "se_cluster", sqrt(diag(vcov(f1))),
           list(clubSandwich = sqrt(diag(cr0))), 1e-5),
  vcompare("slopes gaussian", "sim, 1 slope", "coef", coef(f1), list(nlme = fixef(ml)), 1e-5),
  vcompare("slopes gaussian", "sim, 1 slope", "predict(link)",
           predict(f1, newdata = pat[1:12, ], type = "link"),
           list(nlme = predict(ml, newdata = pat[1:12, ], level = 0)), 1e-5)
)

# copd: FEV1 per visit (the slope variance is at the boundary, a singular fit)
data(copd)
fcopd <- panglm(fev1 ~ visit + treatment + age + smoker, data = copd, index = c("id", "visit"),
                model = "random", random = ~ visit)
mcopd <- suppressMessages(lmer(fev1 ~ visit + treatment + age + smoker + (1 + visit | id),
                               data = copd, REML = FALSE))
g5 <- rbind(
  vcompare("slopes gaussian", "copd (boundary)", "loglik", c(loglik = fcopd$loglik),
           list(lme4 = as.numeric(logLik(mcopd))), 1e-6),
  vcompare("slopes gaussian", "copd (boundary)", "coef", coef(fcopd),
           list(lme4 = fixef(mcopd)), 1e-4)
)

# ---- Poisson, logit, probit: one random slope
set.seed(3)
G2 <- 200; T2 <- 6
d <- data.frame(id = rep(seq_len(G2), each = T2), time = rep(0:(T2 - 1), G2))
d$x <- rnorm(G2 * T2)
Lb <- matrix(c(0.7, 0, -0.15, 0.25), 2, 2, byrow = TRUE)
b <- matrix(rnorm(G2 * 2), G2) %*% t(Lb)
lin <- 0.1 * d$time + 0.4 * d$x + b[d$id, 1] + b[d$id, 2] * d$time
d$yb <- rbinom(G2 * T2, 1, plogis(-0.5 + lin))
d$yp <- rbinom(G2 * T2, 1, pnorm(-0.3 + 0.6 * lin))
d$yc <- rpois(G2 * T2, exp(0.3 + lin))

tmb_vech <- function(m) vech(VarCorr(m)$cond$id)

glm_block <- function(y, family, glm_family, label, with_glmer) {
  fm <- stats::as.formula(paste(y, "~ time + x"))
  fr <- stats::as.formula(paste(y, "~ time + x + (1 + time | id)"))
  f <- panglm(fm, data = d, index = idx, model = "random", family = family,
              random = ~ time, R = 1)
  tm <- glmmTMB(fr, data = d, family = glm_family)
  out <- rbind(
    vcompare("slopes Laplace", label, "coef", coef(f), list(glmmTMB = fixef(tm)$cond), 1e-4),
    vcompare("slopes Laplace", label, "se", sqrt(diag(vcov(f))),
             list(glmmTMB = sqrt(diag(vcov(tm)$cond))), 1e-4),
    vcompare("slopes Laplace", label, "varcomp (vech Sigma)", vech(f$varcomp),
             list(glmmTMB = tmb_vech(tm)), 1e-4),
    vcompare("slopes Laplace", label, "loglik", c(loglik = f$loglik),
             list(glmmTMB = as.numeric(logLik(tm))), 1e-6),
    vcompare("slopes Laplace", label, "ranef (modes)", c(f$ranef),
             list(glmmTMB = c(as.matrix(ranef(tm)$cond$id))), 1e-3)
  )
  if (with_glmer) {
    # glmer stops earlier (bobyqa on a finite-difference deviance): looser tolerance
    m <- glmer(fr, data = d, family = glm_family, nAGQ = 1,
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
    out <- rbind(out,
      vcompare("slopes Laplace", label, "coef", coef(f), list(lme4 = fixef(m)), 1e-3),
      vcompare("slopes Laplace", label, "loglik", c(loglik = f$loglik),
               list(lme4 = as.numeric(logLik(m))), 1e-4))
  }
  out
}

l1 <- glm_block("yc", "poisson", poisson(), "sim Poisson", TRUE)
l2 <- glm_block("yb", "binomial", binomial(), "sim logit", TRUE)
l3 <- glm_block("yp", "probit", binomial("probit"), "sim probit", FALSE)

# adaptive Gauss-Hermite quadrature (7 nodes per dimension) vs GLMMadaptive
aghq_block <- function(y, family, glm_family, label) {
  fm <- stats::as.formula(paste(y, "~ time + x"))
  f <- panglm(fm, data = d, index = idx, model = "random", family = family, random = ~ time)
  ga <- mixed_model(fm, random = ~ time | id, data = d, family = glm_family, nAGQ = 7,
                    control = list(iter_EM = 0, update_GH_every = 1,
                                   tol1 = 1e-10, tol2 = 1e-10, tol3 = 1e-12))
  rbind(
    vcompare("slopes AGHQ R = 7", label, "coef", coef(f), list(GLMMadaptive = fixef(ga)), 1e-4),
    vcompare("slopes AGHQ R = 7", label, "se", sqrt(diag(vcov(f))),
             list(GLMMadaptive = sqrt(diag(vcov(ga)))[1:3]), 1e-3),
    vcompare("slopes AGHQ R = 7", label, "varcomp (vech Sigma)", vech(f$varcomp),
             list(GLMMadaptive = vech(ga$D)), 1e-4),
    vcompare("slopes AGHQ R = 7", label, "loglik", c(loglik = f$loglik),
             list(GLMMadaptive = as.numeric(logLik(ga))), 1e-5)
  )
}
a1 <- aghq_block("yc", "poisson", poisson(), "sim Poisson")
a2 <- aghq_block("yb", "binomial", binomial(), "sim logit")

# quadrature vs Laplace, and the quadrature integral vs a brute-force grid
fp1 <- panglm(yc ~ time + x, data = d, index = idx, model = "random", family = "poisson",
              random = ~ time, R = 1)
fp7 <- panglm(yc ~ time + x, data = d, index = idx, model = "random", family = "poisson",
              random = ~ time)
fp15 <- panglm(yc ~ time + x, data = d, index = idx, model = "random", family = "poisson",
               random = ~ time, R = 15)
fb7 <- panglm(yb ~ time + x, data = d, index = idx, model = "random", family = "binomial",
              random = ~ time)
fb15 <- panglm(yb ~ time + x, data = d, index = idx, model = "random", family = "binomial",
               random = ~ time, R = 15)

# Brute-force integral: trapezoid rule on a 300 x 300 grid over b, centred
# at each individual's mode, +- 6 prior standard deviations wide (the
# trapezoid rule converges geometrically for smooth peaked integrands)
brute_loglik <- function(fit, y, density) {
  S <- fit$varcomp; Si <- solve(S)
  sds <- sqrt(diag(S))
  eta <- as.numeric(model.matrix(~ time + x, d) %*% coef(fit))
  rows_by_id <- split(seq_len(nrow(d)), d$id)
  sum(vapply(seq_along(rows_by_id), function(i) {
    rows <- rows_by_id[[i]]
    g1 <- seq(-6, 6, length.out = 300) * sds[1] + fit$ranef[i, 1]
    g2 <- seq(-6, 6, length.out = 300) * sds[2] + fit$ranef[i, 2]
    bb1 <- rep(g1, 300); bb2 <- rep(g2, each = 300)
    lp <- outer(eta[rows], bb1, "+") + outer(d$time[rows], bb2)
    ll <- colSums(density(y[rows], lp)) - log(2 * pi) - 0.5 * log(det(S)) -
      0.5 * (Si[1, 1] * bb1^2 + 2 * Si[1, 2] * bb1 * bb2 + Si[2, 2] * bb2^2)
    top <- max(ll)
    top + log(sum(exp(ll - top)) * diff(g1[1:2]) * diff(g2[1:2]))
  }, numeric(1)))
}
dpois_log <- function(y, lp) y * lp - exp(lp) - lgamma(y + 1)
dlogit_log <- function(y, lp) plogis((2 * y - 1) * lp, log.p = TRUE)

q1 <- rbind(
  vcompare("slopes quadrature", "sim Poisson", "loglik R = 15 vs dense grid",
           c(loglik = fp15$loglik), list(hand = brute_loglik(fp15, d$yc, dpois_log)), 1e-6),
  vcompare("slopes quadrature", "sim logit", "loglik R = 15 vs dense grid",
           c(loglik = fb15$loglik), list(hand = brute_loglik(fb15, d$yb, dlogit_log)), 1e-6),
  vcompare("slopes quadrature", "sim Poisson", "coef R = 7 vs R = 15", coef(fp7),
           list(`panglm R = 15` = coef(fp15)), 1e-4),
  vcompare("slopes quadrature", "sim logit", "coef R = 7 vs R = 15", coef(fb7),
           list(`panglm R = 15` = coef(fb15)), 1e-3),
  vcompare("slopes quadrature", "sim Poisson", "coef Laplace vs R = 15", coef(fp1),
           list(`panglm R = 15` = coef(fp15)), 1e-2),
  vcompare("slopes quadrature", "sim Poisson", "varcomp Laplace vs R = 15", vech(fp1$varcomp),
           list(`panglm R = 15` = vech(fp15$varcomp)), 1e-2)
)

# population-averaged partial effects and predictions by hand
ape_p <- partialeffects(fp7, "x")
Zp <- cbind(1, d$time)
s2 <- rowSums((Zp %*% fp7$varcomp) * Zp)
eta_p <- predict(fp7, type = "link")
ape_b <- partialeffects(fb7, "x")
gh40 <- panglm:::gauss_hermite_quadrature(40)
eta_b <- predict(fb7, type = "link")
s2b <- rowSums((Zp %*% fb7$varcomp) * Zp)
dens_b <- vapply(seq_along(eta_b), function(i) {
  sum(gh40$weights * dlogis(eta_b[i] + sqrt(2 * s2b[i]) * gh40$nodes))
}, numeric(1))
p1 <- rbind(
  vcompare("slopes APE", "sim Poisson", "APE x (population averaged)",
           c(x = ape_p$estimate), list(hand = mean(coef(fp7)[["x"]] * exp(eta_p + s2 / 2))), 1e-6),
  vcompare("slopes APE", "sim logit", "APE x (population averaged)",
           c(x = ape_b$estimate), list(hand = mean(coef(fb7)[["x"]] * dens_b)), 1e-6),
  vcompare("slopes APE", "sim Poisson", "predict(response)",
           predict(fp7, newdata = d[1:12, ], type = "response"),
           list(hand = exp(eta_p + s2 / 2)[1:12]), 1e-8)
)

# ---- Poisson, logit, probit: two random slopes (three-dimensional effect)
set.seed(2)
G3 <- 300; T3 <- 10
d3 <- data.frame(id = rep(seq_len(G3), each = T3), time = rep(0:(T3 - 1), G3))
d3$x <- rnorm(G3 * T3); d3$w <- rnorm(G3 * T3)
L3 <- matrix(c(0.8, 0, 0, -0.1, 0.25, 0, 0.2, 0.1, 0.5), 3, 3, byrow = TRUE)
b3 <- matrix(rnorm(G3 * 3), G3) %*% t(L3)
lin3 <- 0.1 * d3$time + 0.4 * d3$x - 0.3 * d3$w + b3[d3$id, 1] + b3[d3$id, 2] * d3$time +
  b3[d3$id, 3] * d3$w
d3$yc <- rpois(G3 * T3, exp(0.2 + lin3))
d3$yb <- rbinom(G3 * T3, 1, plogis(-0.2 + lin3))
d3$yp <- rbinom(G3 * T3, 1, pnorm(-0.2 + 0.7 * lin3))

slopes2_block <- function(y, family, glm_family, label) {
  fm <- stats::as.formula(paste(y, "~ time + x + w"))
  f1 <- panglm(fm, data = d3, index = idx, model = "random", family = family,
               random = ~ time + w, R = 1)
  tm <- glmmTMB(stats::update(fm, . ~ . + (1 + time + w | id)), data = d3, family = glm_family)
  f7 <- panglm(fm, data = d3, index = idx, model = "random", family = family,
               random = ~ time + w, R = 7)
  # GLMMadaptive started at the glmmTMB estimates (from its own default
  # start it can stop at a far worse optimum in three dimensions)
  ga <- mixed_model(fm, random = ~ time + w | id, data = d3, family = glm_family, nAGQ = 7,
                    initial_values = list(betas = fixef(tm)$cond,
                                          D = unname(VarCorr(tm)$cond$id[, ])),
                    control = list(iter_EM = 0, update_GH_every = 1,
                                   tol1 = 1e-10, tol2 = 1e-10, tol3 = 1e-12))
  # with three dimensions the two packages scale the nodes with different
  # square roots of the curvature: at R = 7 the rules agree to about 1e-3
  rbind(
    vcompare("slopes Laplace, 2 slopes", label, "coef", coef(f1),
             list(glmmTMB = fixef(tm)$cond), 1e-4),
    vcompare("slopes Laplace, 2 slopes", label, "se", sqrt(diag(vcov(f1))),
             list(glmmTMB = sqrt(diag(vcov(tm)$cond))), 1e-4),
    vcompare("slopes Laplace, 2 slopes", label, "varcomp (vech Sigma)", vech(f1$varcomp),
             list(glmmTMB = tmb_vech(tm)), 1e-4),
    vcompare("slopes Laplace, 2 slopes", label, "loglik", c(loglik = f1$loglik),
             list(glmmTMB = as.numeric(logLik(tm))), 1e-6),
    vcompare("slopes Laplace, 2 slopes", label, "ranef (modes)", c(f1$ranef),
             list(glmmTMB = c(as.matrix(ranef(tm)$cond$id))), 1e-3),
    vcompare("slopes AGHQ R = 7, 2 slopes", label, "coef", coef(f7),
             list(GLMMadaptive = fixef(ga)), 1e-3),
    vcompare("slopes AGHQ R = 7, 2 slopes", label, "se", sqrt(diag(vcov(f7))),
             list(GLMMadaptive = sqrt(diag(vcov(ga)))[1:4]), 1e-3),
    vcompare("slopes AGHQ R = 7, 2 slopes", label, "varcomp (vech Sigma)", vech(f7$varcomp),
             list(GLMMadaptive = vech(ga$D)), 1e-3),
    vcompare("slopes AGHQ R = 7, 2 slopes", label, "loglik", c(loglik = f7$loglik),
             list(GLMMadaptive = as.numeric(logLik(ga))), 1e-5)
  )
}
m1 <- slopes2_block("yc", "poisson", poisson(), "sim Poisson")
m2 <- slopes2_block("yb", "binomial", binomial(), "sim logit")
m3 <- slopes2_block("yp", "probit", binomial("probit"), "sim probit")

# ---- boundary: GLM slope variance at zero (singular covariance)
# Data with a random intercept only: the fitted intercept-slope correlation
# is 1 (a singular Sigma, lme4 reports a boundary fit). The fixed-effect
# standard errors stay finite and match glmmTMB; panglm warns.
set.seed(4)
G4 <- 200; T4 <- 6
d4 <- data.frame(id = rep(seq_len(G4), each = T4), time = rep(0:(T4 - 1), G4))
d4$x <- rnorm(G4 * T4)
u4 <- rnorm(G4, 0, 0.6)
d4$yc <- rpois(G4 * T4, exp(0.3 + 0.1 * d4$time + 0.4 * d4$x + u4[d4$id]))
w4 <- NULL
f4 <- withCallingHandlers(
  panglm(yc ~ time + x, data = d4, index = idx, model = "random", family = "poisson",
         random = ~ time, R = 1),
  warning = function(w) { w4 <<- conditionMessage(w); invokeRestart("muffleWarning") })
stopifnot(grepl("singular", w4), all(is.finite(sqrt(diag(vcov(f4))))))
tm4 <- suppressWarnings(glmmTMB(yc ~ time + x + (1 + time | id), data = d4, family = poisson()))
gm4 <- suppressMessages(glmer(yc ~ time + x + (1 + time | id), data = d4, family = poisson()))
stopifnot(isSingular(gm4))
b1 <- rbind(
  vcompare("slopes boundary", "sim Poisson", "coef", coef(f4),
           list(glmmTMB = fixef(tm4)$cond), 1e-4),
  vcompare("slopes boundary", "sim Poisson", "se", sqrt(diag(vcov(f4))),
           list(glmmTMB = sqrt(diag(vcov(tm4)$cond))), 1e-4),
  vcompare("slopes boundary", "sim Poisson", "loglik", c(loglik = f4$loglik),
           list(glmmTMB = as.numeric(logLik(tm4)), lme4 = as.numeric(logLik(gm4))), 1e-5)
)

# ---- negative binomial (NB2) and gamma (log link) with a random slope
# Lognormal (normal on the log scale) random effects, as glmmTMB's nbinom2
# and Gamma(link = "log"): NB2 theta is glmmTMB's sigma, and the gamma
# shape is 1 / sigma^2 (glmmTMB's sigma is the coefficient of variation).
set.seed(7)
G5 <- 200; T5 <- 6
d5 <- data.frame(id = rep(seq_len(G5), each = T5), time = rep(0:(T5 - 1), G5))
d5$x <- rnorm(G5 * T5)
b5 <- matrix(rnorm(G5 * 2), G5) %*% t(matrix(c(0.6, 0, -0.1, 0.2), 2, 2, byrow = TRUE))
mu5 <- exp(0.5 + 0.1 * d5$time + 0.4 * d5$x + b5[d5$id, 1] + b5[d5$id, 2] * d5$time)
d5$yn <- rnbinom(G5 * T5, mu = mu5, size = 2)
d5$yg <- rgamma(G5 * T5, shape = 3, rate = 3 / mu5)

disp_block <- function(y, family, tmb_family, label, disp, tmb_disp) {
  fm <- stats::as.formula(paste(y, "~ time + x"))
  f <- panglm(fm, data = d5, index = idx, model = "random", family = family,
              random = ~ time, R = 1)
  tm <- glmmTMB(stats::update(fm, . ~ . + (1 + time | id)), data = d5, family = tmb_family)
  rbind(
    vcompare("slopes Laplace", label, "coef", coef(f), list(glmmTMB = fixef(tm)$cond), 1e-4),
    vcompare("slopes Laplace", label, "se", sqrt(diag(vcov(f))),
             list(glmmTMB = sqrt(diag(vcov(tm)$cond))), 1e-4),
    vcompare("slopes Laplace", label, "varcomp (vech Sigma)", vech(f$varcomp),
             list(glmmTMB = tmb_vech(tm)), 1e-4),
    vcompare("slopes Laplace", label, "dispersion parameter", c(disp = disp(f)),
             list(glmmTMB = tmb_disp(sigma(tm))), 1e-4),
    vcompare("slopes Laplace", label, "loglik", c(loglik = f$loglik),
             list(glmmTMB = as.numeric(logLik(tm))), 1e-6),
    vcompare("slopes Laplace", label, "ranef (modes)", c(f$ranef),
             list(glmmTMB = c(as.matrix(ranef(tm)$cond$id))), 1e-3)
  )
}
n1 <- disp_block("yn", "negbin", nbinom2(), "sim NB2", function(f) f$negbin_theta,
                 function(s) s)
n2 <- disp_block("yg", "gamma", Gamma(link = "log"), "sim gamma", function(f) f$shape,
                 function(s) 1 / s^2)

# NB2 adaptive quadrature (R = 7) vs GLMMadaptive (its negative.binomial()
# is NB2 with phi = theta); GLMMadaptive's Gamma.fam() fails to start here
tm5 <- glmmTMB(yn ~ time + x + (1 + time | id), data = d5, family = nbinom2())
fn7 <- panglm(yn ~ time + x, data = d5, index = idx, model = "random", family = "negbin",
              random = ~ time, R = 7)
gan <- mixed_model(yn ~ time + x, random = ~ time | id, data = d5,
                   family = negative.binomial(), nAGQ = 7,
                   initial_values = list(betas = fixef(tm5)$cond,
                                         D = unname(VarCorr(tm5)$cond$id[, ])),
                   control = list(iter_EM = 0, update_GH_every = 1,
                                  tol1 = 1e-10, tol2 = 1e-10, tol3 = 1e-12))
n3 <- rbind(
  vcompare("slopes AGHQ R = 7", "sim NB2", "coef", coef(fn7), list(GLMMadaptive = fixef(gan)), 1e-4),
  vcompare("slopes AGHQ R = 7", "sim NB2", "se", sqrt(diag(vcov(fn7))),
           list(GLMMadaptive = sqrt(diag(vcov(gan)))[1:3]), 1e-4),
  vcompare("slopes AGHQ R = 7", "sim NB2", "varcomp (vech Sigma)", vech(fn7$varcomp),
           list(GLMMadaptive = vech(gan$D)), 1e-4),
  vcompare("slopes AGHQ R = 7", "sim NB2", "dispersion parameter", c(theta = fn7$negbin_theta),
           list(GLMMadaptive = exp(gan$phis)), 1e-4),
  vcompare("slopes AGHQ R = 7", "sim NB2", "loglik", c(loglik = fn7$loglik),
           list(GLMMadaptive = as.numeric(logLik(gan))), 1e-6)
)

vsave(g1, g2, g3, g4, g5, l1, l2, l3, a1, a2, q1, p1, m1, m2, m3, b1, n1, n2, n3,
      name = "slopes")
