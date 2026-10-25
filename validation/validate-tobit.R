# Tobit: panglm(family = tobitfamily()) vs AER::tobit and censReg (pooled),
# and censReg on a pdata.frame (random intercept, Gauss-Hermite quadrature).
# One-limit (left = 0) and two-limit (0, 2) censoring.
source("validation/common.R")
suppressPackageStartupMessages({ library(AER); library(censReg); library(plm); library(sandwich) })

set.seed(11)
G <- 150; T <- 5; n <- G * T
d <- data.frame(id = rep(seq_len(G), each = T), t = rep(seq_len(T), G),
                x1 = rnorm(n), x2 = rbinom(n, 1, 0.4))
ystar <- 0.3 + 0.8 * d$x1 - 0.5 * d$x2 + rep(rnorm(G, 0, 0.8), each = T) + rnorm(n)
d$y1 <- pmax(0, ystar)
d$y2 <- pmin(2, pmax(0, ystar))
pd <- pdata.frame(d, index = c("id", "t"))
idx <- c("id", "t")

pooled <- function(y, right, label) {
  fm <- stats::as.formula(paste(y, "~ x1 + x2"))
  f <- panglm(fm, data = d, index = idx, model = "pooling", family = tobitfamily(0, right))
  a <- tobit(fm, left = 0, right = right, data = d)
  cr <- censReg(fm, left = 0, right = right, data = d)
  k <- 1:3
  # Phi-weighted APE of a continuous regressor, by hand from the AER fit
  xb <- predict(a, type = "lp"); s <- a$scale
  p_in <- pnorm((right - xb) / s) - pnorm((0 - xb) / s)
  ape <- partialeffects(f)
  rbind(
    vcompare("tobit pooled", label, "coef", coef(f), list(AER = coef(a), censReg = coef(cr)[k]), 1e-6),
    vcompare("tobit pooled", label, "log sigma", c(logsigma = log(f$sigma_e)),
             list(AER = log(a$scale), censReg = coef(cr)[4]), 1e-6),
    vcompare("tobit pooled", label, "se", sqrt(diag(vcov(f))),
             list(AER = sqrt(diag(vcov(a)))[k], censReg = sqrt(diag(vcov(cr)))[k]), 1e-5),
    vcompare("tobit pooled", label, "se_HC1", sqrt(diag(vcov(f, type = "HC1"))),
             list(sandwich = sqrt(diag(sandwich(a) * n / (n - 4)))[k]), 1e-5),
    vcompare("tobit pooled", label, "loglik", c(loglik = f$loglik),
             list(AER = as.numeric(logLik(a)), censReg = as.numeric(logLik(cr))), 1e-6),
    vcompare("tobit pooled", label, "APE x1", c(x1 = ape$estimate[ape$term == "x1"]),
             list(hand = mean(p_in) * coef(a)[["x1"]]), 1e-6)
  )
}

random <- function(y, right, label) {
  fm <- stats::as.formula(paste(y, "~ x1 + x2"))
  f <- panglm(fm, data = d, index = idx, model = "random", family = tobitfamily(0, right))
  cr <- censReg(fm, left = 0, right = right, data = pd, nGHQ = 30)
  k <- 1:3
  rbind(
    vcompare("tobit random", label, "coef", coef(f), list(censReg = coef(cr)[k]), 1e-4),
    vcompare("tobit random", label, "log sigma_u, log sigma_e",
             c(logsigma_u = log(f$sigma), logsigma_e = log(f$sigma_e)),
             list(censReg = coef(cr)[4:5]), 1e-4),
    vcompare("tobit random", label, "se", sqrt(diag(vcov(f))),
             list(censReg = sqrt(diag(vcov(cr)))[k]), 1e-4),
    vcompare("tobit random", label, "loglik", c(loglik = f$loglik),
             list(censReg = as.numeric(logLik(cr))), 1e-6)
  )
}

vsave(pooled("y1", Inf, "sim left = 0"), pooled("y2", 2, "sim [0, 2]"),
      random("y1", Inf, "sim left = 0"), random("y2", 2, "sim [0, 2]"),
      name = "tobit")
