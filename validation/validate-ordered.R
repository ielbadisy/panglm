# Validation of the ordered outcome models (family = "ologit" / "oprobit").
#
# - Pooled ordered logit and probit: MASS::polr and ordinal::clm maximize
#   the same likelihood (P(y <= j) = F(c_j - x'beta)). clm has an analytic
#   Hessian, so its standard errors are compared tightly; polr's Hessian
#   comes from optim, so its standard errors are compared at 1e-5.
# - Fixed-effects ordered logit (Baetschmann, Staub and Winkelmann 2015,
#   "blow-up and cluster"): no R package implements the BUC estimator
#   (none on CRAN as of 2026-10), so the reference is the definition
#   itself: the data are blown up by hand (one copy per cutoff, outcome
#   1(y >= k), one stratum per individual and cutoff) and fitted with
#   survival::clogit (exact conditional likelihood). The cluster sandwich
#   is rebuilt from clogit's information and per-individual scores, with
#   panglm's small-sample factor G / (G - 1).
# - Random-intercept ordered logit and probit: ordinal::clmm with adaptive
#   Gauss-Hermite quadrature (nAGQ = 41) against panglm's adaptive
#   Gauss-Hermite quadrature with R = 41 nodes (validate-random.R also
#   compares the default R = 21 and the Laplace fit R = 1 on
#   lme4::schizophrenia). Both approximate the same integrated likelihood to far
#   below 1e-8 at this size, so the remaining difference is clmm's
#   optimizer tolerance (its fit is pushed with rel.tol = 1e-14); compared
#   at 1e-5. clmm's standard errors come from a finite-difference Hessian.
source("validation/common.R")
suppressPackageStartupMessages({
  library(MASS)
  library(ordinal)
  library(survival)
})

set.seed(3)
G <- 200; Tt <- 5
d <- data.frame(id = rep(seq_len(G), each = Tt), t = rep(seq_len(Tt), G))
d$x1 <- rnorm(nrow(d))
d$x2 <- rnorm(nrow(d)) + rep(rnorm(G, sd = 0.5), each = Tt)
d$z <- rep(rbinom(G, 1, 0.5), each = Tt)
u <- rep(rnorm(G), each = Tt)
ys <- 0.8 * d$x1 - 0.4 * d$x2 - 0.5 * d$z + u + rlogis(nrow(d))
d$y <- cut(ys, c(-Inf, -1, 0.5, 2, Inf), labels = c("poor", "fair", "good", "vgood"),
           ordered_result = TRUE)
d$idf <- factor(d$id)
data_name <- "simulated ordered panel (G = 200, T = 5, J = 4)"
f <- y ~ x1 + x2 + z
index <- c("id", "t")
res <- list()

## Pooled ordered logit / probit -------------------------------------------
for (link in c("logit", "probit")) {
  fam <- if (link == "logit") "ologit" else "oprobit"
  label <- paste("pooled", fam)
  p <- panglm(f, d, index = index, model = "pooling", family = fam)
  m <- polr(f, d, Hess = TRUE, method = if (link == "logit") "logistic" else "probit",
            control = list(reltol = 1e-15, maxit = 5000))
  c1 <- clm(f, data = d, link = link, control = clm.control(gradTol = 1e-12))
  cl_order <- c(names(c1$beta), names(c1$alpha))
  res[[length(res) + 1]] <- vcompare(label, data_name, "coefficients", coef(p),
    list(`MASS::polr` = c(coef(m), m$zeta), `ordinal::clm` = coef(c1)[cl_order]), tol = 1e-6)
  res[[length(res) + 1]] <- vcompare(label, data_name, "std.errors", sqrt(diag(vcov(p))),
    list(`ordinal::clm` = sqrt(diag(vcov(c1)))[cl_order]), tol = 1e-6)
  # polr's Hessian is a finite-difference optim Hessian
  res[[length(res) + 1]] <- vcompare(label, data_name, "std.errors (numerical Hessian)",
    sqrt(diag(vcov(p))), list(`MASS::polr` = sqrt(diag(vcov(m)))), tol = 1e-5)
  res[[length(res) + 1]] <- vcompare(label, data_name, "logLik", as.numeric(logLik(p)),
    list(`MASS::polr` = as.numeric(logLik(m)), `ordinal::clm` = as.numeric(logLik(c1))),
    tol = 1e-8)
}

## Fixed-effects ordered logit: blow-up and cluster ------------------------
# z is time invariant and not identified by the FE estimator
f_fe <- y ~ x1 + x2
pf <- panglm(f_fe, d, index = index, model = "within", family = "ologit")
yc <- as.integer(d$y)
J <- nlevels(d$y)
blown <- do.call(rbind, lapply(2:J, function(k) {
  b <- d
  b$dy <- as.numeric(yc >= k)
  b$stratum <- paste(d$id, k)
  b
}))
ctrl <- coxph.control(eps = 1e-11, iter.max = 100)
cf <- clogit(dy ~ x1 + x2 + strata(stratum), data = blown, method = "exact", control = ctrl)
# Cluster sandwich: survival does not give dfbeta residuals for the exact
# conditional likelihood, so the individual scores are central differences
# of each individual's exact clogit log-likelihood (iter.max = 0 evaluates
# it at the given coefficients), and the bread is clogit's inverse
# information.
informative <- tapply(blown$dy, blown$stratum, function(v) length(unique(v)) == 2)
blown_inf <- blown[blown$stratum %in% names(informative)[informative], ]
loglik_at <- function(b, data) {
  clogit(dy ~ x1 + x2 + strata(stratum), data = data, method = "exact", init = b,
         control = coxph.control(iter.max = 0))$loglik[1]
}
scores <- t(vapply(split(blown_inf, blown_inf$id), function(di) {
  vapply(seq_along(coef(cf)), function(j) {
    h <- 1e-5
    up <- dn <- coef(cf); up[j] <- up[j] + h; dn[j] <- dn[j] - h
    (loglik_at(up, di) - loglik_at(dn, di)) / (2 * h)
  }, numeric(1))
}, numeric(length(coef(cf)))))
Gc <- nrow(scores)
V_cl <- vcov(cf) %*% crossprod(scores) %*% vcov(cf) * Gc / (Gc - 1)
res[[length(res) + 1]] <- vcompare("FE ordered logit (BUC)", data_name, "coefficients", coef(pf),
  list(`survival::clogit (hand blow-up)` = coef(cf)), tol = 1e-6)
res[[length(res) + 1]] <- vcompare("FE ordered logit (BUC)", data_name, "classical std.errors",
  sqrt(diag(vcov(pf, type = "classical"))),
  list(`survival::clogit (hand blow-up)` = sqrt(diag(vcov(cf)))), tol = 1e-6)
res[[length(res) + 1]] <- vcompare("FE ordered logit (BUC)", data_name, "cluster std.errors",
  sqrt(diag(vcov(pf))), list(`survival::clogit (hand blow-up)` = sqrt(diag(V_cl))), tol = 1e-6)
res[[length(res) + 1]] <- vcompare("FE ordered logit (BUC)", data_name, "logLik",
  as.numeric(logLik(pf)), list(`survival::clogit (hand blow-up)` = cf$loglik[2]), tol = 1e-8)

## Random-intercept ordered logit / probit ---------------------------------
for (link in c("logit", "probit")) {
  fam <- if (link == "logit") "ologit" else "oprobit"
  label <- paste("RE", fam)
  pr <- panglm(f, d, index = index, model = "random", family = fam, R = 41)
  cm <- clmm(y ~ x1 + x2 + z + (1 | idf), data = d, link = link, nAGQ = 41,
             control = clmm.control(method = "nlminb", maxIter = 1000, gradTol = 1e-10,
                                    rel.tol = 1e-14, eval.max = 5000, iter.max = 5000))
  cm_order <- c(names(cm$beta), names(cm$alpha))
  res[[length(res) + 1]] <- vcompare(label, data_name, "coefficients", coef(pr),
    list(`ordinal::clmm` = c(cm$beta, cm$alpha)), tol = 1e-5)
  res[[length(res) + 1]] <- vcompare(label, data_name, "sigma_u", pr$sigma,
    list(`ordinal::clmm` = as.numeric(cm$ST$idf)), tol = 1e-5)
  res[[length(res) + 1]] <- vcompare(label, data_name, "std.errors", sqrt(diag(vcov(pr))),
    list(`ordinal::clmm` = sqrt(diag(vcov(cm)))[cm_order]), tol = 1e-5)
  res[[length(res) + 1]] <- vcompare(label, data_name, "logLik", as.numeric(logLik(pr)),
    list(`ordinal::clmm` = as.numeric(logLik(cm))), tol = 1e-8)
}

do.call(vsave, c(res, name = "ordered"))
