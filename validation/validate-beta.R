# Beta regression: panglm(family = "beta") against
# * betareg::betareg() for the pooled model (logit and probit links, with
#   and without case weights): coefficients, log precision, standard
#   errors and log-likelihood. panglm's standard errors use the observed
#   information, so betareg is run with hessian = TRUE (its default
#   covariance is the expected information, about 0.5% different here);
# * GLMMadaptive::mixed_model(family = beta.fam(), nAGQ = 21) for the random
#   intercept: the same adaptive Gauss-Hermite rule, so the maximum
#   likelihood estimates and the log-likelihood agree to optimizer
#   precision (GLMMadaptive run with tight tolerances);
# * the log-likelihood at the panglm estimate by one-dimensional numerical
#   integration (stats::integrate) of every individual's contribution;
# * glmmTMB (Laplace approximation) for the random intercept, reported with
#   a loose tolerance: the Laplace approximation differs from 21-node
#   quadrature.
# Data: simulated panel, 400 individuals over 5 periods, logit mean with a
# Gaussian random intercept (sd 0.5) and precision 15.
source("validation/common.R")
suppressPackageStartupMessages({
  library(betareg)
  library(GLMMadaptive)
  library(glmmTMB)
})
options(panglm.threads = 6)

set.seed(7)
G <- 400; TT <- 5
d <- data.frame(id = rep(seq_len(G), each = TT), t = rep(seq_len(TT), G),
                x = rnorm(G * TT), z = rbinom(G * TT, 1, 0.4))
mu <- plogis(0.2 + 0.6 * d$x - 0.4 * d$z + rep(rnorm(G, 0, 0.5), each = TT))
d$y <- rbeta(G * TT, mu * 15, (1 - mu) * 15)
d$w <- rep(sample(1:3, G, TRUE), each = TT)
idx <- c("id", "t")
sname <- "simulated N=400 T=5"
se <- function(v) sqrt(diag(v))
out <- list()

# ---- pooled -----------------------------------------------------------------
pooled <- function(link, weighted) {
  fam <- betafamily(link)
  if (weighted) {
    f <- panglm(y ~ x + z, d, idx, model = "pooling", family = fam, weights = d$w)
    b <- betareg(y ~ x + z, d, link = link, weights = w, hessian = TRUE)
  } else {
    f <- panglm(y ~ x + z, d, idx, model = "pooling", family = fam)
    b <- betareg(y ~ x + z, d, link = link, hessian = TRUE)
  }
  lab <- paste0("beta pooled (", link, if (weighted) ", weighted", ")")
  k <- 1:3
  rbind(
    vcompare(lab, sname, "coef", coef(f), list(betareg = coef(b)[k]), 1e-6),
    vcompare(lab, sname, "log precision", c(logphi = log(f$precision)),
             list(betareg = log(coef(b)[4])), 1e-6),
    vcompare(lab, sname, "se", se(vcov(f)), list(betareg = se(vcov(b))[k]), 1e-5),
    vcompare(lab, sname, "logLik", c(loglik = f$loglik),
             list(betareg = as.numeric(logLik(b))), 1e-8)
  )
}
out$pooled <- rbind(pooled("logit", FALSE), pooled("probit", FALSE), pooled("logit", TRUE))

# ---- random intercept -------------------------------------------------------
m <- panglm(y ~ x + z, d, idx, model = "random", family = "beta")
ga <- mixed_model(y ~ x + z, random = ~ 1 | id, data = d, family = beta.fam(), nAGQ = 21,
                  control = list(iter_EM = 0, iter_qN_outer = 50, tol1 = 1e-12,
                                 tol2 = 1e-12, tol3 = 1e-12, update_GH_every = 1))
lab <- "beta random intercept, R = 21"
out$random <- rbind(
  vcompare(lab, sname, "coef", coef(m), list(GLMMadaptive = fixef(ga)), 1e-6),
  vcompare(lab, sname, "log precision", c(logphi = log(m$precision)),
           list(GLMMadaptive = ga$phis), 1e-6),
  vcompare(lab, sname, "sigma", c(sigma = m$sigma), list(GLMMadaptive = sqrt(ga$D[1])), 1e-6),
  vcompare(lab, sname, "logLik", c(loglik = m$loglik),
           list(GLMMadaptive = as.numeric(logLik(ga))), 1e-6),
  vcompare(lab, sname, "se", se(vcov(m)), list(GLMMadaptive = se(vcov(ga))[1:3]), 1e-3)
)

# log-likelihood at the panglm estimate by numerical integration
th <- m$theta
eta <- as.numeric(m$X %*% th[1:3]); phi <- exp(th[4]); s <- exp(th[5])
grp <- rep(seq_along(m$group_size), m$group_size)
ll_int <- sum(vapply(seq_along(m$group_size), function(i) {
  e <- eta[grp == i]; yy <- m$y[grp == i]
  f <- function(u) vapply(u, function(v) {
    p <- plogis(e + v)
    exp(sum(stats::dbeta(yy, p * phi, (1 - p) * phi, log = TRUE)))
  }, numeric(1)) * stats::dnorm(u, 0, s)
  log(stats::integrate(f, -Inf, Inf, rel.tol = 1e-12)$value)
}, numeric(1)))
out$integrate <- vcompare(lab, sname, "logLik at estimate", c(loglik = m$loglik),
                          list(integrate = ll_int), 1e-8)

# Laplace (glmmTMB): approximate agreement only
tm <- glmmTMB(y ~ x + z + (1 | id), d, family = beta_family())
out$laplace <- rbind(
  vcompare(lab, sname, "coef (vs Laplace)", coef(m), list(glmmTMB = fixef(tm)$cond), 1e-2),
  vcompare(lab, sname, "sigma (vs Laplace)", c(sigma = m$sigma),
           list(glmmTMB = sqrt(VarCorr(tm)$cond$id[1])), 1e-2)
)

do.call(vsave, c(unname(out), name = "beta"))
