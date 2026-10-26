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

vsave(gaussian_re(inv ~ value + capital, Grunfeld, c("firm", "year"), "Grunfeld"),
      gaussian_re(inv ~ value + capital, grun_u, c("firm", "year"), "Grunfeld unbalanced"),
      gaussian_re(hie_fm, hie, c("id", "year"), "HealthIns (unbalanced)"),
      hausman_u,
      name = "random")
