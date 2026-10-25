# Hausman-Taylor (1981): paniv(model = "hausmantaylor") against plm::pht()
# and plm::plm(model = "random", random.method = "ht",
# inst.method = "baltagi"), the Cornwell-Rupert wage equation of Baltagi
# (2013, Table 7.4) on plm::Wages (595 individuals, 7 years). An unbalanced
# version drops 600 random person-years (compared with pht() only: plm's
# random.method = "ht" uses the maximal T in theta, pht() and paniv() use
# the individual T_i and the harmonic mean of T_i).
source("validation/common.R")
suppressPackageStartupMessages(library(plm))
data("Wages", package = "plm")
Wages$id <- rep(1:595, each = 7)
Wages$t <- rep(1:7, 595)
set.seed(20261025)
WagesU <- Wages[-sample(nrow(Wages), 600), ]
idx <- c("id", "t")
f <- lwage ~ wks + south + smsa + married + exp + I(exp^2) + bluecol + ind + union +
  sex + black + ed | sex + black + bluecol + south + smsa + ind
se <- function(v) sqrt(diag(v))
out <- list()
add <- function(...) out[[length(out) + 1L]] <<- vcompare(...)

for (dname in c("Wages", "Wages unbalanced")) {
  d <- if (dname == "Wages") Wages else WagesU
  a <- paniv(f, d, index = idx, model = "hausmantaylor")
  b <- suppressWarnings(plm::pht(f, data = d, index = idx))
  refs_coef <- list(plm.pht = coef(b)[names(coef(a))])
  refs_se <- list(plm.pht = se(vcov(b))[names(coef(a))])
  refs_vc <- list(plm.pht = unlist(b$ercomp$sigma2)[c("idios", "id")])
  if (dname == "Wages") {
    b2 <- plm::plm(lwage ~ wks + south + smsa + married + exp + I(exp^2) + bluecol + ind +
                     union + sex + black + ed | bluecol + south + smsa + ind + sex + black |
                     wks + married + exp + I(exp^2) + union,
                   data = d, index = idx, model = "random", random.method = "ht",
                   inst.method = "baltagi")
    refs_coef$plm.ht <- coef(b2)[names(coef(a))]
    refs_se$plm.ht <- se(vcov(b2))[names(coef(a))]
    refs_vc$plm.ht <- b2$ercomp$sigma2[c("idios", "id")]
  }
  add("paniv hausmantaylor", dname, "coefficients", coef(a), refs_coef, 1e-8)
  add("paniv hausmantaylor", dname, "classical SE", se(vcov(a)), refs_se, 1e-8)
  add("paniv hausmantaylor", dname, "variance components", a$ercomp$sigma2[c("idios", "id")],
      refs_vc, 1e-8)
}
do.call(vsave, c(out, name = "hausmantaylor"))
