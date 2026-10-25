# Between estimator: panglm(model = "between") vs plm(model = "between"),
# and lm() on the individual means (classical and HC1 covariances).
source("validation/common.R")
suppressPackageStartupMessages({ library(plm); library(sandwich) })
data("Grunfeld", package = "plm")
data(copd)

check <- function(formula, data, index, label) {
  f <- panglm(formula, data = data, index = index, model = "between")
  p <- plm(formula, data = data, index = index, model = "between")
  means <- stats::aggregate(stats::model.frame(formula, data),
                            list(.id = data[[index[1]]]), mean)
  m <- lm(formula, data = means)
  rbind(
    vcompare("between", label, "coef", coef(f), list(plm = coef(p), lm_means = coef(m)), 1e-8),
    vcompare("between", label, "se", sqrt(diag(vcov(f))),
             list(plm = sqrt(diag(vcov(p))), lm_means = sqrt(diag(vcov(m)))), 1e-8),
    vcompare("between", label, "se_HC1", sqrt(diag(vcov(f, type = "HC1"))),
             list(sandwich = sqrt(diag(vcovHC(m, type = "HC1")))), 1e-8),
    vcompare("between", label, "se_cluster", sqrt(diag(vcov(f, type = "cluster"))),
             list(sandwich = sqrt(diag(vcovHC(m, type = "HC1")))), 1e-8),
    vcompare("between", label, "loglik", c(loglik = as.numeric(logLik(f))),
             list(lm_means = as.numeric(logLik(m))), 1e-8)
  )
}

# balanced (Grunfeld) and unbalanced (copd: drop a random 15% of rows)
set.seed(3)
copd_u <- copd[sort(sample(nrow(copd), round(0.85 * nrow(copd)))), ]
vsave(check(inv ~ value + capital, Grunfeld, c("firm", "year"), "Grunfeld"),
      check(fev1 ~ crp + age, copd_u, c("id", "visit"), "copd unbalanced"),
      name = "between")
