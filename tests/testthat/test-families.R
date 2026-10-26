test_that("base R family objects map onto panglm families", {
  data(copd)
  ix <- c("id", "visit")
  for (lk in c("logit", "probit")) {
    a <- panglm(hospitalized ~ crp, copd, ix, family = stats::binomial(lk))
    b <- panglm(hospitalized ~ crp, copd, ix, family = binomialfamily(lk))
    expect_equal(coef(a), coef(b))
    expect_identical(a$family$link, lk)
  }
  expect_equal(coef(panglm(exacerbations ~ crp, copd, ix, family = stats::poisson())),
               coef(panglm(exacerbations ~ crp, copd, ix, family = "poisson")))
  expect_equal(coef(panglm(fev1 ~ crp, copd, ix, family = stats::gaussian())),
               coef(panglm(fev1 ~ crp, copd, ix)))
  expect_equal(coef(panglm(fev1 ~ crp, copd, ix, family = stats::Gamma("log"))),
               coef(panglm(fev1 ~ crp, copd, ix, family = "gamma")))
  # the generator itself, not called
  expect_identical(panglm(hospitalized ~ crp, copd, ix, family = stats::binomial)$family$link,
                   "logit")
})

test_that("unsupported base R families stop with a pointer to the panglm constructor", {
  data(copd)
  ix <- c("id", "visit")
  expect_error(panglm(hospitalized ~ crp, copd, ix, family = stats::binomial("cloglog")),
               "binomialfamily")
  expect_error(panglm(fev1 ~ crp, copd, ix, family = stats::Gamma("inverse")), "gammafamily")
  expect_error(panglm(exacerbations ~ crp, copd, ix, family = stats::quasipoisson()),
               "not available in panglm")
  skip_if_missing("MASS")
  expect_error(panglm(exacerbations ~ crp, copd, ix, family = MASS::negative.binomial(2)),
               "negbinfamily")
})
