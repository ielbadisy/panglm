test_that("deprecated family constructors warn once and match their replacements", {
  pairs <- list(
    list(gaussian_family, gaussianfamily, list()),
    list(poisson_family, poissonfamily, list()),
    list(binomial_family, binomialfamily, list()),
    list(binomial_family, binomialfamily, list("probit")),
    list(negbin_family, negbinfamily, list()),
    list(gamma_family, gammafamily, list()),
    list(beta_family, betafamily, list("probit")),
    list(ordered_family, orderedfamily, list("logit")),
    list(tobit_family, tobitfamily, list(0, 2))
  )
  for (p in pairs) {
    panglm:::deprecate_reset()
    expect_warning(old <- do.call(p[[1]], p[[3]]), class = "deprecatedWarning")
    expect_identical(old, do.call(p[[2]], p[[3]]))
    # once per session: a second call is silent
    expect_no_warning(do.call(p[[1]], p[[3]]))
  }
  panglm:::deprecate_reset()
  expect_warning(binomial_family(), "binomialfamily\\(\\)")
})

test_that("a fit with a deprecated family spec is identical to the new one", {
  data(copd, package = "panglm")
  panglm:::deprecate_reset()
  f_new <- panglm(hospitalized ~ crp, data = copd, index = c("id", "visit"),
                  model = "within", family = binomialfamily("logit"))
  expect_warning(
    f_old <- panglm(hospitalized ~ crp, data = copd, index = c("id", "visit"),
                    model = "within", family = binomial_family("logit")),
    class = "deprecatedWarning"
  )
  expect_identical(coef(f_old), coef(f_new))
  expect_identical(vcov(f_old), vcov(f_new))
})
