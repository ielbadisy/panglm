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

test_that("panglm_dispersiontest() is a deprecated alias of pandispersiontest()", {
  data(copd, package = "panglm")
  fit <- panglm(exacerbations ~ treatment + crp, data = copd, index = c("id", "visit"),
                model = "pooling", family = "poisson")
  panglm:::deprecate_reset()
  expect_warning(old <- panglm_dispersiontest(fit), "pandispersiontest\\(\\)",
                 class = "deprecatedWarning")
  expect_identical(old, pandispersiontest(fit))
  expect_s3_class(old, "panglm_dispersiontest")
  expect_no_warning(panglm_dispersiontest(fit))
})

test_that("panglm_hausman() is a deprecated alias of panhausman()", {
  data(copd, package = "panglm")
  fe <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"), model = "within")
  re <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"), model = "random")
  panglm:::deprecate_reset()
  expect_warning(old <- panglm_hausman(fe, re), "panhausman\\(\\)",
                 class = "deprecatedWarning")
  expect_identical(old, panhausman(fe, re))
  expect_s3_class(old, "panglm_hausman")
  expect_identical(
    suppressWarnings(panglm_hausman(fe, method = "mundlak", data = copd, vcov = "HC1")),
    panhausman(fe, method = "mundlak", data = copd, vcov = "HC1")
  )
  expect_output(print(old), "Hausman")
})

test_that("panglm_hurdle() is a deprecated alias of panhurdle(), count_family = family", {
  data(copd, package = "panglm")
  panglm:::deprecate_reset()
  expect_warning(
    old <- panglm_hurdle(exacerbations ~ crp, data = copd, index = c("id", "visit")),
    "panhurdle\\(\\)", class = "deprecatedWarning"
  )
  expect_identical(old, panhurdle(exacerbations ~ crp, data = copd, index = c("id", "visit")))
  old_nb <- panglm_hurdle(exacerbations ~ crp, data = copd, index = c("id", "visit"),
                          count_family = "negbin")
  new_nb <- panhurdle(exacerbations ~ crp, data = copd, index = c("id", "visit"),
                      family = "negbin")
  expect_identical(old_nb, new_nb)
  expect_identical(old_nb$count_family, "negbin")
  expect_s3_class(old_nb, "panglm_hurdle")
  expect_output(print(old_nb), "NB2")
})
