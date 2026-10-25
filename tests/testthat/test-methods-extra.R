test_that("confint.twopart is the Wald interval from coef and vcov", {
  fit <- methods_fit("twopart")
  ci <- confint(fit)
  b <- coef(fit)
  se <- sqrt(diag(vcov(fit)))
  expect_equal(rownames(ci), names(b))
  expect_equal(colnames(ci), c("2.5 %", "97.5 %"))
  expect_equal(unname(ci[, 1]), unname(b - qnorm(0.975) * se))
  expect_equal(unname(ci[, 2]), unname(b + qnorm(0.975) * se))
  expect_equal(ci, confint.default(fit))
  expect_equal(unname(confint(fit, parm = "positive:crp", level = 0.9)),
               unname(confint.default(fit, parm = "positive:crp", level = 0.9)))
  expect_equal(colnames(confint(fit, level = 0.9)), c("5.0 %", "95.0 %"))
  expect_equal(confint(fit, parm = 2:3), ci[2:3, ])
  expect_error(confint(fit, parm = "nope"), "coefficients")
  expect_error(confint(fit, level = 2), "between 0 and 1")
})

test_that("confint.twopart uses the joint covariance of a correlated fit", {
  copd <- methods_data("copd")
  fit <- twopart(cost ~ treatment + crp, data = copd, index = c("id", "visit"),
                 model = "random", positive = "lognormal", correlated = TRUE)
  ci <- confint(fit)
  expect_equal(nrow(ci), length(coef(fit)))
  expect_equal(unname(ci[, 2] - ci[, 1]), unname(2 * qnorm(0.975) * sqrt(diag(vcov(fit)))))
})

test_that("logLik of an event study is that of its two-way fixed-effects fit", {
  es <- methods_fit("eventstudy")
  expect_equal(logLik(es), logLik(es$fit))
})

test_that("hurdle coefficient methods stack the two parts", {
  copd <- methods_data("copd")
  set.seed(5)
  copd$w <- rnorm(nrow(copd))
  fit <- suppressMessages(panglm_hurdle(exacerbations ~ crp + w, data = copd,
                                        index = c("id", "visit")))
  b <- coef(fit)
  expect_equal(names(b), c("zero:crp", "zero:w", "count:crp", "count:w"))
  expect_equal(unname(b[1:2]), unname(fit$zero$coefficients))
  expect_equal(unname(b[3:4]), unname(fit$count$coefficients))
  V <- vcov(fit)
  expect_equal(unname(V[1:2, 1:2]), unname(fit$zero$vcov))
  expect_equal(unname(V[3:4, 3:4]), unname(fit$count$vcov))
  expect_true(all(V[1:2, 3:4] == 0))
  expect_equal(nobs(fit), nrow(methods_data("copd")))
  expect_equal(confint(fit), confint.default(fit))
  skip_if_missing("generics")
  td <- generics::tidy(fit, conf.int = TRUE)
  expect_equal(td$component, c("zero", "zero", "count", "count"))
  expect_equal(td$conf.low, unname(confint(fit)[, 1]))
  expect_equal(generics::glance(fit)$n.positive, fit$n_positive)
})
