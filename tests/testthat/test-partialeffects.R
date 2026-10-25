test_that("pooled probit APEs match marginaleffects", {
  skip_if_missing("marginaleffects")
  data(copd, package = "panglm")
  f <- panglm(hospitalized ~ treatment + age + crp, copd, c("id", "visit"), family = "probit")
  g <- glm(hospitalized ~ treatment + age + crp, binomial("probit"), copd)
  a <- partialeffects(f, variables = c("age", "crp"))
  m <- marginaleffects::avg_slopes(g, variables = c("age", "crp"))
  m <- as.data.frame(m)
  m <- m[match(a$term, m$term), ]
  expect_equal(a$estimate, m$estimate, tolerance = 1e-4)
  expect_equal(a$std.error, m$std.error, tolerance = 1e-3)
})

test_that("pooled linear and FE Poisson APEs equal analytic values", {
  data(copd, package = "panglm")
  f <- panglm(fev1 ~ treatment + crp, copd, c("id", "visit"))
  a <- partialeffects(f)
  expect_equal(a$estimate, unname(coef(f)[a$term]), tolerance = 1e-6)
  expect_equal(a$std.error, unname(sqrt(diag(vcov(f)))[a$term]), tolerance = 1e-4)

  fe <- panglm(exacerbations ~ crp, copd, c("id", "visit"), model = "within", family = "poisson")
  a <- partialeffects(fe)
  expect_equal(a$estimate, mean(fitted(fe)) * unname(coef(fe)), tolerance = 1e-5)
  expect_s3_class(a, "panglm_ape")
  expect_output(print(a), "Average partial effects")
})

test_that("APEs for random CRE and gamma fits run and are finite", {
  d <- make_cost_data()
  f <- panglm(hospitalized ~ treatment + age + crp, d, c("id", "visit"),
              model = "random", family = "probit", cre = TRUE)
  a <- partialeffects(f)
  expect_true(all(is.finite(a$estimate)) && all(a$std.error > 0))
  expect_false(any(grepl("^mean_", a$term)))
  g <- panglm(cost ~ treatment + crp, d, c("id", "visit"), model = "random", family = "gamma")
  expect_true(all(is.finite(partialeffects(g)$estimate)))
})

test_that("conditional logit APEs are refused", {
  data(copd, package = "panglm")
  f <- panglm(hospitalized ~ crp, copd, c("id", "visit"), model = "within", family = "binomial")
  expect_error(partialeffects(f), "not identified")
})
