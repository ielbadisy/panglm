test_that("between estimator matches plm and lm on the individual means", {
  skip_if_missing("plm")
  data("Grunfeld", package = "plm")
  f <- panglm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"),
              model = "between")
  p <- plm::plm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"),
                model = "between")
  expect_equal(coef(f), coef(p), tolerance = 1e-10)
  expect_equal(unname(vcov(f)), unname(vcov(p)), tolerance = 1e-10)
  expect_equal(f$df.residual, 7L)
  means <- aggregate(Grunfeld[c("inv", "value", "capital")], list(firm = Grunfeld$firm), mean)
  m <- lm(inv ~ value + capital, data = means)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(m)), tolerance = 1e-10)
})

test_that("between robust covariance is HC1 on the means and handles unbalanced panels", {
  skip_if_missing("sandwich")
  data(copd)
  set.seed(3)
  d <- copd[sort(sample(nrow(copd), 250)), ]
  f <- panglm(fev1 ~ crp + age, data = d, index = c("id", "visit"), model = "between",
              vcov = "HC1")
  means <- aggregate(d[c("fev1", "crp", "age")], list(id = d$id), mean)
  m <- lm(fev1 ~ crp + age, data = means)
  expect_equal(unname(coef(f)), unname(coef(m)), tolerance = 1e-10)
  expect_equal(unname(vcov(f)), unname(sandwich::vcovHC(m, type = "HC1")), tolerance = 1e-10)
  expect_equal(unname(vcov(f, type = "cluster")), unname(vcov(f, type = "HC1")))
  expect_equal(unname(predict(f, newdata = d[1:3, ])),
               unname(predict(m, newdata = d[1:3, ])), tolerance = 1e-10)
  expect_output(print(summary(f)), "between")
})

test_that("between estimator input checks", {
  data(copd)
  expect_error(panglm(exacerbations ~ crp, data = copd, index = c("id", "visit"),
                      model = "between", family = "poisson"), "gaussian")
  expect_error(panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
                      model = "between", cre = TRUE), "cre")
  expect_error(panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
                      model = "between", weights = rep(2, nrow(copd))), "weights")
  f <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"), model = "between")
  expect_error(vcov(f, type = "cluster", cluster = copd$id), "custom clusters")
})
