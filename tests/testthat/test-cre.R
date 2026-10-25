test_that("Mundlak pooled and random fits reproduce the within slope (gaussian)", {
  data(copd, package = "panglm")
  fe <- panglm(fev1 ~ crp, copd, c("id", "visit"), model = "within")
  po <- panglm(fev1 ~ treatment + crp, copd, c("id", "visit"), model = "pooling", cre = TRUE)
  re <- panglm(fev1 ~ treatment + crp, copd, c("id", "visit"), model = "random", cre = TRUE)
  expect_equal(unname(coef(po)["crp"]), unname(coef(fe)), tolerance = 1e-8)
  expect_equal(unname(coef(re)["crp"]), unname(coef(fe)), tolerance = 1e-6)
  expect_true("mean_crp" %in% names(coef(po)))
  expect_false("mean_treatment" %in% names(coef(po)))
})

test_that("CRE works for nonlinear families and reports a Mundlak test", {
  data(copd, package = "panglm")
  f <- panglm(hospitalized ~ treatment + age + crp, copd, c("id", "visit"),
              model = "random", family = "probit", cre = TRUE)
  s <- summary(f)
  expect_false(is.null(s$cre_test))
  expect_true(s$cre_test$p.value >= 0 && s$cre_test$p.value <= 1)
  expect_output(print(s), "Mundlak")
  expect_equal(unname(predict(f, copd, type = "response")[order(copd$id, copd$visit)]),
               unname(fitted(f)), tolerance = 1e-8)
  expect_error(panglm(fev1 ~ crp, copd, c("id", "visit"), model = "within", cre = TRUE))
})
