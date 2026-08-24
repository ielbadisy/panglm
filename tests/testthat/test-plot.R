test_that("plot.panglm supports coefficient and stored-value diagnostics", {
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  data(copd)
  fit <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
                model = "within", family = "gaussian")

  testthat::expect_no_error(plot(fit, which = "coefficients"))
  testthat::expect_no_error(plot(fit, which = "residuals"))
  testthat::expect_no_error(plot(fit, which = "fitted"))
  testthat::expect_identical(invisible(plot(fit)), fit)
})

test_that("plot.panglm reports unavailable diagnostics", {
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  skip_if_missing("pglm")
  data(UnionWage, package = "pglm")
  fit <- suppressMessages(
    panglm(union ~ wage + exper + rural, data = UnionWage,
           index = c("id", "year"), model = "within", family = "binomial")
  )

  expect_error(plot(fit, which = "residuals"), "require stored fitted values")
})

test_that("plot.panglm_hurdle displays either or both components", {
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  data(copd)
  fit <- suppressMessages(
    panglm_hurdle(exacerbations ~ crp, data = copd,
                  index = c("id", "visit"))
  )

  expect_no_error(plot(fit, which = "zero"))
  expect_no_error(plot(fit, which = "count"))
  expect_no_error(plot(fit, which = "both"))
})
