test_that("within gaussian matches plm()", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  f <- panglm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "within", family = "gaussian")
  p <- plm::plm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "within")
  expect_equal(unname(coef(f)), unname(coef(p)), tolerance = 1e-6)
})

test_that("within poisson matches pglm() conditional MLE", {
  # pglm()'s match.call()/parent.frame() NSE chain is fragile inside
  # testthat::local() frames (reproducibly fails there even though the
  # same call works from a plain script), so the reference coefficients
  # below were captured from a direct (non-testthat) comparison run
  # rather than calling pglm() live here.
  skip_if_missing("plm")
  Grunfeld <- .panglm_test_data$Grunfeld_poisson
  f <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "within", family = "poisson")
  pglm_ref <- c(value = -2.805996e-05, capital = 2.672633e-04)
  expect_equal(unname(coef(f)), unname(pglm_ref), tolerance = 1e-4)
})
