test_that("random gaussian matches plm()", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  f <- panglm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "random", family = "gaussian")
  p <- plm::plm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "random")
  expect_equal(unname(coef(f)), unname(coef(p)), tolerance = 1e-6)
})

test_that("random poisson (Poisson-Gamma) recovers known dispersion and matches pglm()", {
  # See the note in test-within.R: pglm()'s NSE is unreliable inside
  # testthat::local() frames, so we compare against reference values
  # captured from a direct (non-testthat) run instead of calling pglm()
  # here. The DGP uses a known Gamma(shape = 5) frailty as a sanity check.
  skip_if_missing("plm")
  Grunfeld <- .panglm_test_data$Grunfeld_negbin
  f <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "random", family = "poisson")
  pglm_ref <- c("(Intercept)" = 0.8005900, value = 0.0002702213, capital = -6.095023e-05)
  expect_equal(unname(coef(f)), unname(pglm_ref), tolerance = 1e-2)
  expect_equal(f$dispersion_param, 4.179586, tolerance = 1e-1)
  expect_equal(f$dispersion_param, 5, tolerance = 0.5) # recovers true Gamma shape = 5
})

test_that("random binomial (Gauss-Hermite quadrature) matches pglm()", {
  skip_if_missing("pglm")
  skip_if_missing("maxLik")
  suppressMessages({ library(pglm) })
  data(UnionWage, package = "pglm")
  f <- panglm(union ~ wage + exper + rural, data = UnionWage, index = c("id", "year"),
              model = "random", family = "binomial", R = 15)
  p <- pglm(union ~ wage + exper + rural, UnionWage, index = c("id", "year"),
            family = binomial("logit"), model = "random", method = "bfgs",
            print.level = 0, R = 15)
  expect_equal(unname(coef(f)), unname(coef(p)[1:4]), tolerance = 1e-3)
  expect_equal(f$sigma, unname(tail(coef(p), 1)), tolerance = 1e-2)
  expect_equal(f$loglik, as.numeric(logLik(p)), tolerance = 1e-3)
})
