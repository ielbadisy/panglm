test_that("pooled negbin (real NB2, not aliased to Poisson) matches MASS::glm.nb exactly", {
  skip_if_missing("plm")
  skip_if_missing("MASS")
  data(Grunfeld, package = "plm")
  set.seed(1)
  Grunfeld$count <- rnbinom(nrow(Grunfeld), size = 3, mu = exp(0.5 + 0.0002 * Grunfeld$value))

  f <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "pooling", family = "negbin")
  g <- MASS::glm.nb(count ~ value + capital, data = Grunfeld)

  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-4)
  expect_equal(f$theta, g$theta, tolerance = 1e-3)
  expect_equal(f$loglik, as.numeric(stats::logLik(g)), tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(vcov(g)))), tolerance = 1e-3)
})

test_that("pooled negbin is NOT the same fit as pooled poisson (regression guard for the aliasing bug)", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  set.seed(1)
  Grunfeld$count <- rnbinom(nrow(Grunfeld), size = 3, mu = exp(0.5 + 0.0002 * Grunfeld$value))

  f_nb <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "pooling", family = "negbin")
  f_pois <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "pooling", family = "poisson")
  expect_false(isTRUE(all.equal(unname(coef(f_nb)), unname(coef(f_pois)))))
  expect_false(isTRUE(all.equal(f_nb$loglik, f_pois$loglik)))
})

test_that("within negbin (Allison-Waterman dummy FE) matches fixest::fenegbin exactly", {
  skip_if_missing("plm")
  skip_if_missing("fixest")
  data(Grunfeld, package = "plm")
  set.seed(2)
  u <- rep(rnorm(10, sd = 0.5), each = 20)
  Grunfeld$count <- rnbinom(nrow(Grunfeld), size = 4, mu = exp(1 + 0.0002 * Grunfeld$value + u))

  f <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "within", family = "negbin")
  fx <- fixest::fenegbin(count ~ value + capital | firm, data = Grunfeld, notes = FALSE)

  expect_equal(unname(coef(f)), unname(coef(fx)), tolerance = 1e-4)
  expect_equal(f$loglik, as.numeric(stats::logLik(fx)), tolerance = 1e-3)
  expect_length(f$individual_effects, 10)
})

test_that("within negbin is NOT the same fit as within poisson (regression guard for the aliasing bug)", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  set.seed(2)
  u <- rep(rnorm(10, sd = 0.5), each = 20)
  Grunfeld$count <- rnbinom(nrow(Grunfeld), size = 4, mu = exp(1 + 0.0002 * Grunfeld$value + u))

  f_nb <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "within", family = "negbin")
  f_pois <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "within", family = "poisson")
  expect_false(isTRUE(all.equal(unname(coef(f_nb)), unname(coef(f_pois)))))
})

test_that("random negbin (2-parameter beta-negative-binomial) recovers known simulated parameters", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  set.seed(3)
  Grunfeld$u <- rep(rgamma(10, shape = 6, rate = 6), each = 20)
  Grunfeld$count <- rnbinom(nrow(Grunfeld), size = 8, mu = Grunfeld$u * exp(1 + 0.0002 * Grunfeld$value))

  f <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "random", family = "negbin")
  expect_true(is.finite(f$loglik))
  expect_true(f$negbin_a > 0 && f$negbin_b > 0)
  # sanity: coefficients in the right ballpark of the simulated truth
  expect_equal(unname(coef(f)["value"]), 0.0002, tolerance = 1)
})

test_that("random negbin is NOT the same fit as random poisson (regression guard for the aliasing bug)", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  set.seed(3)
  Grunfeld$u <- rep(rgamma(10, shape = 6, rate = 6), each = 20)
  Grunfeld$count <- rnbinom(nrow(Grunfeld), size = 8, mu = Grunfeld$u * exp(1 + 0.0002 * Grunfeld$value))

  f_nb <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "random", family = "negbin")
  f_pois <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "random", family = "poisson")
  expect_false(isTRUE(all.equal(f_nb$loglik, f_pois$loglik)))
})
