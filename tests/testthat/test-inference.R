test_that("confint.panglm matches confint.default(glm) for pooled models", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  f <- panglm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "pooling", family = "gaussian")
  g <- glm(inv ~ value + capital, data = Grunfeld, family = gaussian())
  expect_equal(unname(confint(f)), unname(confint.default(g)), tolerance = 1e-6)
})

test_that("panglm_dispersiontest detects simulated overdispersion and equidispersion correctly", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")

  set.seed(1)
  Grunfeld$count_pois <- rpois(nrow(Grunfeld), lambda = exp(0.5 + 0.0002 * Grunfeld$value))
  f_pois <- panglm(count_pois ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "pooling", family = "poisson")
  d_pois <- panglm_dispersiontest(f_pois)
  expect_true(d_pois$ratio < 1.5) # no real overdispersion in the DGP

  set.seed(1)
  Grunfeld$count_nb <- rnbinom(nrow(Grunfeld), size = 0.5, mu = exp(0.5 + 0.0002 * Grunfeld$value))
  f_nb <- panglm(count_nb ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "pooling", family = "poisson")
  d_nb <- panglm_dispersiontest(f_nb)
  expect_true(d_nb$ratio > 2) # strong overdispersion (true model is NB with small size)
  expect_true(d_nb$p.value < 0.01)
})

test_that("panglm_dispersiontest errors informatively when fitted.values are unavailable", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  set.seed(1)
  Grunfeld$count <- rpois(nrow(Grunfeld), lambda = exp(0.5 + 0.0002 * Grunfeld$value))
  f <- suppressWarnings(panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "random", family = "poisson"))
  expect_error(panglm_dispersiontest(f), "fitted.values are not available")
})
