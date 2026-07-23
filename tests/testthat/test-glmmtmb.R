# glmmTMB comparisons split into two kinds, deliberately not conflated:
#
# 1. Exact-match cases, where the two models are provably the same
#    likelihood written differently (e.g. Poisson FE conditional MLE ==
#    Poisson FE dummy-variable MLE, a well-known equivalence that does NOT
#    hold for negative binomial).
# 2. Ballpark/sanity-recovery cases, where panglm's random-effects models
#    use a Gamma/Beta conjugate mixing distribution and glmmTMB uses a
#    Normal random intercept with Laplace approximation -- genuinely
#    different statistical models. Exact agreement is neither expected nor
#    the right bar; both should simply recover the same simulated truth in
#    the same ballpark.

test_that("FE poisson (conditional MLE) matches glmmTMB dummy-variable MLE exactly", {
  skip_if_missing("plm")
  skip_if_missing("glmmTMB")
  data(Grunfeld, package = "plm")
  set.seed(1)
  Grunfeld$count <- rpois(nrow(Grunfeld), lambda = exp(0.5 + 0.01 * Grunfeld$value / 100))

  f <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "within", family = "poisson")
  gt <- suppressWarnings(glmmTMB::glmmTMB(count ~ value + capital + factor(firm), data = Grunfeld, family = poisson()))

  gt_coef <- glmmTMB::fixef(gt)$cond[c("value", "capital")]
  expect_equal(unname(coef(f)), unname(gt_coef), tolerance = 1e-4)
})

test_that("random poisson (Poisson-Gamma) recovers similar slopes to glmmTMB's Poisson-lognormal GLMM", {
  skip_if_missing("plm")
  skip_if_missing("glmmTMB")
  data(Grunfeld, package = "plm")
  set.seed(42)
  Grunfeld$u <- rep(rgamma(10, shape = 5, rate = 5), each = 20)
  Grunfeld$count <- rpois(nrow(Grunfeld), lambda = Grunfeld$u * exp(1 + 0.0002 * Grunfeld$value))

  f <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "random", family = "poisson")
  gt <- suppressWarnings(glmmTMB::glmmTMB(count ~ value + capital + (1 | firm), data = Grunfeld, family = poisson()))

  gt_coef <- glmmTMB::fixef(gt)$cond[c("value", "capital")]
  # Not an exact-match test (different random-effect distributions): both
  # should be in the same ballpark and close to the simulated truth
  # (value = 0.0002, capital = 0).
  expect_equal(unname(coef(f)[c("value", "capital")]), unname(gt_coef), tolerance = 0.3)
  expect_equal(unname(coef(f)["value"]), 0.0002, tolerance = 1)
})

test_that("random negbin (beta-NB mixture) recovers similar slopes to glmmTMB's nbinom2 GLMM", {
  skip_if_missing("plm")
  skip_if_missing("glmmTMB")
  data(Grunfeld, package = "plm")
  set.seed(3)
  Grunfeld$u <- rep(rgamma(10, shape = 6, rate = 6), each = 20)
  Grunfeld$count <- rnbinom(nrow(Grunfeld), size = 8, mu = Grunfeld$u * exp(1 + 0.0002 * Grunfeld$value))

  f <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "random", family = "negbin")
  gt <- suppressWarnings(glmmTMB::glmmTMB(count ~ value + capital + (1 | firm), data = Grunfeld, family = glmmTMB::nbinom2()))

  gt_coef <- glmmTMB::fixef(gt)$cond[c("value", "capital")]
  expect_equal(unname(coef(f)[c("value", "capital")]), unname(gt_coef), tolerance = 0.3)
})
