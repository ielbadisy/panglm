test_that("pooled gaussian matches glm()", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  f <- panglm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "pooling", family = "gaussian")
  g <- glm(inv ~ value + capital, data = Grunfeld, family = gaussian())
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-6)
})

test_that("pooled poisson matches glm()", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  set.seed(1)
  Grunfeld$count <- rpois(nrow(Grunfeld), lambda = exp(0.5 + 0.01 * Grunfeld$value / 100))
  f <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "pooling", family = "poisson")
  g <- glm(count ~ value + capital, data = Grunfeld, family = poisson())
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-6)
})

test_that("pooled binomial matches glm()", {
  skip_if_missing("pglm")
  data(UnionWage, package = "pglm")
  f <- panglm(union ~ wage + exper + rural, data = UnionWage, index = c("id", "year"), model = "pooling", family = "binomial")
  g <- glm(union ~ wage + exper + rural, data = UnionWage, family = binomial())
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-6)
})
