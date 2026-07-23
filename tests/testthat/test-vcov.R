test_that("pooled HC1/cluster vcov match sandwich", {
  skip_if_missing("plm")
  skip_if_missing("sandwich")
  data(Grunfeld, package = "plm")
  f <- panglm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"),
              model = "pooling", family = "gaussian")
  g <- glm(inv ~ value + capital, data = Grunfeld, family = gaussian())

  expect_equal(sqrt(diag(vcov(f, type = "HC1"))),
               sqrt(diag(sandwich::vcovHC(g, type = "HC1"))),
               tolerance = 1e-8, ignore_attr = TRUE)

  expect_equal(sqrt(diag(vcov(f, type = "cluster"))),
               sqrt(diag(sandwich::vcovCL(g, cluster = Grunfeld$firm, type = "HC1"))),
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("within (FE) cluster vcov matches plm's Stata-style (sss) correction", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  f <- panglm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"),
              model = "within", family = "gaussian")
  p <- plm::plm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "within")

  expect_equal(sqrt(diag(vcov(f, type = "cluster"))),
               sqrt(diag(plm::vcovHC(p, method = "arellano", type = "sss", cluster = "group"))),
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("panglm(..., vcov = 'cluster') matches post-hoc vcov(fit, type = 'cluster')", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  f1 <- panglm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"),
               model = "pooling", family = "gaussian", vcov = "cluster")
  f2 <- panglm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"),
               model = "pooling", family = "gaussian")
  expect_equal(f1$vcov, vcov(f2, type = "cluster"))
  expect_identical(f1$vcov_type, "cluster")
})

test_that("two-way FE gaussian matches plm(effect = 'twoways')", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  f <- panglm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"),
              model = "within", family = "gaussian", effect = "twoways")
  p <- plm::plm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"),
                model = "within", effect = "twoways")
  expect_equal(unname(coef(f)), unname(coef(p)), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(vcov(p)))), tolerance = 1e-6)
  expect_equal(f$df.residual, unname(df.residual(p)))
})

test_that("two-way FE gaussian cluster vcov matches plm's Stata-style (sss) correction", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  f <- panglm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"),
              model = "within", family = "gaussian", effect = "twoways")
  p <- plm::plm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"),
                model = "within", effect = "twoways")
  expect_equal(sqrt(diag(vcov(f, type = "cluster"))),
               sqrt(diag(plm::vcovHC(p, method = "arellano", type = "sss", cluster = "group"))),
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("robust/cluster vcov errors informatively for within negbin (not yet implemented)", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  set.seed(2)
  u <- rep(rnorm(10, sd = 0.5), each = 20)
  Grunfeld$count <- rnbinom(nrow(Grunfeld), size = 4, mu = exp(1 + 0.0002 * Grunfeld$value + u))
  f <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "within", family = "negbin")
  expect_error(vcov(f, type = "cluster"), "not yet implemented")
})

test_that("Hausman test matches plm::phtest", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  fe <- panglm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "within", family = "gaussian")
  re <- panglm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "random", family = "gaussian")
  h <- panglm_hausman(fe, re)
  ph <- plm::phtest(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"))
  expect_equal(h$statistic, unname(ph$statistic), tolerance = 1e-6)
  expect_equal(h$parameter, unname(ph$parameter))
  expect_equal(h$p.value, unname(ph$p.value), tolerance = 1e-6)
})

test_that("tidy()/glance() dispatch and return the expected shape", {
  skip_if_missing("plm")
  skip_if_missing("generics")
  data(Grunfeld, package = "plm")
  f <- panglm(inv ~ value + capital, data = Grunfeld, index = c("firm", "year"),
              model = "pooling", family = "gaussian")
  td <- generics::tidy(f)
  expect_named(td, c("term", "estimate", "std.error", "statistic", "p.value"))
  expect_equal(nrow(td), 3)
  gl <- generics::glance(f)
  expect_equal(nrow(gl), 1)
  expect_true(all(c("model", "family", "nobs") %in% names(gl)))
})
