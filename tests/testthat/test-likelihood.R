test_that("pooled Gaussian likelihood and information criteria match glm", {
  set.seed(201)
  d <- data.frame(id = rep(seq_len(20), each = 4), time = rep(seq_len(4), 20))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(nrow(d), sd = 0.8)
  fit <- panglm(y ~ x, d, index = c("id", "time"), model = "pooling",
                family = "gaussian")
  reference <- stats::glm(y ~ x, data = d, family = stats::gaussian())

  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(reference)),
               tolerance = 1e-8)
  expect_equal(attr(logLik(fit), "df"), attr(logLik(reference), "df"))
  expect_equal(stats::AIC(fit), stats::AIC(reference), tolerance = 1e-8)
  expect_equal(stats::BIC(fit), stats::BIC(reference), tolerance = 1e-8)
})

test_that("pooled count likelihoods include estimated shape parameters", {
  set.seed(202)
  d <- data.frame(id = rep(seq_len(20), each = 5), time = rep(seq_len(5), 20))
  d$x <- stats::rnorm(nrow(d))
  d$y_pois <- stats::rpois(nrow(d), exp(0.2 + 0.3 * d$x))
  d$y_nb <- stats::rnbinom(nrow(d), size = 3, mu = exp(0.2 + 0.3 * d$x))

  poisson_fit <- panglm(y_pois ~ x, d, index = c("id", "time"),
                        model = "pooling", family = "poisson")
  negbin_fit <- panglm(y_nb ~ x, d, index = c("id", "time"),
                       model = "pooling", family = "negbin")
  expect_identical(attr(logLik(poisson_fit), "df"), 2L)
  expect_identical(attr(logLik(negbin_fit), "df"), 3L)
})

test_that("conditional likelihoods count only retained regression parameters", {
  set.seed(203)
  d <- data.frame(id = rep(seq_len(30), each = 5), time = rep(seq_len(5), 30))
  d$x <- stats::rnorm(nrow(d))
  d$count <- stats::rpois(nrow(d), exp(0.2 + 0.4 * d$x))
  d$binary <- stats::rbinom(nrow(d), 1, stats::plogis(0.5 * d$x))

  poisson_fit <- panglm(count ~ x, d, index = c("id", "time"),
                        model = "within", family = "poisson")
  binomial_fit <- suppressMessages(
    panglm(binary ~ x, d, index = c("id", "time"),
           model = "within", family = "binomial")
  )
  expect_identical(attr(logLik(poisson_fit), "df"), 1L)
  expect_identical(attr(logLik(binomial_fit), "df"), 1L)
})

test_that("full fixed-effect and random-effect likelihoods count nuisance parameters", {
  set.seed(204)
  N <- 12; Tt <- 5
  d <- data.frame(id = rep(seq_len(N), each = Tt), time = rep(seq_len(Tt), N))
  d$x <- stats::rnorm(nrow(d))
  d$count <- stats::rnbinom(nrow(d), size = 4, mu = exp(0.3 + 0.4 * d$x))
  d$binary <- stats::rbinom(nrow(d), 1, stats::plogis(0.3 * d$x))

  fe_nb <- suppressMessages(
    panglm(count ~ x, d, index = c("id", "time"), model = "within",
           family = "negbin")
  )
  tw_nb <- panglm(count ~ x, d, index = c("id", "time"), model = "within",
                  family = "negbin", effect = "twoways")
  re_bin <- panglm(binary ~ x, d, index = c("id", "time"), model = "random",
                   family = "binomial", R = 9)

  expect_equal(attr(logLik(fe_nb), "df"), 1L + N + 1L)
  expect_equal(attr(logLik(tw_nb), "df"), 1L + N + Tt)
  expect_identical(attr(logLik(re_bin), "df"), 3L)
})

test_that("glance reports likelihood criteria and non-likelihood fits return NA", {
  skip_if_missing("generics")
  data(copd)
  pooled <- panglm(exacerbations ~ crp, copd, index = c("id", "visit"),
                   model = "pooling", family = "poisson")
  within_gaussian <- panglm(fev1 ~ crp, copd, index = c("id", "visit"),
                            model = "within", family = "gaussian")
  pooled_glance <- generics::glance(pooled)
  within_glance <- generics::glance(within_gaussian)

  expect_equal(pooled_glance$AIC, stats::AIC(pooled))
  expect_equal(pooled_glance$BIC, stats::BIC(pooled))
  expect_identical(pooled_glance$df.logLik, attr(logLik(pooled), "df"))
  expect_true(is.na(within_glance$AIC) && is.na(within_glance$BIC))
})
