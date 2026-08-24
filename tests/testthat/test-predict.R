test_that("in-sample fitted values retain the caller's row order", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  set.seed(101)
  shuffled <- Grunfeld[sample(seq_len(nrow(Grunfeld))), ]
  fit <- panglm(inv ~ value + capital, shuffled,
                index = c("firm", "year"), model = "pooling",
                family = "gaussian")
  reference <- stats::glm(inv ~ value + capital, data = shuffled)

  expect_equal(unname(fitted(fit)), unname(fitted(reference)), tolerance = 1e-8)
  expect_equal(predict(fit, type = "response"), fitted(fit))
  expect_equal(unname(residuals(fit)), shuffled$inv - fitted(fit))
})

test_that("one-way fixed-effect prediction includes known individual effects", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  fit <- panglm(inv ~ value + capital, Grunfeld,
                index = c("firm", "year"), model = "within",
                family = "gaussian")

  expect_equal(predict(fit, Grunfeld, type = "response"), fitted(fit),
               tolerance = 1e-8)
  unseen <- Grunfeld[1:2, ]
  unseen$firm <- 999
  expect_error(predict(fit, unseen), "individual levels")
  expected_component <- as.numeric(
    stats::model.matrix(~ value + capital, unseen)[, c("value", "capital")] %*%
      coef(fit)
  )
  expect_equal(
    predict(fit, unseen, type = "link", allow.new.levels = TRUE),
    expected_component
  )
})

test_that("two-way fixed-effect prediction includes individual and time effects", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  fit <- panglm(inv ~ value + capital, Grunfeld,
                index = c("firm", "year"), model = "within",
                family = "gaussian", effect = "twoways")

  expect_equal(predict(fit, Grunfeld, type = "response"), fitted(fit),
               tolerance = 1e-7)
  unseen <- Grunfeld[1, ]
  unseen$year <- 3000
  expect_error(predict(fit, unseen), "time levels")
})

test_that("conditional binomial fitted values preserve group success totals", {
  set.seed(102)
  d <- data.frame(
    id = rep(seq_len(30), each = 5),
    time = rep(seq_len(5), 30),
    x = stats::rnorm(150)
  )
  d$y <- stats::rbinom(nrow(d), 1, stats::plogis(0.8 * d$x))
  fit <- suppressMessages(
    panglm(y ~ x, d, index = c("id", "time"), model = "within",
           family = "binomial")
  )
  probabilities <- fitted(fit)

  expect_true(all(probabilities >= 0 & probabilities <= 1))
  expect_equal(
    as.numeric(rowsum(probabilities, d$id)),
    as.numeric(rowsum(d$y, d$id)), tolerance = 1e-8
  )
  expect_error(predict(fit, d), "not defined for exact conditional")
})

test_that("random Poisson predictions are marginal population means", {
  set.seed(103)
  d <- data.frame(id = rep(seq_len(20), each = 5), time = rep(seq_len(5), 20))
  d$x <- stats::rnorm(nrow(d))
  frailty <- rep(stats::rgamma(20, shape = 4, rate = 4), each = 5)
  d$y <- stats::rpois(nrow(d), frailty * exp(0.2 + 0.3 * d$x))
  fit <- suppressWarnings(
    panglm(y ~ x, d, index = c("id", "time"), model = "random",
           family = "poisson")
  )

  expected <- exp(as.numeric(stats::model.matrix(~ x, d) %*% coef(fit)))
  expect_equal(fitted(fit), expected, tolerance = 1e-8)
  expect_equal(predict(fit, d, type = "response"), expected, tolerance = 1e-8)
})

test_that("random binomial predictions integrate over the random intercept", {
  set.seed(104)
  d <- data.frame(id = rep(seq_len(25), each = 5), time = rep(seq_len(5), 25))
  d$x <- stats::rnorm(nrow(d))
  intercept <- rep(stats::rnorm(25, sd = 0.7), each = 5)
  d$y <- stats::rbinom(nrow(d), 1, stats::plogis(-0.2 + 0.5 * d$x + intercept))
  fit <- panglm(y ~ x, d, index = c("id", "time"), model = "random",
                family = "binomial", R = 11)

  predictions <- predict(fit, d, type = "response")
  expect_equal(predictions, fitted(fit), tolerance = 1e-10)
  expect_true(all(predictions > 0 & predictions < 1))
  expect_false(isTRUE(all.equal(
    predictions,
    stats::plogis(predict(fit, d, type = "link")), tolerance = 1e-6
  )))
})
