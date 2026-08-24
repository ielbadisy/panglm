test_that("factor terms and interactions match pooled glm", {
  set.seed(301)
  d <- data.frame(
    id = rep(seq_len(30), each = 4),
    time = rep(seq_len(4), 30),
    group = factor(rep(c("control", "active"), each = 60)),
    x = stats::rnorm(120)
  )
  d$y <- 1 + 0.4 * d$x + 0.6 * (d$group == "active") -
    0.3 * d$x * (d$group == "active") + stats::rnorm(120)

  fit <- panglm(y ~ x * group, d, index = c("id", "time"),
                model = "pooling", family = "gaussian")
  reference <- stats::glm(y ~ x * group, data = d)
  expect_equal(coef(fit), coef(reference), tolerance = 1e-8)
  expect_equal(predict(fit, d[1:8, ]), unname(predict(reference, d[1:8, ])),
               tolerance = 1e-8)
})

test_that("rank-deficient columns are reported and omitted", {
  set.seed(302)
  d <- data.frame(id = rep(seq_len(20), each = 4), time = rep(seq_len(4), 20))
  d$x <- stats::rnorm(nrow(d))
  d$x_duplicate <- d$x
  d$constant <- rep(stats::rnorm(20), each = 4)
  d$y <- d$x + rep(stats::rnorm(20), each = 4) + stats::rnorm(nrow(d))

  expect_warning(
    fit <- panglm(y ~ x + x_duplicate + constant, d,
                  index = c("id", "time"), model = "within",
                  family = "gaussian"),
    "dropping unidentified"
  )
  expect_identical(names(coef(fit)), "x")
  expect_setequal(fit$aliased, c("x_duplicate", "constant"))

  expect_error(
    panglm(y ~ constant, d, index = c("id", "time"), model = "within",
           family = "gaussian"),
    "no covariates identified"
  )
})

test_that("na.exclude restores fitted values and residuals to input length", {
  set.seed(303)
  d <- data.frame(id = rep(seq_len(15), each = 4), time = rep(seq_len(4), 15))
  d$x <- stats::rnorm(nrow(d))
  d$y <- 0.5 * d$x + stats::rnorm(nrow(d))
  d$x[c(3, 17)] <- NA_real_
  fit <- panglm(y ~ x, d, index = c("id", "time"), model = "pooling",
                family = "gaussian", na.action = stats::na.exclude)

  expect_length(fitted(fit), nrow(d))
  expect_length(residuals(fit), nrow(d))
  expect_true(all(is.na(fitted(fit)[c(3, 17)])))
  expect_true(all(is.na(residuals(fit)[c(3, 17)])))
  expect_identical(nobs(fit), nrow(d) - 2L)
})

test_that("panel indexes are validated before estimation", {
  d <- data.frame(id = c(1, 1, 2, 2), time = c(1, 2, 1, 2), x = 1:4, y = 2:5)
  expect_error(panglm(y ~ x, d, index = "missing"), "not found")
  expect_error(panglm(y ~ x, d, index = c("id", "id")), "distinct")

  missing_index <- d
  missing_index$time[2] <- NA
  expect_error(
    panglm(y ~ x, missing_index, index = c("id", "time")),
    "must not contain missing"
  )

  duplicate_index <- rbind(d, d[1, ])
  expect_error(
    panglm(y ~ x, duplicate_index, index = c("id", "time")),
    "unique observation"
  )
})

test_that("unsupported weights and offsets are rejected explicitly", {
  set.seed(304)
  d <- data.frame(id = rep(seq_len(10), each = 4), time = rep(seq_len(4), 10))
  d$x <- stats::rnorm(nrow(d))
  d$y <- stats::rpois(nrow(d), exp(0.2 + 0.3 * d$x))

  unit_fit <- panglm(y ~ x, d, index = c("id", "time"), model = "pooling",
                     family = "poisson", weights = rep(1, nrow(d)))
  plain_fit <- panglm(y ~ x, d, index = c("id", "time"), model = "pooling",
                      family = "poisson")
  expect_equal(coef(unit_fit), coef(plain_fit))
  expect_error(
    panglm(y ~ x, d, index = c("id", "time"), weights = seq_len(nrow(d))),
    "weights are not supported"
  )
  expect_error(
    panglm(y ~ x + offset(x), d, index = c("id", "time")),
    "offsets are not supported"
  )
  expect_error(
    panglm(y ~ x, d, index = c("id", "time"), offset = d$x),
    "offsets are not supported"
  )
})

test_that("invalid response domains fail before numerical fitting", {
  d <- data.frame(id = rep(1:3, each = 3), time = rep(1:3, 3), x = 1:9)
  d$y <- c(-1, rep(1, 8))
  expect_error(
    panglm(y ~ x, d, index = c("id", "time"), family = "poisson"),
    "must be nonnegative"
  )
})
