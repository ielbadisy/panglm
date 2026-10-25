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
    panglm(y ~ x, d, index = c("id", "time"), model = "within",
           family = "poisson", weights = seq_len(nrow(d))),
    "only supported for model = 'pooling'"
  )
  expect_error(
    panglm(y ~ x, d, index = c("id", "time"), weights = -seq_len(nrow(d))),
    "nonnegative"
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

test_that("a factor individual index gives the same fit as its character codes", {
  set.seed(303)
  d <- data.frame(id = rep(c("b", "a", "c", "d"), each = 5), time = rep(1:5, 4))
  d$x <- stats::rnorm(20)
  d$y <- d$x + rep(stats::rnorm(4), each = 5) + stats::rnorm(20)
  a <- panglm(y ~ x, d, index = c("id", "time"), model = "within")
  d$id <- factor(d$id, levels = c("d", "c", "b", "a"))
  b <- panglm(y ~ x, d, index = c("id", "time"), model = "within")
  expect_equal(coef(b), coef(a), tolerance = 1e-12)
  expect_equal(nobs(b), 20L)
  expect_equal(sum(b$group_size), 20L)
  expect_length(b$group_size, 4L)
})

lag_panel <- function(seed = 304) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(30), each = 6), time = rep(2001:2006, 30))
  d$x <- stats::rnorm(nrow(d))
  d$z <- stats::rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + rep(stats::rnorm(30), each = 6) + stats::rnorm(nrow(d))
  d <- d[-c(3, 40, 41), ]  # gaps inside some panels
  d[sample(nrow(d)), ]     # rows out of panel order
}

# shifted columns built by hand on the grid of sorted distinct times
lag_by_hand <- function(d, v, k) {
  key <- paste(d$id, d$time)
  d[[v]][match(paste(d$id, d$time - k), key)]
}

test_that("lag(), lead() and diff() in a formula are taken within individuals", {
  d <- lag_panel()
  d$x1 <- lag_by_hand(d, "x", 1); d$x2 <- lag_by_hand(d, "x", 2)
  d$zf <- lag_by_hand(d, "z", -1); d$dx <- d$x - d$x1
  for (m in c("pooling", "within", "random")) {
    a <- panglm(y ~ lag(x, 1:2) + lead(z) + diff(x), d, index = c("id", "time"), model = m)
    b <- panglm(y ~ x1 + x2 + zf + dx, d, index = c("id", "time"), model = m)
    expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-10)
    expect_equal(nobs(a), nobs(b))
    expect_equal(unname(sqrt(diag(vcov(a)))), unname(sqrt(diag(vcov(b)))), tolerance = 1e-10)
  }
  expect_true(all(c("lag(x, 1:2)1", "lag(x, 1:2)2") %in% names(coef(a))))
  # predictions on new rows reach back into the estimation data
  a <- panglm(y ~ lag(x) + z, d, index = c("id", "time"), model = "pooling")
  b <- panglm(y ~ x1 + z, d, index = c("id", "time"), model = "pooling")
  nd <- d[d$time == 2006, ]
  expect_equal(predict(a, newdata = nd), predict(b, newdata = nd), tolerance = 1e-10)
  expect_error(predict(a, newdata = nd[, c("x", "z")]), "panel index")
  pa <- partialeffects(a)
  expect_true(all(is.finite(unlist(pa[vapply(pa, is.numeric, logical(1))]))))
  # a lag beside dynamic = TRUE reaches the wave the dynamic model drops
  a <- panglm(y ~ lag(x) + z, d, index = c("id", "time"), model = "random", dynamic = TRUE)
  b <- panglm(y ~ x1 + z, d, index = c("id", "time"), model = "random", dynamic = TRUE)
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-8)
})

test_that("panel lags in formulas match plm", {
  skip_on_cran()
  skip_if_missing("plm")
  data("EmplUK", package = "plm")
  f <- log(emp) ~ lag(log(emp), 1) + lag(log(wage), 0:1) + diff(log(capital))
  p <- plm::plm(f, data = EmplUK, index = c("firm", "year"), model = "within")
  g <- panglm(f, data = EmplUK, index = c("firm", "year"), model = "within")
  expect_equal(unname(coef(g)), unname(coef(p)), tolerance = 1e-8)
  expect_equal(unname(sqrt(diag(vcov(g)))), unname(sqrt(diag(vcov(p)))), tolerance = 1e-8)
  expect_equal(nobs(g), nobs(p))
  data("Produc", package = "plm")
  f <- log(gsp) ~ lag(log(pcap)) + unemp | lag(log(pcap), 2) + unemp
  iv <- paniv(f, data = Produc, index = c("state", "year"), model = "within")
  piv <- plm::plm(f, data = Produc, index = c("state", "year"), model = "within")
  expect_equal(unname(coef(iv)), unname(coef(piv)), tolerance = 1e-8)
})

test_that("estimators that cannot shift within individuals refuse panel calls", {
  d <- lag_panel()
  expect_error(panglm(y ~ x, d, index = "id"), NA)
  expect_error(panglm(y ~ lag(x), d, index = "id"), "index = c")
  q <- panqr(y ~ lag(x) + z, d, index = c("id", "time"), model = "pooling")
  expect_equal(nobs(q), sum(!is.na(lag_by_hand(d, "x", 1))))
  expect_error(stats::model.frame(formula(q), d), "cannot rebuild")
  d$yc <- stats::rpois(nrow(d), exp(0.3 + 0.4 * d$x))
  d$x1 <- lag_by_hand(d, "x", 1)
  a <- panfmm(yc ~ lag(x), d, index = c("id", "time"), k = 2)
  b <- panfmm(yc ~ x1, d, index = c("id", "time"), k = 2)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-8)
  expect_error(predict(a, newdata = d), "cannot rebuild")
  expect_error(panhurdle(yc ~ lag(x), d, index = c("id", "time")), "not supported")
  dd <- d; dd$time[2] <- dd$time[1]; dd$id[2] <- dd$id[1]
  expect_error(panglm(y ~ lag(x), dd, index = c("id", "time")), "one row per individual")
})
