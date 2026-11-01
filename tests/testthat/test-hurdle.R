test_that("fit_truncated_poisson_dummy matches pscl::hurdle's count part (single group == pooled)", {
  skip_if_missing("pscl")
  set.seed(5)
  n <- 2000
  x1 <- rnorm(n)
  mu <- exp(0.3 + 0.6 * x1)
  p_zero <- stats::plogis(-0.5 - 0.4 * x1)
  is_zero <- stats::rbinom(n, 1, p_zero) == 1
  y <- ifelse(is_zero, 0, stats::rpois(n, mu))
  while (any(!is_zero & y == 0)) {
    idx <- which(!is_zero & y == 0)
    y[idx] <- stats::rpois(length(idx), mu[idx])
  }

  X <- matrix(x1, ncol = 1); colnames(X) <- "x1"
  group <- rep(1, n)
  fit <- fit_truncated_poisson_dummy(X[y > 0, , drop = FALSE], y[y > 0], group[y > 0], 100, 1e-10)

  h <- pscl::hurdle(y ~ x1, dist = "poisson", zero.dist = "binomial")
  ref <- stats::coef(h)

  expect_equal(unname(fit$individual_effects), unname(ref["count_(Intercept)"]), tolerance = 1e-6)
  expect_equal(unname(fit$coefficients["x1"]), unname(ref["count_x1"]), tolerance = 1e-6)
})

test_that("panhurdle() fits a panel with structural zeros and recovers the count-part slope", {
  set.seed(9)
  N <- 80; Tt <- 8
  d <- data.frame(id = rep(1:N, each = Tt), time = rep(1:Tt, N))
  d$x1 <- rnorm(N * Tt)
  a_i <- rep(rnorm(N, sd = 1.5), each = Tt)
  mu <- exp(0.3 + 0.5 * d$x1 + a_i)
  p_zero <- stats::plogis(1 - 1.2 * a_i)
  is_zero <- stats::rbinom(N * Tt, 1, p_zero) == 1
  y <- ifelse(is_zero, 0, stats::rpois(N * Tt, mu))
  while (any(!is_zero & y == 0)) {
    idx <- which(!is_zero & y == 0)
    y[idx] <- stats::rpois(length(idx), mu[idx])
  }
  d$y <- y

  fit <- panhurdle(y ~ x1, data = d, index = c("id", "time"))
  expect_s3_class(fit, "panglm_hurdle")
  expect_true(is.finite(fit$count$coefficients["x1"]))
  expect_equal(unname(fit$count$coefficients["x1"]), 0.5, tolerance = 0.3)
  expect_true(fit$n_positive < fit$n)
})

test_that("zero-truncated NB2 count part matches pscl::hurdle for one group", {
  skip_if_missing("pscl")
  set.seed(27)
  n <- 1200
  x1 <- rnorm(n)
  mu <- exp(0.2 + 0.5 * x1)
  y <- stats::rnbinom(n, size = 2.5, mu = mu)
  while (any(y == 0L)) {
    idx <- which(y == 0L)
    y[idx] <- stats::rnbinom(length(idx), size = 2.5, mu = mu[idx])
  }

  is_hurdle_zero <- stats::rbinom(n, 1, 0.2) == 1L
  y_observed <- y
  y_observed[is_hurdle_zero] <- 0L
  keep <- y_observed > 0L
  X <- matrix(x1[keep], ncol = 1, dimnames = list(NULL, "x1"))
  fit <- fit_truncated_negbin_dummy(X, y_observed[keep], rep(1, sum(keep)), 200, 1e-10)
  ref <- pscl::hurdle(y_observed ~ x1, dist = "negbin", zero.dist = "binomial")
  ref_coef <- stats::coef(ref)

  expect_equal(unname(fit$individual_effects),
               unname(ref_coef["count_(Intercept)"]), tolerance = 1e-4)
  expect_equal(unname(fit$coefficients["x1"]),
               unname(ref_coef["count_x1"]), tolerance = 1e-4)
  expect_equal(unname(fit$theta), unname(ref$theta), tolerance = 1e-3)
})

test_that("panhurdle() supports an NB2 positive-count component", {
  set.seed(28)
  N <- 30; Tt <- 6
  d <- data.frame(id = rep(seq_len(N), each = Tt), time = rep(seq_len(Tt), N))
  d$x1 <- rnorm(nrow(d))
  mu <- exp(0.3 + 0.4 * d$x1 + rep(rnorm(N, sd = 0.3), each = Tt))
  d$y <- stats::rnbinom(nrow(d), size = 3, mu = mu)

  fit <- suppressMessages(
    panhurdle(y ~ x1, d, index = c("id", "time"), family = "negbin")
  )
  expect_identical(fit$count_family, "negbin")
  expect_true(is.finite(fit$count$theta) && fit$count$theta > 0)
  expect_true(is.finite(fit$count$coefficients["x1"]))
})

test_that("NB2 count part reaches the Poisson limit when the counts are underdispersed", {
  # theta runs off to infinity here; the old log P(Y = 0) form turned NaN and
  # the fit stopped with "NAs are not allowed in subscripted assignments"
  set.seed(31)
  n <- 300
  x <- rnorm(n)
  y <- pmin(stats::rpois(n, exp(1 + 0.3 * x)), 3)
  keep <- y > 0
  g <- rep(1:30, each = 10)
  X <- matrix(x[keep], dimnames = list(NULL, "x"))
  nb <- suppressWarnings(fit_truncated_negbin_dummy(X, y[keep], g[keep], 100, 1e-10))
  po <- fit_truncated_poisson_dummy(X, y[keep], g[keep], 100, 1e-10)
  expect_true(nb$theta > 1e6)
  expect_equal(unname(nb$coefficients), unname(po$coefficients), tolerance = 1e-4)
  expect_equal(nb$loglik, po$loglik, tolerance = 1e-6)
})

hurdle_panel <- function(seed = 3, N = 40, Tt = 6) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(N), each = Tt), time = rep(seq_len(Tt), N))
  d$x1 <- rnorm(N * Tt)
  d$x2 <- rnorm(N * Tt)
  d$y <- stats::rnbinom(N * Tt, size = 1.5,
                        mu = exp(0.5 + 0.4 * d$x1 - 0.3 * d$x2 + rep(rnorm(N), each = Tt)))
  d
}

test_that("count part matches glmmTMB truncated counts with group dummies", {
  skip_if_missing("glmmTMB")
  d <- hurdle_panel()
  dp <- d[d$y > 0, ]
  ones <- tapply(dp$y, dp$id, function(v) all(v == 1))
  dp <- dp[!ones[as.character(dp$id)], ]
  for (fam in c("poisson", "negbin")) {
    fit <- suppressMessages(panhurdle(y ~ x1 + x2, d, c("id", "time"), family = fam))
    ref <- glmmTMB::glmmTMB(
      y ~ x1 + x2 + factor(id), data = dp,
      family = if (fam == "poisson") glmmTMB::truncated_poisson() else glmmTMB::truncated_nbinom2()
    )
    b <- glmmTMB::fixef(ref)$cond[c("x1", "x2")]
    se <- sqrt(diag(stats::vcov(ref)$cond))[c("x1", "x2")]
    expect_equal(unname(fit$count$coefficients), unname(b), tolerance = 1e-4)
    expect_equal(unname(sqrt(diag(fit$count$vcov))), unname(se), tolerance = 1e-4)
    expect_equal(fit$count$loglik, as.numeric(stats::logLik(ref)), tolerance = 1e-6)
    if (fam == "negbin") expect_equal(fit$count$theta, stats::sigma(ref), tolerance = 1e-4)
  }
})

test_that("zero part matches survival::clogit", {
  skip_if_missing("survival")
  suppressMessages(library(survival))
  d <- hurdle_panel()
  fit <- suppressMessages(panhurdle(y ~ x1 + x2, d, c("id", "time")))
  d$pos <- as.numeric(d$y > 0)
  ref <- clogit(pos ~ x1 + x2 + strata(id), data = d, method = "exact")
  expect_equal(unname(fit$zero$coefficients), unname(stats::coef(ref)), tolerance = 1e-6)
  expect_equal(fit$zero$loglik, ref$loglik[2], tolerance = 1e-6)
})

test_that("groups whose positive counts are all 1 are dropped from the count part", {
  d <- hurdle_panel()
  expect_message(fit <- panhurdle(y ~ x1 + x2, d, c("id", "time")), "all 1 dropped")
  pos <- d[d$y > 0, ]
  n_ones <- sum(tapply(pos$y, pos$id, function(v) all(v == 1)))
  expect_gt(n_ones, 0)
  expect_identical(fit$count$n_dropped_groups, n_ones)
  expect_length(fit$count$fitted.values, fit$n_positive)
  expect_true(all(is.finite(stats::coef(fit))))

  d$y <- as.numeric(d$y > 0)
  expect_error(suppressMessages(panhurdle(y ~ x1, d, c("id", "time"))), "identified")
})

test_that("panhurdle() fits the docvisits panel in seconds", {
  skip_on_cran()
  data(docvisits, package = "panglm", envir = environment())
  elapsed <- system.time(
    fit <- suppressMessages(panhurdle(docvis ~ age + hhninc + hsat, docvisits,
                                      c("id", "year"), family = "negbin"))
  )[["elapsed"]]
  expect_lt(elapsed, 20)
  expect_identical(fit$count$convergence, 0L)
  expect_true(all(is.finite(sqrt(diag(stats::vcov(fit))))))
  expect_true(fit$count$theta > 1 && fit$count$theta < 10)
})

test_that("time-invariant regressors are dropped with a warning", {
  data(copd, package = "panglm", envir = environment())
  ref <- suppressMessages(panhurdle(exacerbations ~ crp + fev1, copd, c("id", "visit")))
  expect_warning(
    fit <- suppressMessages(panhurdle(exacerbations ~ crp + age + fev1, copd, c("id", "visit"))),
    "time-invariant regressor\\(s\\) dropped.*age"
  )
  expect_equal(stats::coef(fit), stats::coef(ref))
  expect_equal(stats::vcov(fit), stats::vcov(ref))
  expect_error(
    suppressWarnings(suppressMessages(
      panhurdle(exacerbations ~ age + smoker, copd, c("id", "visit"))
    )),
    "no time-varying regressors left"
  )
})

test_that("a regressor varying only where a part has no information is dropped from that part", {
  d <- hurdle_panel()
  pos_n <- tapply(d$y > 0 & d$y != 1, d$id, sum)
  # z varies only in individuals with no count-part information
  quiet_ids <- as.integer(names(pos_n)[pos_n == 0])
  expect_gt(length(quiet_ids), 0)
  d$z <- ifelse(d$id %in% quiet_ids, d$time, 0)
  fit <- NULL
  expect_warning(fit <- suppressMessages(panhurdle(y ~ x1 + z, d, c("id", "time"))),
                 "not identified in the count part dropped from it: z")
  expect_true("zero:z" %in% names(stats::coef(fit)))
  expect_false("count:z" %in% names(stats::coef(fit)))
})
