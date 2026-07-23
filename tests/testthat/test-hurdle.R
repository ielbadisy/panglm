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

test_that("panglm_hurdle() fits a panel with structural zeros and recovers the count-part slope", {
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

  fit <- panglm_hurdle(y ~ x1, data = d, index = c("id", "time"))
  expect_s3_class(fit, "panglm_hurdle")
  expect_true(is.finite(fit$count$coefficients["x1"]))
  expect_equal(unname(fit$count$coefficients["x1"]), 0.5, tolerance = 0.3)
  expect_true(fit$n_positive < fit$n)
})
