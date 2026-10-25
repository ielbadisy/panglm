test_that("pooled negbin (real NB2, not aliased to Poisson) matches MASS::glm.nb exactly", {
  skip_if_missing("plm")
  skip_if_missing("MASS")
  data(Grunfeld, package = "plm")
  set.seed(1)
  Grunfeld$count <- rnbinom(nrow(Grunfeld), size = 3, mu = exp(0.5 + 0.0002 * Grunfeld$value))

  f <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "pooling", family = "negbin")
  g <- MASS::glm.nb(count ~ value + capital, data = Grunfeld)

  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-4)
  expect_equal(f$theta, g$theta, tolerance = 1e-3)
  expect_equal(f$loglik, as.numeric(stats::logLik(g)), tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(vcov(g)))), tolerance = 1e-3)
})

test_that("pooled negbin is NOT the same fit as pooled poisson (regression guard for the aliasing bug)", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  set.seed(1)
  Grunfeld$count <- rnbinom(nrow(Grunfeld), size = 3, mu = exp(0.5 + 0.0002 * Grunfeld$value))

  f_nb <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "pooling", family = "negbin")
  f_pois <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "pooling", family = "poisson")
  expect_false(isTRUE(all.equal(unname(coef(f_nb)), unname(coef(f_pois)))))
  expect_false(isTRUE(all.equal(f_nb$loglik, f_pois$loglik)))
})

test_that("within negbin (Allison-Waterman dummy FE) matches fixest::fenegbin exactly", {
  skip_if_missing("plm")
  skip_if_missing("fixest")
  data(Grunfeld, package = "plm")
  set.seed(2)
  u <- rep(rnorm(10, sd = 0.5), each = 20)
  Grunfeld$count <- rnbinom(nrow(Grunfeld), size = 4, mu = exp(1 + 0.0002 * Grunfeld$value + u))

  f <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "within", family = "negbin")
  fx <- fixest::fenegbin(count ~ value + capital | firm, data = Grunfeld, notes = FALSE)

  expect_equal(unname(coef(f)), unname(coef(fx)), tolerance = 1e-4)
  expect_equal(f$loglik, as.numeric(stats::logLik(fx)), tolerance = 1e-3)
  expect_length(f$individual_effects, 10)
})

test_that("within negbin is NOT the same fit as within poisson (regression guard for the aliasing bug)", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  set.seed(2)
  u <- rep(rnorm(10, sd = 0.5), each = 20)
  Grunfeld$count <- rnbinom(nrow(Grunfeld), size = 4, mu = exp(1 + 0.0002 * Grunfeld$value + u))

  f_nb <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "within", family = "negbin")
  f_pois <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "within", family = "poisson")
  expect_false(isTRUE(all.equal(unname(coef(f_nb)), unname(coef(f_pois)))))
})

test_that("random negbin (2-parameter beta-negative-binomial) recovers known simulated parameters", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  set.seed(3)
  Grunfeld$u <- rep(rgamma(10, shape = 6, rate = 6), each = 20)
  Grunfeld$count <- rnbinom(nrow(Grunfeld), size = 8, mu = Grunfeld$u * exp(1 + 0.0002 * Grunfeld$value))

  f <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "random", family = "negbin")
  expect_true(is.finite(f$loglik))
  expect_true(f$negbin_a > 0 && f$negbin_b > 0)
  # sanity: coefficients in the right ballpark of the simulated truth
  expect_equal(unname(coef(f)["value"]), 0.0002, tolerance = 1)
})

test_that("random negbin is NOT the same fit as random poisson (regression guard for the aliasing bug)", {
  skip_if_missing("plm")
  data(Grunfeld, package = "plm")
  set.seed(3)
  Grunfeld$u <- rep(rgamma(10, shape = 6, rate = 6), each = 20)
  Grunfeld$count <- rnbinom(nrow(Grunfeld), size = 8, mu = Grunfeld$u * exp(1 + 0.0002 * Grunfeld$value))

  f_nb <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "random", family = "negbin")
  f_pois <- panglm(count ~ value + capital, data = Grunfeld, index = c("firm", "year"), model = "random", family = "poisson")
  expect_false(isTRUE(all.equal(f_nb$loglik, f_pois$loglik)))
})

test_that("random negbin at the a -> Inf boundary warns and reproduces the gamma Poisson slopes", {
  f <- exacerbations ~ treatment + age + smoker + crp
  expect_warning(nb <- panglm(f, copd, index = c("id", "visit"), model = "random",
                              family = "negbin"), "diverges")
  pois <- panglm(f, copd, index = c("id", "visit"), model = "random", family = "poisson")
  expect_gt(nb$negbin_a, 1e4)
  expect_equal(coef(nb)[-1], coef(pois)[-1], tolerance = 1e-3)
  expect_equal(nb$loglik, pois$loglik, tolerance = 1e-5)
  # the slopes keep usable standard errors
  expect_true(all(diag(vcov(nb))[-1] > 0))
  # weighted (generic optimizer) path
  expect_warning(panglm(f, copd, index = c("id", "visit"), model = "random",
                        family = "negbin", weights = rep(2, nrow(copd))), "diverges")
})

test_that("random negbin away from the boundary is silent", {
  set.seed(11)
  G <- 200; Tt <- 5
  d <- data.frame(id = rep(seq_len(G), each = Tt), t = rep(seq_len(Tt), G), x = rnorm(G * Tt))
  d$y <- rnbinom(G * Tt, mu = exp(0.5 + 0.4 * d$x + rep(rnorm(G, 0, 0.5), each = Tt)), size = 1.5)
  m <- expect_silent(panglm(y ~ x, d, index = c("id", "t"), model = "random", family = "negbin"))
  expect_lt(m$negbin_a, 1e4)
  expect_true(all(diag(vcov(m)) > 0))
})

test_that("negbin summaries report theta, not a dispersion of 1", {
  data(copd, package = "panglm", envir = environment())
  pooled <- panglm(exacerbations ~ crp + treatment, data = copd,
                   index = c("id", "visit"), model = "pooling", family = "negbin")
  s <- summary(pooled)
  expect_null(s$dispersion)
  expect_equal(s$negbin_theta, pooled$theta)
  out <- paste(capture.output(s), collapse = "\n")
  expect_no_match(out, "Dispersion parameter", fixed = TRUE)
  expect_match(out, "Negative binomial theta (NB2):", fixed = TRUE)
  expect_match(out, "(SE ", fixed = TRUE)
  expect_equal(generics::glance(pooled)$theta, pooled$theta)
  skip_if_not_installed("MASS")
  ref <- MASS::glm.nb(exacerbations ~ crp + treatment, data = copd)
  expect_equal(s$negbin_theta, ref$theta, tolerance = 1e-5)
  expect_equal(s$negbin_theta_se, ref$SE.theta, tolerance = 1e-4)

  within <- suppressMessages(panglm(exacerbations ~ crp, data = copd,
                                    index = c("id", "visit"), model = "within",
                                    family = "negbin"))
  sw <- summary(within)
  expect_null(sw$dispersion)
  expect_true(is.finite(sw$negbin_theta) && is.finite(sw$negbin_theta_se))
  expect_no_match(paste(capture.output(sw), collapse = "\n"), "Dispersion parameter",
                  fixed = TRUE)
  expect_equal(generics::glance(within)$theta, within$theta)
  expect_true(is.na(generics::glance(panglm(exacerbations ~ crp, data = copd,
    index = c("id", "visit"), model = "pooling", family = "poisson"))$theta))
})

test_that("random-intercept negbin summaries report the beta mixing parameters", {
  data(copd, package = "panglm", envir = environment())
  fit <- suppressWarnings(panglm(exacerbations ~ crp, data = copd, index = c("id", "visit"),
                                 model = "random", family = "negbin"))
  s <- suppressWarnings(summary(fit))
  expect_equal(s$negbin_ab, c(a = fit$negbin_a, b = fit$negbin_b))
  expect_match(paste(suppressWarnings(capture.output(s)), collapse = "\n"), "Beta mixing parameters", fixed = TRUE)
  expect_true(is.na(suppressWarnings(generics::glance(fit))$theta))
})
