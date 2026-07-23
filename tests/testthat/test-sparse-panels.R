# These DGPs are deliberately adversarial (heavy zero-inflation, extreme
# skew across units, unbalanced group sizes) -- the kind of real-world
# count panel that first surfaced a singular-Hessian crash in the
# random-effects Poisson-Gamma fit. Keep this file even when it looks
# redundant with test-random.R's "well-behaved" DGPs: it's the permanent
# regression guard for that failure mode, not an ad hoc one-off.

test_that("random poisson (Poisson-Gamma) does not crash on a heavily zero-inflated, skewed panel", {
  set.seed(99)
  n_id <- 40; n_t <- 5
  d <- data.frame(id = rep(1:n_id, each = n_t), time = rep(1:n_t, n_id))
  d$x1 <- rnorm(nrow(d))
  # 35 near-dead units, 5 highly active units -- extreme cross-unit skew
  u <- rep(c(rep(-5, 35), rep(2, 5)), each = n_t)
  d$y <- rpois(nrow(d), lambda = exp(-1 + 0.3 * d$x1 + u))

  expect_true(mean(d$y == 0) > 0.7) # confirm the DGP really is heavily zero-inflated
  fit <- suppressWarnings(panglm(y ~ x1, data = d, index = c("id", "time"), model = "random", family = "poisson"))
  expect_true(is.finite(fit$loglik))
  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(diag(vcov(fit)))))
})

test_that("random poisson on a sparse panel warns (not errors) when LM damping engages", {
  set.seed(99)
  n_id <- 40; n_t <- 5
  d <- data.frame(id = rep(1:n_id, each = n_t), time = rep(1:n_t, n_id))
  d$x1 <- rnorm(nrow(d))
  u <- rep(c(rep(-5, 35), rep(2, 5)), each = n_t)
  d$y <- rpois(nrow(d), lambda = exp(-1 + 0.3 * d$x1 + u))

  # Either it converges without damping (no warning) or it warns explicitly
  # about damping -- it must never throw an uncaught error.
  expect_no_error(
    withCallingHandlers(
      panglm(y ~ x1, data = d, index = c("id", "time"), model = "random", family = "poisson"),
      warning = function(w) invokeRestart("muffleWarning")
    )
  )
})

test_that("within poisson handles unbalanced group sizes without error", {
  set.seed(7)
  n_id <- 30
  sizes <- sample(2:8, n_id, replace = TRUE)
  ids <- rep(seq_len(n_id), sizes)
  d <- data.frame(id = ids, time = ave(ids, ids, FUN = seq_along))
  d$x1 <- rnorm(nrow(d))
  a <- rep(rnorm(n_id), sizes)
  d$y <- rpois(nrow(d), lambda = exp(0.2 + 0.4 * d$x1 + a))

  fit <- panglm(y ~ x1, data = d, index = c("id", "time"), model = "within", family = "poisson")
  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(fitted(fit))))
})

test_that("within negbin (Allison-Waterman) handles a sparse, unbalanced count panel", {
  set.seed(11)
  n_id <- 25
  sizes <- sample(3:6, n_id, replace = TRUE)
  ids <- rep(seq_len(n_id), sizes)
  d <- data.frame(id = ids, time = ave(ids, ids, FUN = seq_along))
  d$x1 <- rnorm(nrow(d))
  a <- rep(rnorm(n_id, sd = 0.7), sizes)
  d$y <- rnbinom(nrow(d), size = 3, mu = exp(-0.5 + 0.3 * d$x1 + a))

  fit <- panglm(y ~ x1, data = d, index = c("id", "time"), model = "within", family = "negbin")
  expect_true(is.finite(fit$loglik))
  expect_true(fit$theta > 0)
  # All-zero groups are dropped before fitting (their fixed effect is
  # unbounded under the dummy-variable NB2 MLE); individual_effects should
  # cover exactly the retained groups.
  expect_length(fit$individual_effects, fit$n_used_groups)
  expect_equal(fit$n_used_groups + fit$n_dropped_groups, n_id)
})
