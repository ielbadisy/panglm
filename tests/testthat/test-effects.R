group_of <- function(fit) rep(seq_along(fit$group_size), fit$group_size)

# Posterior mean of a scalar random intercept by brute-force integration.
brute_mean <- function(fit, i, logf_i, sigma) {
  g <- group_of(fit)
  eta <- fit$linear.predictors[g == i]
  yy <- fit$y[g == i]
  f <- function(u) vapply(u, function(v) exp(logf_i(yy, eta + v)) * stats::dnorm(v, 0, sigma), 0)
  stats::integrate(function(u) u * f(u), -Inf, Inf, rel.tol = 1e-12)$value /
    stats::integrate(f, -Inf, Inf, rel.tol = 1e-12)$value
}

test_that("Gaussian within effects match fixest::fixef()", {
  skip_if_missing("fixest")
  data(copd)
  fit <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"), model = "within")
  ref <- fixest::fixef(fixest::feols(fev1 ~ crp | id, data = copd))$id
  a <- fixef(fit)
  expect_setequal(names(a), names(ref))
  expect_equal(unname(a), unname(ref[names(a)]), tolerance = 1e-8)
})

test_that("Poisson within effects match fixest::fepois()", {
  skip_if_missing("fixest")
  data(copd)
  fit <- panglm(exacerbations ~ crp, data = copd, index = c("id", "visit"),
                model = "within", family = "poisson")
  ref <- fixest::fixef(fixest::fepois(exacerbations ~ crp | id, data = copd))$id
  a <- fixef(fit)
  expect_setequal(names(a), names(ref))
  expect_equal(unname(a), unname(ref[names(a)]), tolerance = 1e-6)
  # closed form
  g <- group_of(fit)
  xb <- drop(fit$X %*% coef(fit))
  cf <- log(rowsum(fit$y, g)[, 1] / rowsum(exp(xb), g)[, 1])
  names(cf) <- fit$group_levels
  expect_equal(a, cf[is.finite(cf)], tolerance = 1e-10)
})

test_that("fixef() refuses fits that do not estimate individual effects", {
  data(copd)
  cl <- panglm(hospitalized ~ crp, data = copd, index = c("id", "visit"),
               model = "within", family = "binomial")
  expect_error(fixef(cl), "conditional logit")
  po <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"), model = "pooling")
  expect_error(fixef(po), "within")
  expect_error(ranef(po), "random")
  re <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"), model = "random")
  expect_error(fixef(re), "ranef")
})

test_that("two-way fits return individual and time effects", {
  data(copd)
  fit <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"), model = "within",
                effect = "twoways")
  fe <- fixef(fit)
  expect_named(fe, c("individual", "time"))
  expect_equal(mean(fe$time), 0, tolerance = 1e-10)
})

test_that("Gaussian random-intercept BLUPs match lme4 with plugged components", {
  skip_if_missing("lme4")
  data(copd)
  fit <- panglm(fev1 ~ crp + treatment, data = copd, index = c("id", "visit"),
                model = "random")
  r <- ranef(fit)
  expect_identical(dim(r), c(fit$n_groups, 1L))
  expect_identical(rownames(r), as.character(fit$group_levels))
  expect_equal(ranef(fit, type = "mode"), structure(r, type = "mode"))
  m <- lme4::lmer(fev1 ~ crp + treatment + (1 | id), data = copd, REML = FALSE)
  ref <- lme4::ranef(m)$id
  ref <- stats::setNames(ref[, 1], rownames(ref))[rownames(r)]
  # Swamy-Arora vs ML components: close but not identical
  expect_gt(stats::cor(r[, 1], ref), 0.9999)
  expect_lt(max(abs(r[, 1] - ref)), 1e-2)
  # with lme4's beta and variance components plugged in: identical
  vc <- as.data.frame(lme4::VarCorr(m))$vcov
  plug <- fit
  plug$linear.predictors <- drop(fit$X %*% lme4::fixef(m)[colnames(fit$X)])
  plug$sigma_mu2 <- vc[1]
  plug$sigma_v2 <- vc[2]
  expect_equal(ranef(plug)[, 1], ref, tolerance = 1e-8)
})

test_that("random-intercept logit posterior means agree with brute force, glmer and GLMMadaptive", {
  data(copd)
  fit <- panglm(hospitalized ~ crp + treatment, data = copd, index = c("id", "visit"),
                model = "random", family = "binomial")
  r <- ranef(fit)
  lf <- function(y, eta) sum(stats::dbinom(y, 1, stats::plogis(eta), log = TRUE))
  bf <- vapply(1:4, function(i) brute_mean(fit, i, lf, fit$sigma), 0)
  expect_equal(unname(r[1:4, 1]), bf, tolerance = 1e-7)
  skip_if_missing("lme4")
  m <- lme4::glmer(hospitalized ~ crp + treatment + (1 | id), data = copd,
                   family = stats::binomial, nAGQ = 21)
  ref <- lme4::ranef(m)$id
  ref <- stats::setNames(ref[, 1], rownames(ref))[rownames(r)]
  # glmer reports conditional modes
  rm <- ranef(fit, type = "mode")
  expect_gt(stats::cor(rm[, 1], ref), 0.9999)
  expect_lt(max(abs(rm[, 1] - ref)), 1e-3)
  expect_gt(stats::cor(r[, 1], ref), 0.999)
  skip_if_missing("GLMMadaptive")
  ga <- GLMMadaptive::mixed_model(hospitalized ~ crp + treatment, random = ~ 1 | id,
                                  data = copd, family = stats::binomial(), nAGQ = 21)
  ra <- GLMMadaptive::ranef(ga, type = "mean")
  expect_gt(stats::cor(r[, 1], ra[rownames(r), 1]), 0.999)
})

test_that("posterior means of other normal-mixing families agree with brute force", {
  data(copd)
  ix <- c("id", "visit")
  gm <- panglm(fev1 ~ crp, data = copd, index = ix, model = "random", family = "gamma")
  lf <- function(y, eta) sum(stats::dgamma(y, shape = gm$shape, rate = gm$shape / exp(eta), log = TRUE))
  expect_equal(unname(ranef(gm)[1:3, 1]),
               vapply(1:3, function(i) brute_mean(gm, i, lf, gm$sigma), 0), tolerance = 1e-7)

  copd$y01 <- (copd$fev1 - min(copd$fev1) + 0.05) / (diff(range(copd$fev1)) + 0.1)
  bt <- panglm(y01 ~ crp, data = copd, index = ix, model = "random", family = "beta")
  lf <- function(y, eta) {
    mu <- stats::plogis(eta)
    sum(stats::dbeta(y, mu * bt$precision, (1 - mu) * bt$precision, log = TRUE))
  }
  expect_equal(unname(ranef(bt)[1:3, 1]),
               vapply(1:3, function(i) brute_mean(bt, i, lf, bt$sigma), 0), tolerance = 1e-7)
})

test_that("Poisson-gamma random effects have the closed-form posterior", {
  data(copd)
  fit <- panglm(exacerbations ~ crp, data = copd, index = c("id", "visit"),
                model = "random", family = "poisson")
  g <- group_of(fit)
  A <- fit$dispersion_param + rowsum(fit$y, g)[, 1]
  B <- fit$dispersion_param + rowsum(exp(fit$linear.predictors), g)[, 1]
  expect_equal(unname(ranef(fit)[, 1]), unname(digamma(A) - log(B)), tolerance = 1e-12)
  expect_equal(unname(ranef(fit, type = "mode")[, 1]), unname(log(A / B)), tolerance = 1e-12)
})

test_that("panglm's generics dispatch through nlme's", {
  skip_if_missing("nlme")
  data(copd)
  fit <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"), model = "within")
  expect_equal(nlme::fixef(fit), fixef(fit))
})
