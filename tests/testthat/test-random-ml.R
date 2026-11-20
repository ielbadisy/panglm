test_that("Gaussian random intercept by ML and REML matches lme4::lmer()", {
  skip_if_missing("lme4")
  data(copd)
  d <- copd[!(copd$id %% 7 == 0 & copd$visit > 2), ] # unbalanced
  for (reml in c(FALSE, TRUE)) {
    fit <- panglm(fev1 ~ visit + treatment + age, d, index = c("id", "visit"),
                  model = "random", method = if (reml) "reml" else "ml")
    ref <- lme4::lmer(fev1 ~ visit + treatment + age + (1 | id), d, REML = reml)
    vc <- as.data.frame(lme4::VarCorr(ref))$vcov
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-7)
    expect_equal(attr(logLik(fit), "df"), attr(logLik(ref), "df"))
    expect_equal(coef(fit), lme4::fixef(ref), tolerance = 1e-6)
    expect_equal(unname(sqrt(diag(vcov(fit)))), unname(sqrt(diag(as.matrix(vcov(ref))))),
                 tolerance = 1e-5)
    expect_equal(c(fit$sigma_mu2, fit$sigma_v2), vc, tolerance = 1e-5)
    expect_identical(fit$method, if (reml) "reml" else "ml")
  }
})

test_that("ML random intercept is nested in the random-slopes fit (LR test)", {
  skip_if_missing("lme4")
  data(copd)
  ml <- panglm(fev1 ~ visit + treatment, copd, index = c("id", "visit"), model = "random",
               method = "ml")
  # the slope variance is near zero in copd: both fits warn of a singular fit
  rs <- suppressWarnings(panglm(fev1 ~ visit + treatment, copd, index = c("id", "visit"),
                                model = "random", random = ~ visit))
  r0 <- lme4::lmer(fev1 ~ visit + treatment + (1 | id), copd, REML = FALSE)
  r1 <- suppressMessages(lme4::lmer(fev1 ~ visit + treatment + (visit | id), copd, REML = FALSE))
  expect_equal(as.numeric(logLik(rs) - logLik(ml)), as.numeric(logLik(r1) - logLik(r0)),
               tolerance = 1e-4)
  expect_identical(rs$method, "ml")
})

test_that("a zero random-intercept variance is found on the boundary", {
  set.seed(11)
  d <- data.frame(id = rep(1:30, each = 4), t = rep(1:4, 30))
  d$x <- rnorm(120); d$y <- 1 + d$x + rnorm(120)
  ml <- panglm(y ~ x, d, index = c("id", "t"), model = "random", method = "ml")
  pool <- panglm(y ~ x, d, index = c("id", "t"), model = "pooling")
  if (ml$sigma_mu2 < 1e-8) expect_equal(as.numeric(logLik(ml)), as.numeric(logLik(pool)))
  expect_gte(as.numeric(logLik(ml)), as.numeric(logLik(pool)) - 1e-8)
})

test_that("method keeps the Swamy-Arora default and rejects other models", {
  data(copd)
  a <- panglm(fev1 ~ crp, copd, index = c("id", "visit"), model = "random")
  b <- panglm(fev1 ~ crp, copd, index = c("id", "visit"), model = "random", method = "swar")
  expect_identical(coef(a), coef(b)); expect_identical(a$method, "swar")
  expect_error(panglm(fev1 ~ crp, copd, index = c("id", "visit"), method = "ml"),
               "applies to model = 'random' with family = 'gaussian'")
  expect_error(panglm(exacerbations ~ crp, copd, index = c("id", "visit"), model = "random",
                      family = "poisson", method = "reml"), "family = 'gaussian'")
  expect_error(panglm(fev1 ~ visit, copd, index = c("id", "visit"), model = "random",
                      random = ~ visit, method = "reml"), "unweighted random intercept")
})

test_that("cluster covariance of the ML fit is the GLS cluster sandwich", {
  data(copd)
  fit <- panglm(fev1 ~ crp + treatment, copd, index = c("id", "visit"), model = "random",
                method = "ml", vcov = "cluster")
  qd <- quasi_demean_cpp(fit$X, fit$y, fit$group_start, fit$group_size, fit$theta)
  B <- solve(crossprod(qd$X)); e <- as.numeric(qd$y - qd$X %*% coef(fit))
  S <- rowsum(e * qd$X, rep(seq_along(fit$group_size), fit$group_size))
  G <- nrow(S)
  expect_equal(unname(vcov(fit)), unname(B %*% crossprod(S) %*% B * G / (G - 1)),
               tolerance = 1e-8)
})
