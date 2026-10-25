fmm_sim <- function(N = 200, TT = 3, seed = 1) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
  female <- rbinom(N, 1, 0.5)
  d$female <- rep(female, each = TT)
  high <- rep(rbinom(N, 1, plogis(-1 + female)), each = TT)
  d$x <- rnorm(N * TT)
  eta <- ifelse(high == 1, 1.5 + 0.5 * d$x, -0.3 + 0.2 * d$x)
  d$y <- rpois(N * TT, exp(eta))
  d$ynb <- rnbinom(N * TT, size = 2, mu = exp(eta))
  d$yb <- as.numeric(d$y > 1)
  d$yg <- eta + rnorm(N * TT, sd = 0.5)
  d
}

test_that("k = 1 equals the pooled GLM", {
  d <- fmm_sim()
  a <- panfmm(y ~ x, d, c("id", "t"), k = 1, family = "poisson")
  g <- glm(y ~ x, poisson, d, control = glm.control(epsilon = 1e-14))
  expect_equal(unname(coef(a)), unname(coef(g)), tolerance = 1e-8)
  expect_equal(a$loglik, as.numeric(logLik(g)), tolerance = 1e-10)
  expect_equal(unname(sqrt(diag(vcov(a)))), unname(sqrt(diag(vcov(g)))), tolerance = 1e-5)
  a <- panfmm(ynb ~ x, d, c("id", "t"), k = 1, family = "negbin")
  p <- panglm(ynb ~ x, d, c("id", "t"), model = "pooling", family = "negbin")
  expect_equal(unname(a$classcoef[, 1]), unname(coef(p)), tolerance = 1e-6)
  expect_equal(unname(a$theta), p$theta, tolerance = 1e-5)
  a <- panfmm(yg ~ x, d, c("id", "t"), k = 1, family = "gaussian")
  g <- glm(yg ~ x, gaussian, d)
  expect_equal(a$loglik, as.numeric(logLik(g)), tolerance = 1e-10)
})

test_that("panel Poisson mixture matches flexmix at the same optimum", {
  skip_if_missing("flexmix")
  d <- fmm_sim()
  a <- panfmm(y ~ x, d, c("id", "t"), k = 2, nstart = 3, seed = 1)
  b <- suppressWarnings(flexmix::flexmix(
    y ~ x | id, data = d, k = 2, cluster = a$posterior[as.character(d$id), ],
    model = flexmix::FLXMRglm(family = "poisson"),
    control = list(tolerance = 1e-14, iter.max = 5000)))
  expect_equal(a$loglik, b@logLik, tolerance = 1e-8)
  expect_equal(as.numeric(a$classcoef),
               as.numeric(flexmix::parameters(b)), tolerance = 1e-5)
  expect_equal(unname(a$prior_mean), unname(flexmix::prior(b)), tolerance = 1e-6)
  expect_equal(BIC(a), -2 * b@logLik + b@df * log(nrow(d)), tolerance = 1e-8)
  # classes ordered by ascending mean outcome
  expect_true(a$classmean[1] < a$classmean[2])
  expect_true(all(abs(rowSums(a$posterior) - 1) < 1e-12))
})

test_that("concomitant variables give a multinomial logit prior", {
  skip_if_missing("flexmix")
  d <- fmm_sim()
  a <- panfmm(y ~ x, d, c("id", "t"), k = 2, concomitant = ~ female, nstart = 3, seed = 1)
  b <- suppressWarnings(flexmix::flexmix(
    y ~ x | id, data = d, k = 2, cluster = a$posterior[as.character(d$id), ],
    model = flexmix::FLXMRglm(family = "poisson"),
    concomitant = flexmix::FLXPmultinom(~ female),
    control = list(tolerance = 1e-14, iter.max = 5000)))
  expect_equal(a$loglik, b@logLik, tolerance = 1e-7)
  # flexmix's concomitant M step (nnet::multinom, reltol 1e-8) stops early on
  # this flat likelihood; panfmm is at the optimum (zero score)
  expect_equal(as.numeric(a$prior_coef[, 2]),
               as.numeric(flexmix::parameters(b, which = "concomitant")[, 2]), tolerance = 1e-3)
  expect_lt(max(abs(colSums(a$score))), 1e-5)
  expect_identical(names(coef(a))[5:6], c("prior2:(Intercept)", "prior2:female"))
  expect_error(panfmm(y ~ x, d, c("id", "t"), k = 2, concomitant = ~ x),
               "constant within individuals")
})

test_that("NB2, logit and Gaussian mixtures fit and the score vanishes", {
  d <- fmm_sim()
  for (fam in c("negbin", "binomial", "gaussian")) {
    y <- switch(fam, negbin = "ynb", binomial = "yb", gaussian = "yg")
    f <- stats::as.formula(paste(y, "~ x"))
    a <- panfmm(f, d, c("id", "t"), k = 2, family = fam, nstart = 2, seed = 2)
    expect_s3_class(a, "panfmm")
    expect_lt(max(abs(colSums(a$score))), 1e-4)
    expect_equal(dim(vcov(a)), rep(length(coef(a)), 2))
  }
})

test_that("methods: predict, summary, cluster vcov, seed, panfmmselect", {
  d <- fmm_sim(N = 120)
  a <- panfmm(y ~ x, d, c("id", "t"), k = 2, nstart = 2, seed = 3)
  b <- panfmm(y ~ x, d, c("id", "t"), k = 2, nstart = 2, seed = 3)
  expect_identical(coef(a), coef(b))
  pc <- predict(a, type = "class")
  expect_equal(dim(pc), c(nrow(d), 2L))
  expect_equal(predict(a), as.numeric(pc %*% a$prior_mean), tolerance = 1e-12)
  expect_length(predict(a, weights = "posterior"), nrow(d))
  nd <- data.frame(x = c(0, 1))
  expect_equal(predict(a, newdata = nd, type = "class")[, 1],
               exp(a$classcoef[1, 1] + a$classcoef[2, 1] * nd$x), tolerance = 1e-12)
  expect_error(predict(a, newdata = nd, type = "posterior"), "estimation sample")
  vc <- vcov(a, type = "cluster")
  expect_true(all(is.finite(diag(vc))))
  expect_identical(nobs(a), nrow(d))
  expect_equal(dim(confint(a)), c(length(coef(a)), 2L))
  expect_output(print(summary(a)), "Class 2")
  expect_output(print(a), "latent class Poisson")
  sel <- panfmmselect(y ~ x, d, c("id", "t"), k = 1:2, nstart = 2, seed = 3)
  expect_equal(sel$BIC[2], BIC(a), tolerance = 1e-8)
  expect_true(sel$BIC[2] < sel$BIC[1])
})

test_that("individual weights: weighted GLM at k = 1, replicated individuals at k = 2", {
  d <- fmm_sim()
  d$sw <- rep(runif(200, 0.5, 2), each = 3)
  a <- panfmm(y ~ x, d, c("id", "t"), k = 1, weights = d$sw)
  g <- suppressWarnings(glm(y ~ x, poisson, d, weights = sw,
                            control = glm.control(epsilon = 1e-14)))
  expect_equal(unname(coef(a)), unname(coef(g)), tolerance = 1e-8)
  expect_equal(a$loglik, as.numeric(logLik(g)), tolerance = 1e-10)
  expect_equal(unname(sqrt(diag(vcov(a)))), unname(sqrt(diag(vcov(g)))), tolerance = 1e-5)
  skip_if_missing("sandwich")
  expect_equal(unname(sqrt(diag(vcov(a, type = "cluster")))),
               unname(sqrt(diag(sandwich::vcovCL(g, cluster = ~ id, type = "HC0")))),
               tolerance = 1e-6)

  d$w <- rep(rep(1:2, 100), each = 3)
  a <- panfmm(y ~ x, d, c("id", "t"), k = 2, weights = d$w, nstart = 3, seed = 1)
  dr <- d[rep(seq_len(nrow(d)), d$w), ]
  dr$rid <- paste(dr$id, ave(dr$t, dr$id, dr$t, FUN = seq_along))
  b <- panfmm(y ~ x, dr, c("rid", "t"), k = 2, nstart = 3, seed = 1)
  expect_equal(coef(a), coef(b), tolerance = 1e-7)
  expect_equal(a$loglik, b$loglik, tolerance = 1e-10)
  expect_equal(sqrt(diag(vcov(a))), sqrt(diag(vcov(b))), tolerance = 1e-6)
  expect_lt(max(abs(colSums(a$weights * a$score))), 1e-5)
  expect_length(a$weights, 200L)
  expect_null(panfmm(y ~ x, d, c("id", "t"), k = 1, weights = rep(1, nrow(d)))$weights)

  expect_error(panfmm(y ~ x, d, c("id", "t"), k = 2, weights = runif(nrow(d))),
               "constant within individuals")
  expect_error(panfmm(y ~ x, d, c("id", "t"), k = 2, weights = -d$w), "nonnegative")
  expect_error(panfmm(y ~ x, d, c("id", "t"), k = 2, weights = d$w[-1]), "nrow")
})
