slopes_sim <- function(G = 80, T = 5, seed = 7) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(G), each = T), time = rep(0:(T - 1), G))
  d$x <- stats::rnorm(G * T)
  b1 <- stats::rnorm(G, 0, 0.5)
  b2 <- -0.2 * b1 + stats::rnorm(G, 0, 0.15)
  lin <- 0.1 * d$time + 0.4 * d$x + b1[d$id] + b2[d$id] * d$time
  d$y <- 1 + lin + stats::rnorm(G * T, 0, 0.3)
  d$yc <- stats::rpois(G * T, exp(0.2 + lin))
  d$yb <- stats::rbinom(G * T, 1, stats::plogis(-0.3 + lin))
  d
}

test_that("random = ~ 1 leaves the random-intercept fit unchanged", {
  data(copd)
  a <- panglm(fev1 ~ crp + age, data = copd, index = c("id", "visit"), model = "random")
  b <- panglm(fev1 ~ crp + age, data = copd, index = c("id", "visit"), model = "random",
              random = ~ 1)
  a$call <- b$call <- NULL
  expect_identical(a, b)
  expect_null(a$varcomp)
})

test_that("Gaussian random slopes match lme4::lmer(REML = FALSE)", {
  skip_if_missing("lme4")
  d <- slopes_sim()
  f <- panglm(y ~ time + x, data = d, index = c("id", "time"), model = "random",
              random = ~ time)
  m <- lme4::lmer(y ~ time + x + (1 + time | id), data = d, REML = FALSE)
  expect_equal(unname(coef(f)), unname(lme4::fixef(m)), tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(as.matrix(vcov(m))))),
               tolerance = 1e-5)
  vc <- lme4::VarCorr(m)$id
  expect_equal(unname(f$varcomp), unname(vc[, ]), tolerance = 1e-4)
  expect_equal(f$sigma_v2, stats::sigma(m)^2, tolerance = 1e-5)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(m)), tolerance = 1e-8)
  expect_equal(attr(logLik(f), "df"), attr(logLik(m), "df"))
  expect_equal(unname(f$ranef), unname(as.matrix(lme4::ranef(m)$id)), tolerance = 1e-4)
  expect_identical(rownames(f$ranef), as.character(sort(unique(d$id))))
})

test_that("Gaussian random slopes: two slopes, cre, cluster vcov, methods", {
  d <- slopes_sim()
  d$x2 <- stats::rnorm(nrow(d))
  f2 <- panglm(y ~ time + x, data = d, index = c("id", "time"), model = "random",
               random = ~ time + x)
  expect_equal(dim(f2$varcomp), c(3L, 3L))
  expect_equal(dim(f2$ranef), c(80L, 3L))

  fc <- panglm(y ~ x + x2, data = d, index = c("id", "time"), model = "random",
               random = ~ time, cre = TRUE, vcov = "cluster")
  expect_true(all(c("mean_x", "mean_x2") %in% names(coef(fc))))
  expect_identical(fc$vcov_type, "cluster")
  expect_false(isTRUE(all.equal(vcov(fc), vcov(fc, type = "classical"))))
  expect_error(vcov(fc, type = "HC1"), "HC1")
  expect_output(print(summary(fc)), "Random-effects correlation")
  expect_output(print(fc), "panglm")

  f <- panglm(y ~ time + x, data = d, index = c("id", "time"), model = "random",
              random = ~ time)
  nd <- d[1:6, ]
  expect_equal(predict(f, newdata = nd, type = "link"),
               as.numeric(cbind(1, nd$time, nd$x) %*% coef(f)))
  expect_equal(predict(f, newdata = nd, type = "response"),
               predict(f, newdata = nd, type = "link"))
  ape <- partialeffects(f, "x")
  expect_equal(ape$estimate, unname(coef(f)["x"]), tolerance = 1e-8)
})

test_that("Poisson and logit random slopes (Laplace) match glmmTMB", {
  skip_on_cran()
  skip_if_missing("glmmTMB")
  d <- slopes_sim()
  for (spec in list(list("yc", "poisson", stats::poisson()),
                    list("yb", "binomial", stats::binomial()))) {
    fm <- stats::as.formula(paste(spec[[1]], "~ time + x"))
    f <- panglm(fm, data = d, index = c("id", "time"), model = "random",
                family = spec[[2]], random = ~ time, R = 1)
    tm <- glmmTMB::glmmTMB(stats::update(fm, . ~ . + (1 + time | id)), data = d,
                           family = spec[[3]])
    expect_equal(unname(coef(f)), unname(glmmTMB::fixef(tm)$cond), tolerance = 1e-4)
    expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(vcov(tm)$cond))),
                 tolerance = 1e-4)
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(tm)), tolerance = 1e-6)
    expect_equal(unname(f$varcomp), unname(glmmTMB::VarCorr(tm)$cond$id[, ]),
                 tolerance = 1e-3)
  }
})

test_that("GLM random slopes with two slopes match glmmTMB (Laplace)", {
  skip_on_cran()
  skip_if_missing("glmmTMB")
  set.seed(2)
  G <- 120; T <- 8
  d <- data.frame(id = rep(seq_len(G), each = T), time = rep(0:(T - 1), G))
  d$x <- stats::rnorm(G * T); d$w <- stats::rnorm(G * T)
  L <- matrix(c(0.8, 0, 0, -0.1, 0.25, 0, 0.2, 0.1, 0.5), 3, 3, byrow = TRUE)
  b <- matrix(stats::rnorm(G * 3), G) %*% t(L)
  lin <- 0.1 * d$time + 0.4 * d$x + b[d$id, 1] + b[d$id, 2] * d$time + b[d$id, 3] * d$w
  d$yc <- stats::rpois(G * T, exp(0.2 + lin))
  f <- panglm(yc ~ time + x + w, data = d, index = c("id", "time"), model = "random",
              family = "poisson", random = ~ time + w, R = 1)
  tm <- glmmTMB::glmmTMB(yc ~ time + x + w + (1 + time + w | id), data = d,
                         family = stats::poisson())
  expect_equal(dim(f$varcomp), c(3L, 3L))
  expect_equal(unname(coef(f)), unname(glmmTMB::fixef(tm)$cond), tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(vcov(tm)$cond))),
               tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(tm)), tolerance = 1e-6)
  expect_equal(unname(f$varcomp), unname(glmmTMB::VarCorr(tm)$cond$id[, ]),
               tolerance = 1e-3)
  expect_equal(attr(logLik(f), "df"), 10)
})

test_that("a GLM slope variance at the boundary keeps finite standard errors", {
  skip_on_cran()
  skip_if_missing("glmmTMB")
  set.seed(4)
  G <- 200; T <- 6
  d <- data.frame(id = rep(seq_len(G), each = T), time = rep(0:(T - 1), G))
  d$x <- stats::rnorm(G * T)
  u <- stats::rnorm(G, 0, 0.6)
  d$yc <- stats::rpois(G * T, exp(0.3 + 0.1 * d$time + 0.4 * d$x + u[d$id]))
  expect_warning(
    f <- panglm(yc ~ time + x, data = d, index = c("id", "time"), model = "random",
                family = "poisson", random = ~ time, R = 1),
    "singular")
  se <- sqrt(diag(vcov(f)))
  expect_true(all(is.finite(se)))
  tm <- suppressWarnings(glmmTMB::glmmTMB(yc ~ time + x + (1 + time | id), data = d,
                                          family = stats::poisson()))
  expect_equal(unname(se), unname(sqrt(diag(vcov(tm)$cond))), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(tm)), tolerance = 1e-6)
  # a singular information matrix falls back on its generalized inverse
  H <- -diag(c(4, 1, 0))
  V <- slopes_theta_vcov(H)
  expect_true(all(is.finite(V)))
  expect_equal(V, diag(c(0.25, 1, 0)))
})

test_that("analytic scores of the quadrature likelihood match numerical ones", {
  set.seed(5)
  G <- 25; T <- 5
  d <- data.frame(id = rep(seq_len(G), each = T), time = rep(0:(T - 1), G))
  d$x <- stats::rnorm(G * T); d$w <- stats::rnorm(G * T)
  lin <- 0.3 * d$x + rep(stats::rnorm(G, 0, 0.6), each = T) +
    rep(stats::rnorm(G, 0, 0.2), each = T) * d$time
  ys <- list(poisson = stats::rpois(G * T, exp(0.2 + lin)),
             logit = stats::rbinom(G * T, 1, stats::plogis(lin)),
             probit = stats::rbinom(G * T, 1, stats::pnorm(lin)))
  fams <- list(poisson = stats::poisson(), logit = stats::binomial(),
               probit = stats::binomial("probit"))
  X <- cbind(1, d$time, d$x)
  for (nm in names(ys)) for (q in 2:3) for (R in c(1L, 3L)) {
    Z <- cbind(1, d$time, d$w)[, seq_len(q), drop = FALSE]
    contrib <- slopes_glm_contrib(fams[[nm]], X, ys[[nm]], Z, d$id, R)
    theta <- c(0.1, 0.05, 0.3, if (q == 2) c(0.6, 0.05, 0.2) else c(0.6, 0.05, 0.1, 0.2, 0.05, 0.3))
    contrib(theta, FALSE)
    expect_equal(contrib(theta, TRUE)$score, ml_numeric_scores(contrib, theta),
                 tolerance = 1e-6, info = paste(nm, q, R))
  }
})

test_that("the default number of quadrature nodes adapts to the dimension", {
  expect_identical(panglm:::slopes_default_nodes(2L), 7L)
  expect_identical(panglm:::slopes_default_nodes(3L), 5L)
  expect_identical(panglm:::slopes_default_nodes(5L), 3L)
  expect_identical(panglm:::slopes_default_nodes(8L), 1L)
  d <- slopes_sim(G = 30)
  d$w <- stats::rnorm(nrow(d))
  # w carries no random effect in these data: a boundary fit
  f <- suppressWarnings(panglm(yc ~ time + x, data = d, index = c("id", "time"),
                               model = "random", family = "poisson", random = ~ time + w))
  expect_identical(f$quadrature_nodes, 5L)
  expect_warning(panglm:::slopes_nodes(13L, 3L), "nodes per individual")
  expect_silent(panglm:::slopes_nodes(7L, 3L))
})

test_that("adaptive quadrature converges and population-averaged methods work", {
  d <- slopes_sim()
  f7 <- panglm(yc ~ time + x, data = d, index = c("id", "time"), model = "random",
               family = "poisson", random = ~ time)
  f11 <- panglm(yc ~ time + x, data = d, index = c("id", "time"), model = "random",
                family = "poisson", random = ~ time, R = 11)
  expect_identical(f7$quadrature_nodes, 7L)
  expect_equal(coef(f7), coef(f11), tolerance = 1e-4)
  expect_equal(f7$loglik, f11$loglik, tolerance = 1e-6)

  Z <- cbind(1, d$time)
  s2 <- rowSums((Z %*% f7$varcomp) * Z)
  eta <- predict(f7, type = "link")
  expect_equal(fitted(f7), exp(eta + s2 / 2))
  ape <- partialeffects(f7, "x")
  expect_equal(ape$estimate, mean(coef(f7)[["x"]] * exp(eta + s2 / 2)), tolerance = 1e-7)
  vc <- vcov(f7, type = "cluster")
  expect_true(all(is.finite(vc)))
  expect_equal(dim(vc), c(3L, 3L))
  expect_equal(attr(logLik(f7), "df"), 6)
})

test_that("random slopes give clear errors outside their scope", {
  d <- slopes_sim(G = 30)
  idx <- c("id", "time")
  expect_error(panglm(y ~ x, data = d, index = idx, model = "pooling", random = ~ time),
               "model = 'random'")
  expect_error(panglm(y ~ x, data = d, index = idx, model = "random", random = y ~ time),
               "one-sided")
  expect_error(panglm(y ~ x, data = d, index = idx, model = "random", random = ~ 0 + time),
               "intercept")
  expect_error(panglm(yc ~ x, data = d, index = idx, model = "random", family = "negbin",
                      random = ~ time), "random slopes are available")
  expect_error(panglm(y ~ x, data = d, index = idx, model = "random", random = ~ time,
                      weights = rep(2, nrow(d))), "weights")
  d$time_na <- d$time; d$time_na[3] <- NA
  expect_error(panglm(y ~ x, data = d, index = idx, model = "random", random = ~ time_na),
               "missing values")
})
