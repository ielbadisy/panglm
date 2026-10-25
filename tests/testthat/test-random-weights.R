make_re_data <- function(G = 150, Tt = 4, seed = 4) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(G), each = Tt), t = rep(seq_len(Tt), G), x = rnorm(G * Tt))
  u <- rep(rnorm(G, 0, 0.6), each = Tt)
  d$yc <- rpois(G * Tt, exp(0.3 + 0.5 * d$x + u))
  d$yb <- rbinom(G * Tt, 1, plogis(-0.2 + 0.8 * d$x + u))
  d$yg <- rgamma(G * Tt, shape = 2, rate = 2 / exp(1 + 0.4 * d$x + u))
  d$yn <- 1 + 0.5 * d$x + u + rnorm(G * Tt)
  d$w <- rep(sample(1:3, G, TRUE), each = Tt)
  d
}

duplicate_individuals <- function(d) {
  dd <- d[rep(seq_len(nrow(d)), d$w), ]
  dd$id <- paste(dd$id, stats::ave(dd$t, dd$id, dd$t, FUN = seq_along))
  dd
}

# the counts are Poisson given the effect, so the random-effects negative
# binomial is at its boundary (a -> Inf) and warns
fit_case <- function(family, ...) {
  if (family == "negbin") {
    expect_warning(m <- panglm(..., family = family), "diverges")
    m
  } else {
    panglm(..., family = family)
  }
}

cases <- list(c("yc", "poisson"), c("yc", "negbin"), c("yb", "binomial"),
              c("yb", "probit"), c("yg", "gamma"), c("yn", "gaussian"))

test_that("R individual likelihoods reproduce the compiled random-effects fits", {
  d <- make_re_data()
  for (cs in cases) {
    m <- fit_case(cs[2], stats::as.formula(paste(cs[1], "~ x")), d, index = c("id", "t"),
                  model = "random")
    grp <- rep(seq_along(m$group_size), m$group_size)
    contrib <- panglm:::random_contrib(m$family, m$X, m$y, grp, m$quadrature_nodes %||% 21L)
    ll <- sum(contrib(panglm:::random_theta(m, m$family), FALSE)$ll)
    expect_equal(ll, m$loglik, tolerance = 1e-8, info = cs[2])
  }
})

test_that("integer individual weights equal duplicated individuals", {
  d <- make_re_data()
  dd <- duplicate_individuals(d)
  for (cs in cases) {
    f <- stats::as.formula(paste(cs[1], "~ x"))
    mw <- fit_case(cs[2], f, d, index = c("id", "t"), model = "random", weights = d$w)
    md <- fit_case(cs[2], f, dd, index = c("id", "t"), model = "random",
                   weights = rep(1 + 1e-9, nrow(dd)))
    expect_equal(coef(mw)[-1], coef(md)[-1], tolerance = 1e-5, info = cs[2])
    expect_true(all(is.finite(vcov(mw, type = "cluster"))), info = cs[2])
  }
})

test_that("weighted Gaussian random effects is maximum likelihood (lmer)", {
  skip_if_not_installed("lme4")
  d <- make_re_data()
  m <- panglm(yn ~ x, d, index = c("id", "t"), model = "random", family = "gaussian",
              weights = rep(1 + 1e-9, nrow(d)))
  l <- lme4::lmer(yn ~ x + (1 | id), d, REML = FALSE)
  expect_equal(unname(coef(m)), unname(lme4::fixef(l)), tolerance = 1e-5)
  expect_equal(m$loglik, as.numeric(stats::logLik(l)), tolerance = 1e-5)
})

test_that("Swamy-Arora cluster vcov is the quasi-demeaned sandwich", {
  skip_if_not_installed("plm")
  d <- make_re_data()
  m <- panglm(yn ~ x, d, index = c("id", "t"), model = "random", family = "gaussian")
  p <- plm::plm(yn ~ x, d, index = c("id", "t"), model = "random")
  n <- nrow(d); k <- 2
  ref <- plm::vcovHC(p, method = "arellano", type = "sss") / ((n - 1) / (n - k))
  expect_equal(unname(vcov(m, type = "cluster")), matrix(ref, 2, 2), tolerance = 1e-6)
})

test_that("random-effects vcov rules", {
  d <- make_re_data()
  m <- panglm(yc ~ x, d, index = c("id", "t"), model = "random", family = "poisson",
              vcov = "cluster")
  expect_identical(m$vcov_type, "cluster")
  expect_error(vcov(m, type = "HC1"), "HC1")
  expect_error(vcov(m, type = "cluster", cluster = seq_len(nrow(d))), "whole individuals")
  d$w2 <- d$w + (d$t == 2)
  expect_error(panglm(yc ~ x, d, index = c("id", "t"), model = "random",
                      family = "poisson", weights = d$w2), "constant within individuals")
})

test_that("zero-weight individuals are dropped from estimation", {
  d <- make_re_data()
  d$w0 <- ifelse(d$id <= 20, 0, 1.5)
  m0 <- panglm(yg ~ x, d, index = c("id", "t"), model = "random", family = "gamma", weights = d$w0)
  m1 <- panglm(yg ~ x, d[d$id > 20, ], index = c("id", "t"), model = "random", family = "gamma",
               weights = rep(1.5, sum(d$id > 20)))
  expect_equal(coef(m0), coef(m1), tolerance = 1e-6)
})

test_that("negbin random effects falls back to the generic optimizer", {
  set.seed(1)
  d <- make_re_data(G = 200, Tt = 5)
  d$y <- rnbinom(nrow(d), mu = exp(0.3 + 0.5 * d$x + rep(rnorm(200, 0, 0.6), each = 5)), size = 1.5)
  m <- panglm(y ~ x, d, index = c("id", "t"), model = "random", family = "negbin")
  expect_true(is.finite(m$loglik) && m$loglik < 0)
  expect_true(all(is.finite(sqrt(diag(vcov(m))))))
})
