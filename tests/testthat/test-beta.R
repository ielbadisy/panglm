make_beta_data <- function(G = 150, Tt = 4, seed = 7) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(G), each = Tt), t = rep(seq_len(Tt), G),
                  x = rnorm(G * Tt), z = rbinom(G * Tt, 1, 0.4))
  mu <- plogis(0.2 + 0.6 * d$x - 0.4 * d$z + rep(rnorm(G, 0, 0.5), each = Tt))
  d$y <- rbeta(G * Tt, mu * 15, (1 - mu) * 15)
  d$w <- rep(sample(1:3, G, TRUE), each = Tt)
  d
}

test_that("pooled beta regression matches betareg", {
  skip_if_not_installed("betareg")
  d <- make_beta_data()
  m <- panglm(y ~ x + z, d, index = c("id", "t"), model = "pooling", family = "beta")
  b <- betareg::betareg(y ~ x + z, d)
  expect_equal(unname(coef(m)), unname(coef(b)[1:3]), tolerance = 1e-6)
  expect_equal(m$precision, unname(coef(b)[4]), tolerance = 1e-6)
  expect_equal(m$loglik, as.numeric(stats::logLik(b)), tolerance = 1e-8)
  mw <- panglm(y ~ x + z, d, index = c("id", "t"), model = "pooling", family = "beta",
               weights = d$w)
  bw <- betareg::betareg(y ~ x + z, d, weights = w)
  expect_equal(unname(coef(mw)), unname(coef(bw)[1:3]), tolerance = 1e-6)
  pr <- betareg::betareg(y ~ x + z, d, link = "probit")
  mp <- panglm(y ~ x + z, d, index = c("id", "t"), family = beta_family("probit"))
  expect_equal(unname(coef(mp)), unname(coef(pr)[1:3]), tolerance = 1e-6)
})

test_that("random-intercept beta regression is exact quadrature and close to glmmTMB", {
  d <- make_beta_data()
  m <- panglm(y ~ x + z, d, index = c("id", "t"), model = "random", family = "beta")
  th <- m$theta
  eta <- as.numeric(m$X %*% th[1:3]); phi <- exp(th[4]); s <- exp(th[5])
  grp <- rep(seq_along(m$group_size), m$group_size)
  ll <- sum(vapply(seq_along(m$group_size), function(i) {
    e <- eta[grp == i]; yy <- m$y[grp == i]
    f <- function(u) vapply(u, function(v) {
      p <- plogis(e + v)
      exp(sum(stats::dbeta(yy, p * phi, (1 - p) * phi, log = TRUE)))
    }, numeric(1)) * stats::dnorm(u, 0, s)
    log(stats::integrate(f, -Inf, Inf, rel.tol = 1e-10)$value)
  }, numeric(1)))
  expect_equal(m$loglik, ll, tolerance = 1e-6)
  skip_if_not_installed("glmmTMB")
  g <- glmmTMB::glmmTMB(y ~ x + z + (1 | id), d, family = glmmTMB::beta_family())
  expect_equal(unname(coef(m)), unname(glmmTMB::fixef(g)$cond), tolerance = 2e-3)
  expect_equal(m$sigma, sqrt(glmmTMB::VarCorr(g)$cond$id[1]), tolerance = 1e-2)
})

test_that("beta fits support weights, cluster vcov, CRE, predictions and APEs", {
  d <- make_beta_data()
  dd <- d[rep(seq_len(nrow(d)), d$w), ]
  dd$id <- paste(dd$id, stats::ave(dd$t, dd$id, dd$t, FUN = seq_along))
  mw <- panglm(y ~ x + z, d, index = c("id", "t"), model = "random", family = "beta",
               weights = d$w, vcov = "cluster")
  md <- panglm(y ~ x + z, dd, index = c("id", "t"), model = "random", family = "beta")
  expect_equal(coef(mw), coef(md), tolerance = 1e-5)
  expect_true(all(is.finite(vcov(mw))))
  expect_error(vcov(mw, type = "HC1"), "HC1")
  mc <- panglm(y ~ x + z, d, index = c("id", "t"), model = "random", family = "beta", cre = TRUE)
  expect_true(all(c("mean_x", "mean_z") %in% names(coef(mc))))
  p <- predict(md, newdata = d[1:5, ], type = "response")
  expect_true(all(p > 0 & p < 1))
  ape <- partialeffects(md)
  expect_equal(nrow(ape), 2L)
  expect_error(panglm(y ~ x, transform(d, y = round(y)), index = c("id", "t"),
                      family = "beta"), "strictly between 0 and 1")
  expect_error(panglm(y ~ x, d, index = c("id", "t"), model = "within", family = "beta"),
               "no fixed-effects beta")
})
