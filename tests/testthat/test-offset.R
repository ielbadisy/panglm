rate_data <- function() {
  set.seed(21)
  d <- data.frame(id = rep(1:80, each = 3), t = rep(1:3, 80))
  d$x <- rnorm(240); d$g <- factor(sample(c("a", "b", "c"), 240, TRUE))
  d$py <- runif(240, 0.2, 3)
  d$n <- rpois(240, d$py * exp(-0.5 + 0.4 * d$x + 0.3 * (d$g == "b")))
  d$lpy <- log(d$py)
  d
}
ctl <- stats::glm.control(epsilon = 1e-14, maxit = 100)

test_that("pooled Poisson rate model with an offset matches glm()", {
  d <- rate_data()
  ref <- glm(n ~ x + g, poisson, d, offset = log(py), control = ctl)
  fits <- list(panglm(n ~ x + g, d, index = c("id", "t"), family = "poisson", offset = "lpy"),
               panglm(n ~ x + g, d, index = c("id", "t"), family = "poisson", offset = d$lpy),
               panglm(n ~ x + g + offset(log(py)), d, index = c("id", "t"), family = "poisson"))
  for (f in fits) {
    expect_equal(coef(f), coef(ref), tolerance = 1e-8)
    expect_equal(vcov(f), vcov(ref), tolerance = 1e-8)
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(ref)), tolerance = 1e-10)
    expect_equal(unname(fitted(f)), unname(fitted(ref)), tolerance = 1e-8)
    expect_equal(unname(predict(f)), unname(predict(ref)), tolerance = 1e-8)
  }
  nd <- d[c(5, 50, 100), ]
  for (f in fits[-2]) {
    expect_equal(unname(predict(f, newdata = nd, type = "response")),
                 unname(predict(ref, nd, type = "response")), tolerance = 1e-8)
  }
  expect_error(predict(fits[[2]], newdata = nd), "numeric vector")
  skip_if_missing("sandwich")
  h <- panglm(n ~ x + g, d, index = c("id", "t"), family = "poisson", offset = "lpy",
              vcov = "HC1")
  expect_equal(vcov(h), sandwich::vcovHC(ref, type = "HC1"), tolerance = 1e-8)
  cl <- panglm(n ~ x + g, d, index = c("id", "t"), family = "poisson", offset = "lpy",
               vcov = "cluster")
  expect_equal(unname(vcov(cl)), unname(sandwich::vcovCL(ref, cluster = d$id, type = "HC1")),
               tolerance = 1e-8)
})

test_that("offsets match glm() for the gaussian, binomial and gamma families", {
  d <- rate_data(); d$o <- d$x / 3
  d$yb <- rbinom(240, 1, plogis(0.2 + d$x)); d$yg <- rgamma(240, 2, 2 / exp(1 + 0.2 * d$x))
  d$yn <- 2 + d$x + rnorm(240)
  cases <- list(list("gaussian", gaussian(), "yn"), list("binomial", binomial(), "yb"),
                list("probit", binomial("probit"), "yb"), list("gamma", Gamma("log"), "yg"))
  for (cs in cases) {
    f <- reformulate(c("x", "g"), cs[[3]])
    a <- panglm(f, d, index = c("id", "t"), family = cs[[1]], offset = "o")
    b <- glm(f, cs[[2]], d, offset = o, control = ctl)
    expect_equal(coef(a), coef(b), tolerance = 1e-7, info = cs[[1]])
    expect_equal(vcov(a), vcov(b), tolerance = 1e-6, info = cs[[1]])
    expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-8, info = cs[[1]])
  }
})

test_that("negative binomial with an offset matches MASS::glm.nb()", {
  skip_if_missing("MASS")
  d <- rate_data()
  d$n <- MASS::rnegbin(240, d$py * exp(-0.5 + 0.4 * d$x), theta = 1.5)
  ref <- MASS::glm.nb(n ~ x + g + offset(lpy), d, control = glm.control(1e-14, 100))
  a <- panglm(n ~ x + g, d, index = c("id", "t"), family = "negbin", offset = "lpy")
  # glm.nb estimates theta to about 1e-4 (theta.ml's default tolerance)
  expect_equal(coef(a), coef(ref), tolerance = 1e-5)
  expect_equal(a$theta, ref$theta, tolerance = 1e-4)
  expect_gte(as.numeric(logLik(a)), as.numeric(logLik(ref)) - 1e-8)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(ref)), tolerance = 1e-8)
})

test_that("partial effects of an offset model hold the offset", {
  d <- rate_data()
  f <- panglm(n ~ x + g, d, index = c("id", "t"), family = "poisson", offset = "lpy")
  ref <- glm(n ~ x + g, poisson, d, offset = lpy, control = ctl)
  ape <- partialeffects(f)
  expect_equal(ape$estimate[ape$term == "x"], coef(ref)[["x"]] * mean(fitted(ref)),
               tolerance = 1e-6)
})

test_that("a zero offset leaves the compiled fit unchanged", {
  d <- rate_data(); d$z <- 0
  a <- panglm(n ~ x, d, index = c("id", "t"), family = "poisson", offset = "z")
  b <- panglm(n ~ x, d, index = c("id", "t"), family = "poisson")
  expect_identical(coef(a), coef(b)); expect_null(a$offset)
})
