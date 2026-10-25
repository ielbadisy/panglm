make_tobit_data <- function(G = 120, Tt = 4, seed = 11) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(G), each = Tt), t = rep(seq_len(Tt), G),
                  x = rnorm(G * Tt), z = rbinom(G * Tt, 1, 0.4))
  ystar <- 0.3 + 0.8 * d$x - 0.5 * d$z + rep(rnorm(G, 0, 0.8), each = Tt) + rnorm(G * Tt)
  d$y <- pmax(0, ystar)
  d$y2 <- pmin(2, pmax(0, ystar))
  d$w <- rep(sample(1:3, G, TRUE), each = Tt)
  d
}

test_that("pooled Tobit matches AER::tobit with one and two limits", {
  skip_if_not_installed("AER")
  d <- make_tobit_data()
  m <- panglm(y ~ x + z, d, index = c("id", "t"), model = "pooling", family = tobit_family(0))
  a <- AER::tobit(y ~ x + z, left = 0, data = d)
  expect_equal(unname(coef(m)), unname(coef(a)), tolerance = 1e-6)
  expect_equal(m$sigma_e, a$scale, tolerance = 1e-6)
  expect_equal(m$loglik, as.numeric(stats::logLik(a)), tolerance = 1e-8)
  expect_equal(unname(sqrt(diag(vcov(m)))), unname(sqrt(diag(vcov(a)))[1:3]), tolerance = 1e-5)
  m2 <- panglm(y2 ~ x + z, d, index = c("id", "t"), family = tobit_family(0, 2))
  a2 <- AER::tobit(y2 ~ x + z, left = 0, right = 2, data = d)
  expect_equal(unname(coef(m2)), unname(coef(a2)), tolerance = 1e-6)
  expect_equal(m2$loglik, as.numeric(stats::logLik(a2)), tolerance = 1e-8)
  expect_equal(unname(m2$censoring), c(sum(d$y2 == 0), sum(d$y2 > 0 & d$y2 < 2), sum(d$y2 == 2)))
})

test_that("weighted pooled Tobit matches survreg with case weights", {
  skip_if_not_installed("AER")
  d <- make_tobit_data()
  m <- panglm(y ~ x + z, d, index = c("id", "t"), family = tobit_family(0), weights = d$w)
  a <- AER::tobit(y ~ x + z, left = 0, data = d, weights = w)
  expect_equal(unname(coef(m)), unname(coef(a)), tolerance = 1e-6)
})

test_that("random-intercept Tobit matches censReg", {
  skip_if_not_installed("censReg")
  skip_if_not_installed("plm")
  d <- make_tobit_data()
  m <- panglm(y ~ x + z, d, index = c("id", "t"), model = "random", family = tobit_family(0))
  pd <- plm::pdata.frame(d, index = c("id", "t"))
  r <- censReg::censReg(y ~ x + z, left = 0, data = pd, nGHQ = 30)
  expect_equal(unname(coef(m)), unname(coef(r)[1:3]), tolerance = 1e-4)
  expect_equal(c(log(m$sigma), log(m$sigma_e)), unname(coef(r)[4:5]), tolerance = 1e-4)
  expect_equal(m$loglik, as.numeric(stats::logLik(r)), tolerance = 1e-6)
  expect_true(all(is.finite(sqrt(diag(vcov(m, type = "cluster"))))))
  expect_error(vcov(m, type = "HC1"), "cluster")
})

test_that("Tobit predictions and partial effects are on the censored-mean scale", {
  d <- make_tobit_data()
  m <- panglm(y ~ x + z, d, index = c("id", "t"), family = tobit_family(0))
  eta <- predict(m, type = "link")
  s <- m$sigma_e
  expect_equal(predict(m, type = "response"),
               pnorm(eta / s) * eta + s * dnorm(eta / s), tolerance = 1e-10)
  ape <- partialeffects(m)
  expect_equal(ape$estimate[ape$term == "x"], mean(pnorm(eta / s)) * coef(m)[["x"]],
               tolerance = 1e-6)
})

test_that("Tobit input checks", {
  d <- make_tobit_data()
  expect_error(tobit_family(2, 1), "left < right")
  expect_error(tobit_family(-Inf, Inf), "finite")
  expect_error(panglm(y ~ x, d, index = c("id", "t"), model = "within", family = tobit_family(0)),
               "fixed-effects Tobit")
  expect_error(panglm(y ~ x, d, index = c("id", "t"), family = tobit_family(1)), "within")
  expect_s3_class(panglm(y ~ x, d, index = c("id", "t"), family = "tobit"), "panglm")
})
