with_pdf <- function(code) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  force(code)
}

ix <- c("id", "visit")

test_that("the effects view draws fixed and random effects", {
  data(copd)
  fe <- panglm(fev1 ~ crp, data = copd, index = ix, model = "within")
  re <- panglm(fev1 ~ crp, data = copd, index = ix, model = "random")
  hb <- panglm(hospitalized ~ crp, data = copd, index = ix, model = "random",
               family = "binomial", R = 7)
  pr <- panglm(exacerbations ~ crp, data = copd, index = ix, model = "random",
               family = "poisson")
  po <- panglm(fev1 ~ crp, data = copd, index = ix, model = "pooling")
  cl <- panglm(hospitalized ~ crp, data = copd, index = ix, model = "within",
               family = "binomial")
  with_pdf({
    old <- graphics::par(c("mfrow", "mar"))
    expect_invisible(plot(fe, which = "effects"))
    expect_identical(plot(re, which = "effects", main = "RE"), re)
    expect_identical(plot(hb, which = "effects"), hb)
    expect_identical(plot(pr, which = "effects"), pr)
    expect_equal(graphics::par(c("mfrow", "mar")), old)
    expect_error(plot(po, which = "effects"), "within")
    expect_error(plot(cl, which = "effects"), "conditional logit")
  })
})

test_that("random-slope effects get one panel per effect plus scatters", {
  set.seed(3)
  G <- 60; Tn <- 6
  d <- data.frame(id = rep(seq_len(G), each = Tn), t = rep(seq_len(Tn), G))
  b0 <- rnorm(G, 0, 0.8); b1 <- 0.5 * b0 + rnorm(G, 0, 0.2)
  d$x <- rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + b0[d$id] + (0.3 + b1[d$id]) * d$t / Tn + rnorm(nrow(d), 0, 0.3)
  d$tt <- d$t / Tn
  fit <- panglm(y ~ x + tt, data = d, index = c("id", "t"), model = "random", random = ~ tt)
  with_pdf({
    old <- graphics::par("mfrow")
    expect_identical(plot(fit, which = "effects"), fit)
    expect_equal(graphics::par("mfrow"), old)
  })
})

test_that("assumed mixing densities integrate to one", {
  data(copd)
  pr <- panglm(exacerbations ~ crp, data = copd, index = ix, model = "random",
               family = "poisson")
  p <- panglm:::panglm_effect_prior(pr)
  expect_equal(stats::integrate(p$density, -Inf, Inf)$value, 1, tolerance = 1e-6)
  fake <- list(family = list(family = "negbin"), negbin_a = 3, negbin_b = 4)
  p <- panglm:::panglm_effect_prior(fake)
  expect_equal(stats::integrate(p$density, -Inf, Inf)$value, 1, tolerance = 1e-6)
  m <- stats::integrate(function(v) v * p$density(v), -Inf, Inf)$value
  expect_equal(m, digamma(4) - digamma(3), tolerance = 1e-6)
})
