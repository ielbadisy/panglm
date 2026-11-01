with_pdf <- function(code) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  force(code)
}

count_formula <- exacerbations ~ crp + treatment + age

test_that("the coefficient plot selects terms and keeps the margin", {
  data(copd)
  fit <- panglm(fev1 ~ crp + age + treatment, data = copd, index = c("id", "visit"),
                model = "random")
  with_pdf({
    old <- graphics::par("mai")
    expect_invisible(plot(fit, parm = c("crp", "treatment")))
    expect_identical(plot(fit, parm = 2:3), fit)
    expect_equal(graphics::par("mai"), old)
    expect_error(plot(fit, parm = "zz"), "parm")
  })
})

test_that("Pearson residuals match stats::glm on pooled fits", {
  data(copd)
  p <- panglm(count_formula, data = copd, index = c("id", "visit"), model = "pooling",
              family = "poisson")
  g <- glm(count_formula, data = copd, family = poisson)
  expect_equal(unname(residuals(p, type = "pearson")), unname(residuals(g, type = "pearson")),
               tolerance = 1e-6)
  expect_equal(residuals(p), residuals(p, type = "response"))

  b <- panglm(I(exacerbations > 0) ~ crp + treatment, data = copd,
              index = c("id", "visit"), model = "pooling", family = "binomial")
  gb <- glm(I(exacerbations > 0) ~ crp + treatment, data = copd, family = binomial)
  expect_equal(unname(residuals(b, type = "pearson")), unname(residuals(gb, type = "pearson")),
               tolerance = 1e-6)

  gam <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"), model = "pooling",
                family = "gamma")
  gg <- glm(fev1 ~ crp, data = copd, family = Gamma(link = "log"))
  expect_equal(unname(residuals(gam, type = "pearson")), unname(residuals(gg, type = "pearson")),
               tolerance = 1e-5)

  skip_if_missing("MASS")
  nb <- panglm(count_formula, data = copd, index = c("id", "visit"), model = "pooling",
               family = "negbin")
  gn <- MASS::glm.nb(count_formula, data = copd)
  expect_equal(unname(residuals(nb, type = "pearson")), unname(residuals(gn, type = "pearson")),
               tolerance = 1e-4)
})

test_that("rootogram expected frequencies match an independent computation", {
  data(copd)
  p <- panglm(count_formula, data = copd, index = c("id", "visit"), model = "pooling",
              family = "poisson")
  g <- glm(count_formula, data = copd, family = poisson)
  rg <- panglm:::panglm_rootogram(p, max = 8)
  ref <- vapply(0:8, function(k) sum(dpois(k, fitted(g))), numeric(1))
  expect_equal(rg$expected, ref, tolerance = 1e-6)
  expect_equal(rg$observed, tabulate(copd$exacerbations + 1, 9))

  skip_if_missing("MASS")
  nb <- panglm(count_formula, data = copd, index = c("id", "visit"), model = "pooling",
               family = "negbin")
  gn <- MASS::glm.nb(count_formula, data = copd)
  rg <- panglm:::panglm_rootogram(nb, max = 8)
  ref <- vapply(0:8, function(k) sum(dnbinom(k, size = gn$theta, mu = fitted(gn))), numeric(1))
  expect_equal(rg$expected, ref, tolerance = 1e-4)

  if (requireNamespace("topmodels", quietly = TRUE)) {
    tm <- topmodels::rootogram(gn, breaks = -0.5:8.5, plot = FALSE)
    expect_equal(rg$expected, as.numeric(tm$expected), tolerance = 1e-4)
  }
})

test_that("random-intercept rootograms use proper marginal distributions", {
  data(copd)
  re <- suppressWarnings(panglm(count_formula, data = copd, index = c("id", "visit"),
                                model = "random", family = "poisson"))
  rg <- panglm:::panglm_rootogram(re, max = 200)
  # the marginal probabilities of every observation sum to one, mean = fitted
  expect_equal(sum(rg$expected), nobs(re), tolerance = 1e-6)
  expect_equal(sum(rg$count * rg$expected), sum(re$fitted.values), tolerance = 1e-4)

  renb <- suppressWarnings(panglm(count_formula, data = copd, index = c("id", "visit"),
                                  model = "random", family = "negbin"))
  rg <- panglm:::panglm_rootogram(renb, max = 400)
  expect_equal(sum(rg$expected), nobs(renb), tolerance = 1e-4)
  if (renb$negbin_a > 1) {
    expect_equal(sum(rg$count * rg$expected), sum(renb$fitted.values), tolerance = 1e-3)
  }
  if (renb$negbin_a > 2) {
    # the beta-negative-binomial variance used by the Pearson residuals
    lambda <- exp(renb$linear.predictors[1])
    k <- 0:2000
    pk <- exp(lgamma(lambda + k) - lgamma(lambda) - lgamma(k + 1) +
                lbeta(renb$negbin_a + lambda, renb$negbin_b + k) -
                lbeta(renb$negbin_a, renb$negbin_b))
    m <- sum(k * pk)
    v <- panglm:::panglm_variance(renb, m)
    expect_equal(v, sum((k - m)^2 * pk), tolerance = 1e-3)
  }
})

test_that("every diagnostic view draws and returns the fit", {
  data(copd)
  gau <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"), model = "within")
  pois <- panglm(count_formula, data = copd, index = c("id", "visit"), model = "pooling",
                 family = "poisson")
  fe <- suppressMessages(panglm(exacerbations ~ crp, data = copd, index = c("id", "visit"),
                                model = "within", family = "negbin"))
  bin <- panglm(I(exacerbations > 0) ~ crp + treatment, data = copd,
                index = c("id", "visit"), model = "random", family = "probit")
  with_pdf({
    for (fit in list(gau, pois, fe, bin)) {
      expect_identical(plot(fit, which = "residuals"), fit)
      expect_invisible(plot(fit, which = "residuals", type = "response"))
      expect_invisible(plot(fit, which = "qq"))
      expect_invisible(plot(fit, which = "fitted"))
    }
    expect_invisible(plot(pois, which = "rootogram"))
    expect_invisible(plot(fe, which = "rootogram", max = 5))
    expect_identical(plot(bin, which = "calibration"), bin)
    expect_invisible(plot(bin, which = "calibration", bins = 3, level = 0.9))
    expect_error(plot(gau, which = "calibration"), "binomial")
    expect_error(plot(gau, which = "rootogram"), "Poisson")
    expect_error(plot(bin, which = "calibration", bins = 0), "bins")
  })
})

test_that("calibration bins agree with a direct computation", {
  set.seed(3)
  p <- runif(500)
  y <- rbinom(500, 1, p)
  cal <- panglm:::panglm_calibration(y, p, bins = 5)
  bin <- cut(p, quantile(p, seq(0, 1, 0.2)), include.lowest = TRUE, labels = FALSE)
  expect_equal(cal$observed, as.numeric(tapply(y, bin, mean)))
  expect_equal(cal$predicted, as.numeric(tapply(p, bin, mean)))
  ci <- binom.test(sum(y[bin == 1]), sum(bin == 1))$conf.int
  expect_equal(c(cal$lower[1], cal$upper[1]), as.numeric(ci))
})

test_that("families without Pearson residuals fall back with a message", {
  data(copd)
  tob <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"), model = "pooling",
                family = tobitfamily(left = min(copd$fev1)))
  expect_error(residuals(tob, type = "pearson"), "not defined")
  with_pdf({
    expect_message(plot(tob, which = "residuals"), "response residuals")
    expect_message(plot(tob, which = "qq"), "response residuals")
  })
})

test_that("random-slope fits have no rootogram", {
  data(copd)
  fit <- suppressWarnings(panglm(exacerbations ~ crp + visit, data = copd,
                                 index = c("id", "visit"), model = "random",
                                 family = "poisson", random = ~ visit))
  expect_error(panglm:::panglm_rootogram(fit), "random-slope")
  expect_true(all(is.finite(residuals(fit, type = "pearson"))))
})
