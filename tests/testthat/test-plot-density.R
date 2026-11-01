with_pdf <- function(code) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  force(code)
}

ix <- c("id", "visit")

total_mass <- function(d, lower, upper) {
  stats::integrate(d$f, lower, upper, rel.tol = 1e-10, subdivisions = 1000L)$value +
    sum(d$masses$predicted)
}

copd_censored <- function() {
  data(copd, envir = environment())
  copd$yc <- pmax(copd$fev1, 1.5)
  copd
}

test_that("the pooled Gaussian implied density is the average normal density", {
  data(copd)
  fit <- panglm(fev1 ~ crp + age, data = copd, index = ix, model = "pooling")
  d <- panglm:::panglm_outcome_density(fit)
  grid <- seq(0.5, 4, length.out = 50)
  sigma <- summary(stats::lm(fev1 ~ crp + age, data = copd))$sigma
  ref <- vapply(grid, function(t) mean(stats::dnorm(t, fit$fitted.values, sigma)), 0)
  expect_equal(d$f(grid), ref, tolerance = 1e-12)
  expect_equal(total_mass(d, -Inf, Inf), 1, tolerance = 1e-3)
})

test_that("implied densities integrate to one with their point masses", {
  data(copd)
  re <- panglm(fev1 ~ crp, data = copd, index = ix, model = "random")
  expect_equal(total_mass(panglm:::panglm_outcome_density(re), -Inf, Inf), 1, tolerance = 1e-3)
  fe <- panglm(fev1 ~ crp, data = copd, index = ix, model = "within")
  expect_equal(total_mass(panglm:::panglm_outcome_density(fe), -Inf, Inf), 1, tolerance = 1e-3)
  gp <- panglm(fev1 ~ crp, data = copd, index = ix, model = "pooling", family = "gamma")
  expect_equal(total_mass(panglm:::panglm_outcome_density(gp), 0, Inf), 1, tolerance = 1e-3)
  gr <- panglm(fev1 ~ crp, data = copd, index = ix, model = "random", family = "gamma")
  expect_equal(total_mass(panglm:::panglm_outcome_density(gr), 0, Inf), 1, tolerance = 1e-3)
  copd$y01 <- (copd$fev1 - min(copd$fev1) + 0.05) / (diff(range(copd$fev1)) + 0.1)
  br <- panglm(y01 ~ crp, data = copd, index = ix, model = "random", family = "beta")
  expect_equal(total_mass(panglm:::panglm_outcome_density(br), 0, 1), 1, tolerance = 1e-3)

  cc <- copd_censored()
  tp <- panglm(yc ~ crp, data = cc, index = ix, model = "pooling",
               family = tobitfamily(left = 1.5))
  dt <- panglm:::panglm_outcome_density(tp)
  expect_equal(total_mass(dt, 1.5, Inf), 1, tolerance = 1e-3)
  tr <- panglm(yc ~ crp, data = cc, index = ix, model = "random",
               family = tobitfamily(left = 1.5))
  expect_equal(total_mass(panglm:::panglm_outcome_density(tr), 1.5, Inf), 1, tolerance = 1e-3)
  cc$yc2 <- pmin(cc$yc, 3)
  t2 <- panglm(yc2 ~ crp, data = cc, index = ix, model = "pooling",
               family = tobitfamily(left = 1.5, right = 3))
  expect_equal(total_mass(panglm:::panglm_outcome_density(t2), 1.5, 3), 1, tolerance = 1e-3)
})

test_that("the Tobit point mass is the average censoring probability", {
  cc <- copd_censored()
  fit <- panglm(yc ~ crp + age, data = cc, index = ix, model = "pooling",
                family = tobitfamily(left = 1.5))
  d <- panglm:::panglm_outcome_density(fit)
  expect_equal(d$masses$predicted,
               mean(stats::pnorm((1.5 - fit$linear.predictors) / fit$sigma_e)), tolerance = 1e-12)
  expect_equal(d$masses$observed, mean(cc$yc <= 1.5))
  skip_if_missing("AER")
  ref <- AER::tobit(yc ~ crp + age, left = 1.5, data = cc)
  lp <- stats::predict(ref, type = "lp")
  expect_equal(d$masses$predicted, mean(stats::pnorm((1.5 - lp) / ref$scale)), tolerance = 1e-6)
})

test_that("the pooled gamma density uses the dispersion of stats::glm()", {
  data(copd)
  fit <- panglm(fev1 ~ crp + age, data = copd, index = ix, model = "pooling", family = "gamma")
  ref <- stats::glm(fev1 ~ crp + age, data = copd, family = stats::Gamma("log"))
  phi <- summary(ref)$dispersion
  grid <- seq(0.5, 4, length.out = 40)
  expect_equal(panglm:::panglm_outcome_density(fit)$f(grid),
               vapply(grid, function(t) mean(stats::dgamma(t, shape = 1 / phi,
                                                             rate = 1 / (phi * stats::fitted(ref)))), 0),
               tolerance = 1e-6)
})

test_that("two-part implied distributions integrate to one", {
  data(copd)
  set.seed(1)
  copd$cost <- ifelse(copd$exacerbations > 0,
                      exp(6 + 0.1 * copd$crp + stats::rnorm(nrow(copd), 0, 0.5)), 0)
  for (pt in c("gamma", "lognormal")) {
    fit <- twopart(cost ~ treatment + crp, data = copd, index = ix, positive = pt)
    d <- panglm:::twopart_outcome_density(fit)
    expect_equal(d$masses$predicted, mean(1 - predict(fit, type = "prob")), tolerance = 1e-12)
    expect_equal(d$masses$observed, mean(copd$cost == 0))
    expect_equal(total_mass(d, 0, Inf), 1, tolerance = 1e-3)
  }
  re <- twopart(cost ~ treatment + crp, data = copd, index = ix, model = "random",
                positive = "lognormal", R = 7)
  expect_equal(total_mass(panglm:::twopart_outcome_density(re), 0, Inf), 1, tolerance = 1e-3)
  # correlated intercepts: with rho = 0 the joint density is the independent one
  cr <- re
  cr$correlated <- TRUE
  cr$joint <- list(sigma_binary = re$binary$sigma, sigma_positive = sqrt(re$positive$sigma_mu2),
                   sigma_e = sqrt(re$positive$sigma_v2), rho = 0)
  grid <- c(50, 200, 800, 2000)
  expect_equal(panglm:::twopart_outcome_density(cr)$f(grid),
               panglm:::twopart_outcome_density(re)$f(grid), tolerance = 1e-10)
  cr$joint$rho <- 0.6
  expect_equal(total_mass(panglm:::twopart_outcome_density(cr), 0, Inf), 1, tolerance = 1e-3)
})

test_that("the density view draws and refuses discrete outcomes", {
  data(copd)
  g <- panglm(fev1 ~ crp, data = copd, index = ix, model = "random")
  cc <- copd_censored()
  t1 <- panglm(yc ~ crp, data = cc, index = ix, model = "pooling", family = tobitfamily(left = 1.5))
  pois <- panglm(exacerbations ~ crp, data = copd, index = ix, model = "pooling",
                 family = "poisson")
  bin <- panglm(hospitalized ~ crp, data = copd, index = ix, model = "pooling",
                family = "binomial")
  bt <- panglm(fev1 ~ crp, data = copd, index = ix, model = "between")
  with_pdf({
    old <- graphics::par("mar")
    expect_invisible(plot(g, which = "density"))
    expect_identical(plot(t1, which = "density", bw = 0.1), t1)
    expect_equal(graphics::par("mar"), old)
    expect_error(plot(pois, which = "density"), "rootogram")
    expect_error(plot(bin, which = "density"), "calibration")
    expect_error(plot(bt, which = "density"), "between")
  })
  set.seed(1)
  copd$cost <- ifelse(copd$exacerbations > 0, stats::rgamma(nrow(copd), 2, 0.01), 0)
  tp <- twopart(cost ~ crp, data = copd, index = ix)
  with_pdf({
    expect_identical(plot(tp, which = "density"), tp)
    d <- panglm:::twopart_outcome_density(tp)
    expect_true(panglm:::.panglm_plot_density(d, log = "x")$log)
    out <- panglm:::.panglm_plot_density(d, log = "")
    expect_false(out$log)
    expect_equal(out$implied, d$f(out$x), tolerance = 1e-4)
    expect_error(panglm:::.panglm_plot_density(panglm:::panglm_outcome_density(g), log = "x"),
                 "positive")
  })
})
