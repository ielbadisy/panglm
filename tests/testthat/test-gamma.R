test_that("pooled gamma (weighted and robust) matches glm()", {
  d <- make_cost_data()
  set.seed(2); w <- stats::runif(nrow(d), 0.5, 2)
  f <- panglm(cost ~ treatment + age + crp, d, c("id", "visit"), family = "gamma", weights = w,
              tol = 1e-14)
  g <- glm(cost ~ treatment + age + crp, Gamma("log"), d, weights = w,
           control = glm.control(epsilon = 1e-14))
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-7)
  expect_equal(unname(vcov(f)), unname(vcov(g)), tolerance = 1e-7)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-7)
  expect_equal(AIC(f), AIC(g), tolerance = 1e-7)

  skip_if_missing("sandwich")
  f <- panglm(cost ~ treatment + age + crp, d, c("id", "visit"), family = "gamma", vcov = "HC1")
  g <- glm(cost ~ treatment + age + crp, Gamma("log"), d)
  expect_equal(unname(vcov(f)), unname(sandwich::vcovHC(g, "HC1")), tolerance = 1e-6)
})

test_that("gamma requires a positive outcome", {
  d <- make_cost_data()
  expect_error(panglm(cost0 ~ crp, d, c("id", "visit"), family = "gamma"), "positive")
})

test_that("fixed-effects gamma matches glm() with individual dummies", {
  d <- make_cost_data()
  f <- panglm(cost ~ crp, d, c("id", "visit"), model = "within", family = "gamma",
              tol = 1e-14)
  g <- glm(cost ~ crp + factor(id), Gamma("log"), d, control = glm.control(epsilon = 1e-14))
  expect_equal(unname(coef(f)), unname(coef(g)["crp"]), tolerance = 1e-7)
  expect_equal(as.numeric(sqrt(vcov(f, type = "classical"))),
               unname(sqrt(vcov(g)["crp", "crp"])), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-6)
  expect_equal(unname(fitted(f)[order(order(d$id, d$visit))]), unname(fitted(g)),
               tolerance = 1e-6)
})

test_that("random-intercept gamma with R = 1 matches glmmTMB (Laplace)", {
  skip_if_missing("glmmTMB")
  skip_on_cran()
  d <- make_cost_data()
  f <- panglm(cost ~ treatment + age + crp, d, c("id", "visit"), model = "random",
              family = "gamma", R = 1)
  g <- glmmTMB::glmmTMB(cost ~ treatment + age + crp + (1 | id), d, family = Gamma("log"))
  expect_equal(unname(coef(f)), unname(glmmTMB::fixef(g)$cond), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-5)
  expect_equal(f$sigma, sqrt(glmmTMB::VarCorr(g)$cond$id[1]), tolerance = 1e-3)
  expect_equal(f$shape, 1 / sigma(g)^2, tolerance = 1e-3)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(vcov(g)$cond))), tolerance = 1e-3)
})

test_that("random-intercept gamma predictions include the lognormal factor", {
  d <- make_cost_data()
  f <- panglm(cost ~ treatment + crp, d, c("id", "visit"), model = "random", family = "gamma")
  eta <- as.numeric(cbind(1, d$treatment, d$crp) %*% coef(f))
  expect_equal(unname(predict(f, d, type = "response")), exp(eta + f$sigma^2 / 2),
               tolerance = 1e-8)
  expect_output(print(summary(f)), "shape")
})
