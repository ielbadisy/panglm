test_that("pooled two-part gamma prediction equals the product of glm parts", {
  d <- make_cost_data()
  f <- twopart(cost0 ~ treatment + crp, d, c("id", "visit"), positive = "gamma",
               tol = 1e-14)
  pos <- d$cost0 > 0
  g1 <- glm(pos ~ treatment + crp, binomial("probit"), d, control = glm.control(epsilon = 1e-14))
  g2 <- glm(cost0 ~ treatment + crp, Gamma("log"), d[pos, ], control = glm.control(epsilon = 1e-14))
  expect_equal(unname(coef(f)[grepl("^binary:", names(coef(f)))]), unname(coef(g1)),
               tolerance = 1e-6)
  expect_equal(unname(coef(f)[grepl("^positive:", names(coef(f)))]), unname(coef(g2)),
               tolerance = 1e-6)
  manual <- predict(g1, d, type = "response") * predict(g2, d, type = "response")
  expect_equal(unname(predict(f, d, type = "response")), unname(manual), tolerance = 1e-3)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g1) + logLik(g2)), tolerance = 1e-6)
})

test_that("lognormal two-part uses smearing retransformation", {
  d <- make_cost_data()
  f <- twopart(cost0 ~ treatment + crp, d, c("id", "visit"), positive = "lognormal")
  pos <- d$cost0 > 0
  g2 <- lm(log(cost0) ~ treatment + crp, d[pos, ])
  g1 <- glm(pos ~ treatment + crp, binomial("probit"), d)
  smear <- mean(exp(residuals(g2)))
  manual <- predict(g1, d, type = "response") * exp(predict(g2, d)) * smear
  expect_equal(unname(predict(f, d, type = "response")), unname(manual), tolerance = 1e-6)
})

test_that("random, CRE and weighted two-part variants run", {
  d <- make_cost_data()
  f <- twopart(cost0 ~ treatment + crp, d, c("id", "visit"), model = "random",
               positive = "gamma", cre = TRUE)
  expect_true(all(is.finite(predict(f, d))))
  expect_output(print(summary(f)), "Mundlak")
  a <- partialeffects(f)
  expect_true(all(is.finite(a$estimate)))
  set.seed(3)
  fw <- twopart(cost0 ~ treatment + crp, d, c("id", "visit"),
                weights = stats::runif(nrow(d), 0.5, 2), vcov = "cluster")
  expect_equal(dim(vcov(fw)), rep(length(coef(fw)), 2))
})
