weights_data <- function() {
  data(copd, package = "panglm", envir = environment())
  set.seed(1)
  w <- stats::runif(nrow(copd), 0.2, 3)
  list(d = copd, w = w)
}

test_that("weighted pooled probit and gaussian match glm()", {
  x <- weights_data(); d <- x$d; w <- x$w
  f <- panglm(hospitalized ~ treatment + age + crp, d, c("id", "visit"),
              family = "probit", weights = w)
  g <- suppressWarnings(glm(hospitalized ~ treatment + age + crp,
                            binomial("probit"), d, weights = w, control = glm.control(epsilon = 1e-14)))
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-6)
  expect_equal(unname(vcov(f)), unname(vcov(g)), tolerance = 1e-5)

  f <- panglm(fev1 ~ treatment + age + crp, d, c("id", "visit"), weights = w)
  g <- glm(fev1 ~ treatment + age + crp, data = d, weights = w)
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-8)
  expect_equal(unname(vcov(f)), unname(vcov(g)), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-8)
})

test_that("zero weights drop observations", {
  x <- weights_data(); d <- x$d; w <- x$w
  w[1:5] <- 0
  f <- panglm(exacerbations ~ treatment + crp, d, c("id", "visit"),
              family = "poisson", weights = w)
  g <- glm(exacerbations ~ treatment + crp, poisson, d[-(1:5), ], weights = w[-(1:5)])
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-7)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-7)
})

test_that("weighted robust and cluster covariances match sandwich", {
  skip_if_missing("sandwich")
  x <- weights_data(); d <- x$d; w <- x$w
  g <- glm(exacerbations ~ treatment + age + crp, poisson, d, weights = w)
  f <- panglm(exacerbations ~ treatment + age + crp, d, c("id", "visit"),
              family = "poisson", weights = w, vcov = "HC1")
  expect_equal(unname(vcov(f)), unname(sandwich::vcovHC(g, "HC1")), tolerance = 1e-6)
  f <- panglm(exacerbations ~ treatment + age + crp, d, c("id", "visit"),
              family = "poisson", weights = w, vcov = "cluster")
  expect_equal(unname(vcov(f)),
               unname(sandwich::vcovCL(g, cluster = ~id, type = "HC1")), tolerance = 1e-6)
})

test_that("weighted pooled negative binomial matches MASS::glm.nb()", {
  skip_if_missing("MASS")
  skip_if_missing("sandwich")
  x <- weights_data(); d <- x$d
  w <- round(x$w * 2) + 1
  f <- panglm(exacerbations ~ treatment + crp, d, c("id", "visit"),
              family = "negbin", weights = w, vcov = "HC1")
  g <- MASS::glm.nb(exacerbations ~ treatment + crp, d, weights = w)
  expect_equal(unname(coef(f)), unname(coef(g)), tolerance = 1e-5)
  expect_equal(f$theta, g$theta, tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(g)), tolerance = 1e-6)
  expect_equal(unname(vcov(f)), unname(sandwich::vcovHC(g, "HC1")), tolerance = 1e-4)
})

test_that("unweighted pooled negative binomial robust vcov matches sandwich", {
  skip_if_missing("MASS")
  skip_if_missing("sandwich")
  x <- weights_data(); d <- x$d
  f <- panglm(exacerbations ~ treatment + crp, d, c("id", "visit"),
              family = "negbin", vcov = "HC1")
  g <- MASS::glm.nb(exacerbations ~ treatment + crp, d)
  expect_equal(unname(vcov(f)), unname(sandwich::vcovHC(g, "HC1")), tolerance = 1e-4)
})
