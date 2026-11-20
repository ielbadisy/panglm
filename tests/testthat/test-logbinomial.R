test_that("log-binomial model matches glm(binomial('log'))", {
  skip_if_missing("MASS")
  bw <- MASS::birthwt; bw$id <- seq_len(nrow(bw))
  ref <- glm(low ~ smoke + ht + ui, binomial("log"), bw, start = c(log(mean(bw$low)), 0, 0, 0),
             control = glm.control(1e-14, 100))
  a <- panglm(low ~ smoke + ht + ui, bw, index = "id", family = binomialfamily("log"))
  b <- panglm(low ~ smoke + ht + ui, bw, index = "id", family = binomial("log"), vcov = "HC1")
  expect_true(ref$converged)
  expect_equal(coef(a), coef(ref), tolerance = 1e-7)
  expect_equal(vcov(a), vcov(ref), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(ref)), tolerance = 1e-10)
  expect_equal(coef(b), coef(ref), tolerance = 1e-7)
  skip_if_missing("sandwich")
  expect_equal(vcov(b), sandwich::vcovHC(ref, type = "HC1"), tolerance = 1e-6)
})

test_that("log-binomial boundary maximum matches glm() and logbin", {
  skip_if_missing("MASS")
  bw <- MASS::birthwt; bw$race <- factor(bw$race); bw$id <- seq_len(nrow(bw))
  f <- low ~ smoke + lwt + age + race + ht + ui
  expect_warning(a <- panglm(f, bw, index = "id", family = binomialfamily("log"), maxit = 5000),
                 "maximum on the parameter boundary (fitted risk at 1)", fixed = TRUE)
  expect_true(a$converged); expect_lt(a$iterations, 50); expect_equal(max(fitted(a)), 1)
  st <- c(log(mean(bw$low)), rep(0, 7))
  ref <- suppressWarnings(glm(f, binomial("log"), bw, start = st, control = glm.control(1e-14, 1e5)))
  expect_equal(coef(a), coef(ref), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(ref)), tolerance = 1e-8)
  skip_if_missing("logbin")
  # logbin evaluates its formula in the caller, so pass the formula object itself
  lb <- suppressWarnings(do.call(logbin::logbin, list(f, data = bw, method = "glm2",
    start = st - c(0.05, rep(0, 7)), control = list(maxit = 1e5, epsilon = 1e-12))))
  expect_true(lb$boundary)
  expect_equal(coef(a), coef(lb), tolerance = 1e-3)
  expect_equal(as.numeric(logLik(a)), lb$loglik, tolerance = 1e-6)
})

test_that("log-binomial interior maximum near the boundary converges and maxit is honoured", {
  skip_if_missing("MASS")
  bw <- MASS::birthwt; bw$race <- factor(bw$race); bw$id <- seq_len(nrow(bw))
  f <- low ~ age + lwt + race + smoke + ptl + ht + ui + ftv
  expect_warning(m <- panglm(f, bw, index = "id", family = binomialfamily("log"), maxit = 5),
                 "did not converge in 5 Newton iterations")
  expect_identical(m$iterations, 5L)
  expect_no_warning(a <- panglm(f, bw, index = "id", family = binomialfamily("log")))
  expect_true(a$converged); expect_lt(max(fitted(a)), 0.999)
  X <- a$X; mu <- fitted(a)
  expect_lt(max(abs(crossprod(X, (bw$low - mu) / (1 - mu)))), 1e-6)
  # glm() stops short from its default start; from the panglm estimate it stays there
  st <- c(log(mean(bw$low)), rep(0, 9))
  g0 <- suppressWarnings(glm(f, binomial("log"), bw, start = st, control = glm.control(1e-14, 1e5)))
  expect_gt(as.numeric(logLik(a)), as.numeric(logLik(g0)))
  ref <- glm(f, binomial("log"), bw, start = coef(a), control = glm.control(1e-14, 1e5))
  expect_equal(coef(a), coef(ref), tolerance = 1e-8)
  expect_equal(vcov(a), vcov(ref), tolerance = 1e-6)
  skip_if_missing("logbin")
  # logbin evaluates its formula in the caller, so pass the formula object itself
  lb <- suppressWarnings(do.call(logbin::logbin, list(f, data = bw, method = "glm2",
    start = st - c(0.05, rep(0, 9)), control = list(maxit = 1e5, epsilon = 1e-12))))
  expect_false(lb$boundary)
  expect_equal(coef(a), coef(lb), tolerance = 1e-3)
  expect_equal(as.numeric(logLik(a)), lb$loglik, tolerance = 1e-6)
})

test_that("the log-binomial model is pooled only", {
  data(copd)
  expect_error(panglm(hospitalized ~ crp, copd, index = c("id", "visit"), model = "random",
                      family = binomialfamily("log")), "model = 'pooling' only")
})
