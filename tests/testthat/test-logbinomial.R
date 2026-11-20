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

test_that("log-binomial fitted risks stay below 1 where glm() fails", {
  skip_if_missing("MASS")
  bw <- MASS::birthwt; bw$id <- seq_len(nrow(bw))
  f <- low ~ smoke + lwt + age + race + ht + ui
  expect_warning(panglm(f, bw, index = "id", family = binomialfamily("log")),
                 "did not converge in 100 IRLS iterations; the largest fitted risk")
  a <- panglm(f, bw, index = "id", family = binomialfamily("log"), maxit = 2000)
  expect_lt(max(fitted(a)), 1)
  # the maximum is interior: the score is zero
  X <- a$X; mu <- fitted(a)
  expect_lt(max(abs(crossprod(X, (bw$low - mu) / (1 - mu)))), 1e-5)
})

test_that("the log-binomial model is pooled only", {
  data(copd)
  expect_error(panglm(hospitalized ~ crp, copd, index = c("id", "visit"), model = "random",
                      family = binomialfamily("log")), "model = 'pooling' only")
})
