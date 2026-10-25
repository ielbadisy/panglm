sel_sim <- function(N = 200, TT = 4, rho = 0.7, attrition = FALSE, seed = 1) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
  d$x <- rnorm(N * TT) + rep(rnorm(N), each = TT)
  d$z <- rnorm(N * TT)
  d$w <- rep(rbinom(N, 1, 0.5), each = TT)
  xbar <- ave(d$x, d$id)
  alpha <- 0.8 * xbar + rep(rnorm(N, sd = 0.3), each = TT)
  v <- rnorm(N * TT)
  u <- rho * v + sqrt(1 - rho^2) * rnorm(N * TT)
  d$s <- as.numeric(0.3 + 0.5 * d$x + d$z + 0.4 * xbar - 0.3 * d$w + v > 0)
  if (attrition) {
    d$s[d$t == 1] <- 1
    d$s <- ave(d$s, d$id, FUN = cumprod)
  }
  d$y <- ifelse(d$s == 1, 1 + 0.5 * d$x + 0.3 * d$w + alpha + u, NA)
  d
}

test_that("per-period probits and inverse Mills ratios match glm by hand", {
  d <- sel_sim()
  fit <- panselect(s ~ x + z + w, y ~ x + w, data = d, index = c("id", "t"), boot = 0)
  d$mean_x <- ave(d$x, d$id); d$mean_z <- ave(d$z, d$id)
  for (p in 1:4) {
    r <- d$t == p
    g <- glm(s ~ x + z + w + mean_x + mean_z, data = d[r, ],
             family = binomial(link = "probit"),
             control = glm.control(epsilon = 1e-12, maxit = 100))
    expect_equal(unname(fit$probit[[as.character(p)]]), unname(coef(g)), tolerance = 1e-7)
    xb <- predict(g, type = "link")
    expect_equal(fit$imr[r], unname(dnorm(xb) / pnorm(xb)), tolerance = 1e-7)
  }
  expect_identical(fit$lambda_terms, paste0("lambda:", 1:4))
  expect_identical(fit$vcov_type, "cluster")
})

test_that("the within second step equals plm on the hand-built design", {
  skip_if_missing("plm")
  d <- sel_sim(attrition = TRUE, seed = 2)
  fit <- panselect(s ~ x + z + w, y ~ x, data = d, index = c("id", "t"),
                   model = "within", boot = 0)
  # wave 1: everyone observed, no probit and no correction term
  expect_identical(fit$lambda_terms, paste0("lambda:", 2:4))
  expect_null(fit$probit[["1"]])
  for (p in 2:4) d[[paste0("lam", p)]] <- fit$imr * (d$t == p)
  ref <- plm::plm(y ~ x + lam2 + lam3 + lam4, data = d[d$s == 1, ],
                  index = c("id", "t"), model = "within")
  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-10)
  expect_equal(unname(vcov(fit, type = "classical")), unname(vcov(ref)), tolerance = 1e-10)
  vc <- plm::vcovHC(ref, method = "arellano", type = "sss", cluster = "group")
  expect_equal(unname(vcov(fit, type = "cluster")), unname(as.matrix(unclass(vc)[, ])), tolerance = 1e-10)
  b <- coef(ref)[2:4]
  w <- as.numeric(t(b) %*% solve(vc[2:4, 2:4], b))
  expect_equal(fit$selection_test$cluster$statistic, w, tolerance = 1e-10)
  expect_equal(fit$selection_test$cluster$df, 3)
})

test_that("the pooled second step with a common lambda equals lm + vcovCL", {
  skip_if_missing("sandwich")
  d <- sel_sim(seed = 3)
  fit <- panselect(s ~ x + z + w, y ~ x + w, data = d, index = c("id", "t"),
                   model = "pooling", lambda = "common", boot = 0)
  d$lambda <- fit$imr
  d$mean_x <- ave(d$x, d$id); d$mean_z <- ave(d$z, d$id)
  ref <- lm(y ~ x + w + mean_x + mean_z + lambda, data = d[d$s == 1, ])
  expect_equal(coef(fit), coef(ref), tolerance = 1e-10)
  expect_equal(vcov(fit, type = "cluster"),
               sandwich::vcovCL(ref, cluster = ~id, type = "HC1"), tolerance = 1e-10)
})

test_that("one period reproduces the two-step heckit", {
  skip_if_missing("sampleSelection")
  data("Mroz87", package = "sampleSelection", envir = environment())
  Mroz87$kids <- as.numeric(Mroz87$kids5 + Mroz87$kids618 > 0)
  Mroz87$id <- seq_len(nrow(Mroz87)); Mroz87$t <- 1
  hk <- sampleSelection::heckit(lfp ~ age + I(age^2) + faminc + kids + educ,
                                wage ~ exper + I(exper^2) + educ + city, data = Mroz87)
  fit <- panselect(lfp ~ age + I(age^2) + faminc + kids + educ,
                   wage ~ exper + I(exper^2) + educ + city, data = Mroz87,
                   index = c("id", "t"), model = "pooling", boot = 0)
  expect_equal(unname(coef(fit)), unname(coef(hk)[7:12]), tolerance = 1e-6)
})

test_that("the panel bootstrap is reproducible and leaves the RNG state alone", {
  d <- sel_sim(N = 120, seed = 4)
  set.seed(10); a <- runif(1); set.seed(10)
  f1 <- panselect(s ~ x + z + w, y ~ x, data = d, index = c("id", "t"), boot = 10, seed = 1)
  expect_identical(runif(1), a)
  f2 <- panselect(s ~ x + z + w, y ~ x, data = d, index = c("id", "t"), boot = 10, seed = 1)
  expect_identical(vcov(f1), vcov(f2))
  expect_identical(f1$vcov_type, "boot")
  expect_equal(vcov(f1, type = "boot"), vcov(f1))
  expect_equal(dim(f1$boot$draws), c(10L, length(coef(f1))))
  expect_false(is.null(f1$selection_test$boot))
  expect_equal(nobs(f1), sum(d$s))
  ci <- confint(f1)
  expect_equal(unname(ci[, 1]), unname(coef(f1) - qnorm(0.975) * sqrt(diag(vcov(f1)))))
  expect_output(print(f1), "Wooldridge selection correction")
  expect_output(print(summary(f1)), "Test of selection bias")
  expect_output(print(summary(f1)), "panel bootstrap, 10 replicates")
})

test_that("panselect corrects the bias of the naive within estimator", {
  skip_on_cran()
  skip_if_missing("plm")
  d <- sel_sim(N = 3000, TT = 4, rho = 0.9, seed = 5)
  naive <- coef(plm::plm(y ~ x, d[d$s == 1, ], index = c("id", "t"), model = "within"))[["x"]]
  fit <- panselect(s ~ x + z + w, y ~ x, data = d, index = c("id", "t"), boot = 0)
  expect_gt(abs(naive - 0.5), 0.07)
  expect_lt(abs(coef(fit)[["x"]] - 0.5), 3 * sqrt(vcov(fit)["x", "x"]))
  expect_lt(fit$selection_test$cluster$p.value, 0.001)
})

test_that("panselect() reports invalid input clearly", {
  d <- sel_sim(N = 50)
  idx <- c("id", "t")
  expect_error(panselect(s ~ x + z, y ~ x, data = d), "'index' is required")
  expect_error(panselect(s ~ x + z, y ~ x, data = d, index = "id"), "two distinct")
  expect_error(panselect(s ~ x + z, y ~ x, data = d, index = c("id", "nope")), "not found")
  expect_error(panselect(~ x + z, y ~ x, data = d, index = idx), "two-sided")
  expect_error(panselect(s ~ x + z, y ~ x, data = d, index = idx, boot = 1), "'boot'")
  d$s2 <- d$s * 2
  expect_error(panselect(s2 ~ x + z, y ~ x, data = d, index = idx, boot = 0), "0/1")
  expect_warning(panselect(s ~ x, y ~ x, data = d, index = idx, boot = 0), "nonlinearity")
  d$all <- 1
  expect_error(panselect(all ~ x + z, y ~ x, data = d, index = idx, boot = 0),
               "no selection to correct")
  expect_error(vcov(panselect(s ~ x + z, y ~ x, data = d, index = idx, boot = 0), type = "boot"),
               "boot = 0")
})
