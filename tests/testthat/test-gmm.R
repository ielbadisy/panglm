gmm_ab_formula <- function() {
  log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) + log(capital) +
    lag(log(output), 0:1) | lag(log(emp), 2:99)
}

gmm_bb_formula <- function() {
  log(emp) ~ lag(log(emp), 1) + lag(log(wage), 0:1) + lag(log(capital), 0:1) |
    lag(log(emp), 2:99) + lag(log(wage), 2:99) + lag(log(capital), 2:99)
}

# plm::pgmm() evaluates a call to plm() in its parent frame, so it needs
# plm() visible there without attaching plm
gmm_pgmm <- function(formula, data, index = NULL, ...) {
  plm <- plm::plm
  plm::pgmm(formula, data = data, index = index, ...)
}

gmm_sim <- function(N = 500, TT = 6, rho = 0.5, seed = 1) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
  mu <- rnorm(N)
  d$x <- rnorm(N * TT)
  d$y <- NA_real_
  for (i in seq_len(N)) {
    r <- (i - 1) * TT + seq_len(TT)
    y <- mu[i] / (1 - rho) + rnorm(1)
    for (s in seq_len(TT)) {
      y <- rho * y + 0.4 * d$x[r[s]] + mu[i] + rnorm(1)
      d$y[r[s]] <- y
    }
  }
  d
}

test_that("difference GMM matches plm::pgmm (Arellano-Bond employment equation)", {
  skip_if_missing("plm")
  data("EmplUK", package = "plm")
  f <- gmm_ab_formula()
  a <- pangmm(f, EmplUK, index = c("firm", "year"), model = "twostep")
  b <- gmm_pgmm(f, EmplUK, model = "twosteps")
  vb <- plm::vcovHC(b)
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-8)
  expect_equal(unname(sqrt(diag(vcov(a)))), unname(sqrt(diag(vb))), tolerance = 1e-8)
  expect_equal(unname(sqrt(diag(vcov(a, type = "classical")))), unname(sqrt(diag(vcov(b)))),
               tolerance = 1e-8)
  expect_equal(a$artest$statistic,
               unname(c(plm::mtest(b, 1, vb)$statistic, plm::mtest(b, 2, vb)$statistic)),
               tolerance = 1e-8)
  expect_equal(a$hansen$statistic, unname(plm::sargan(b)$statistic), tolerance = 1e-8)
  expect_equal(a$hansen$df, unname(plm::sargan(b)$parameter))
  expect_equal(a$n_instruments, ncol(b$W[[1]]))
  expect_equal(a$timedummies, as.character(1979:1984))

  # one step: robust SE, and classical SE / Sargan after scaling plm's B1 and A1
  a1 <- pangmm(f, EmplUK, index = c("firm", "year"))
  b1 <- gmm_pgmm(f, EmplUK, model = "onestep")
  N <- length(b1$W)
  expect_equal(unname(coef(a1)), unname(coef(b1)), tolerance = 1e-8)
  expect_equal(unname(vcov(a1)), unname(plm::vcovHC(b1)), tolerance = 1e-8)
  expect_equal(unname(vcov(a1, type = "classical")), unname(vcov(b1)) * N * a1$sigma2,
               tolerance = 1e-8)
  expect_equal(a1$sargan$statistic,
               unname(plm::sargan(b1, weights = "onestep")$statistic) / (N * a1$sigma2),
               tolerance = 1e-8)
})

test_that("collapsed instruments, lag limits and standard instruments match plm", {
  skip_if_missing("plm")
  data("EmplUK", package = "plm")
  f <- gmm_ab_formula()
  a <- pangmm(f, EmplUK, index = c("firm", "year"), collapse = TRUE, model = "twostep")
  b <- gmm_pgmm(f, EmplUK, collapse = TRUE, model = "twosteps")
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-8)
  expect_equal(unname(vcov(a)), unname(plm::vcovHC(b)), tolerance = 1e-8)
  expect_equal(a$n_instruments, ncol(b$W[[1]]))

  a4 <- pangmm(f, EmplUK, index = c("firm", "year"), maxlag = 3)
  f4 <- log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) + log(capital) +
    lag(log(output), 0:1) | lag(log(emp), 2:3)
  a4f <- pangmm(f4, EmplUK, index = c("firm", "year"))
  b4 <- gmm_pgmm(f4, EmplUK, model = "onestep")
  expect_equal(coef(a4), coef(a4f))
  expect_equal(unname(coef(a4)), unname(coef(b4)), tolerance = 1e-8)

  f3 <- log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) + log(capital) |
    lag(log(emp), 2:99) | lag(log(wage), 1:2) + log(capital)
  a3 <- pangmm(f3, EmplUK, index = c("firm", "year"), effect = "individual")
  b3 <- gmm_pgmm(f3, EmplUK, effect = "individual", model = "onestep")
  expect_equal(unname(coef(a3)), unname(coef(b3)), tolerance = 1e-8)
  expect_equal(unname(vcov(a3)), unname(plm::vcovHC(b3)), tolerance = 1e-8)
})

test_that("system GMM matches plm::pgmm with transformation = 'ld'", {
  skip_if_missing("plm")
  data("EmplUK", package = "plm")
  f <- gmm_bb_formula()
  a <- pangmm(f, EmplUK, index = c("firm", "year"), transformation = "ld", model = "twostep")
  b <- gmm_pgmm(f, EmplUK, transformation = "ld", model = "twosteps")
  vb <- plm::vcovHC(b)
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(a)))), unname(sqrt(diag(vb))), tolerance = 1e-6)
  expect_equal(a$artest$statistic,
               unname(c(plm::mtest(b, 1, vb)$statistic, plm::mtest(b, 2, vb)$statistic)),
               tolerance = 1e-6)
  expect_equal(a$hansen$statistic, unname(plm::sargan(b)$statistic), tolerance = 1e-6)
  expect_true("(Intercept)" %in% a$timedummies)

  # effect = "individual" keeps a constant in the level equation
  ai <- pangmm(f, EmplUK, index = c("firm", "year"), transformation = "ld",
               effect = "individual")
  expect_identical(ai$timedummies, "(Intercept)")
  expect_true(is.finite(coef(ai)[["(Intercept)"]]))
})

test_that("GMM recovers the persistence parameter in a simulated panel", {
  d <- gmm_sim(N = 1500, TT = 6)
  ab <- pangmm(y ~ lag(y, 1) + x | lag(y, 2:99), data = d, index = c("id", "t"),
               model = "twostep")
  bb <- pangmm(y ~ lag(y, 1) + x | lag(y, 2:99), data = d, index = c("id", "t"),
               transformation = "ld", model = "twostep", collapse = TRUE)
  expect_lt(abs(coef(ab)[["lag(y, 1)"]] - 0.5), 0.1)
  expect_lt(abs(coef(bb)[["lag(y, 1)"]] - 0.5), 0.1)
  expect_lt(abs(coef(bb)[["x"]] - 0.4), 0.05)
  # AR(1) in differences by construction, no AR(2)
  expect_lt(ab$artest$p.value[1], 0.01)
  expect_gt(ab$artest$p.value[2], 0.001)
  expect_gt(ab$hansen$p.value, 0.001)
  # Windmeijer SEs exceed the conventional two-step SEs here
  expect_true(all(sqrt(diag(vcov(ab))) > sqrt(diag(vcov(ab, type = "classical")))))
})

test_that("unbalanced panels with gaps and data order do not matter", {
  d <- gmm_sim(N = 300, TT = 6, seed = 2)
  set.seed(3)
  d <- d[-sample(nrow(d), 150), ]
  f <- y ~ lag(y, 1) + x | lag(y, 2:99)
  a <- pangmm(f, d, index = c("id", "t"))
  a_shuffled <- pangmm(f, d[sample(nrow(d)), ], index = c("id", "t"))
  expect_equal(coef(a), coef(a_shuffled), tolerance = 1e-10)
  skip_if_missing("plm")
  b <- gmm_pgmm(f, d, index = c("id", "t"), model = "onestep")
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-8)
  expect_equal(unname(vcov(a)), unname(plm::vcovHC(b)), tolerance = 1e-8)
})

test_that("methods return the expected structures", {
  d <- gmm_sim(N = 200, TT = 5, seed = 4)
  fit <- pangmm(y ~ lag(y, 1) + x | lag(y, 2:99), data = d, index = c("id", "t"),
                model = "twostep", vcov = "classical")
  expect_s3_class(fit, "pangmm")
  expect_identical(fit$vcov_type, "classical")
  expect_equal(vcov(fit), vcov(fit, type = "classical"))
  expect_false(isTRUE(all.equal(vcov(fit), vcov(fit, type = "robust"))))
  expect_equal(nobs(fit), fit$nobs_diff)
  expect_length(residuals(fit), 200)
  ci <- confint(fit)
  expect_equal(dim(ci), c(length(coef(fit)), 2L))
  s <- summary(fit)
  expect_s3_class(s, "summary.pangmm")
  expect_false(any(rownames(s$coefficients) %in% fit$timedummies))
  expect_true(all(fit$timedummies %in% rownames(summary(fit, time.dummies = TRUE)$coefficients)))
  expect_output(print(fit), "difference GMM")
  expect_output(print(s), "Hansen J")
  expect_output(print(s), "AR\\(2\\)")
})

test_that("input errors are reported", {
  d <- gmm_sim(N = 50, TT = 5, seed = 5)
  expect_error(pangmm(y ~ lag(y, 1) + x, d, index = c("id", "t")), "two parts")
  expect_error(pangmm(y ~ lag(y, 1) | lag(y, 2:99), d, index = "id"), "two distinct")
  expect_error(pangmm(y ~ lag(y, 1) | lag(y, 2:99), d, index = c("id", "zz")), "not found")
  expect_error(pangmm(y ~ lag(y, 1) | lag(y, 3:99), d, index = c("id", "t"), maxlag = 2),
               "maxlag")
  expect_error(pangmm(y ~ lag(y, 1) | lag(y, 2:99), d, index = c("id", "t"),
                      collapse = c(z = TRUE)), "unknown names")
  dd <- rbind(d, d[1, ])
  expect_error(pangmm(y ~ lag(y, 1) | lag(y, 2:99), dd, index = c("id", "t")), "unique")
})

test_that("pangmm warns about too many instruments and notes it in the summary", {
  set.seed(5)
  N <- 30; TT <- 10
  d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
  mu <- rep(rnorm(N), each = TT)
  d$y <- ave(mu + rnorm(N * TT), d$id, FUN = function(e) stats::filter(e, 0.5,
    method = "recursive"))
  w <- character()
  fit <- withCallingHandlers(
    pangmm(y ~ lag(y, 1) | lag(y, 2:99), data = d, index = c("id", "t"), model = "twostep"),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("too many instruments", w)))
  expect_true(any(grepl("generalized inverse", w)))
  expect_gte(fit$n_instruments, N)
  expect_output(print(summary(fit)), "biased toward OLS")
  expect_no_warning(fc <- pangmm(y ~ lag(y, 1) | lag(y, 2:99), data = d, index = c("id", "t"),
                                 collapse = TRUE))
  expect_false(any(is.nan(fc$artest$statistic)))
  skip_if_not_installed("plm")
  p <- suppressWarnings(gmm_pgmm(y ~ lag(y, 1) | lag(y, 2:99), data = d, index = c("id", "t"),
                                 effect = "twoways", model = "twosteps"))
  expect_equal(unname(coef(fit)[1]), unname(coef(p)[1]), tolerance = 1e-6)
})

test_that("pangmm() accepts vcov = 'cluster' as a synonym of 'robust'", {
  d <- gmm_sim(N = 150, TT = 6, seed = 4)
  f <- y ~ lag(y, 1) + x | lag(y, 2:99)
  a <- pangmm(f, d, index = c("id", "t"), model = "twostep")
  b <- pangmm(f, d, index = c("id", "t"), model = "twostep", vcov = "cluster")
  expect_identical(vcov(b), vcov(a))
  expect_identical(b$vcov_type, "robust")
  expect_identical(vcov(a, type = "cluster"), vcov(a, type = "robust"))
})

test_that("pangmm drops a regressor collinear with the time effects with a warning", {
  set.seed(3)
  N <- 150; TT <- 6
  d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
  d$age <- rep(runif(N, 20, 60), each = TT) + d$t  # first difference is 1
  d$x <- rnorm(nrow(d))
  mu <- rnorm(N)
  d$y <- 0
  for (i in seq_len(N)) {
    r <- (i - 1) * TT + seq_len(TT)
    y <- mu[i]
    for (s in seq_len(TT)) { y <- 0.5 * y + 0.3 * d$x[r[s]] + mu[i] + rnorm(1); d$y[r[s]] <- y }
  }
  expect_warning(
    a <- pangmm(y ~ lag(y, 1) + age + x | lag(y, 2:99), data = d, index = c("id", "t")),
    "dropping unidentified column\\(s\\): age")
  b <- pangmm(y ~ lag(y, 1) + x | lag(y, 2:99) | age + x, data = d, index = c("id", "t"))
  expect_equal(coef(a), coef(b))
  expect_equal(vcov(a), vcov(b))
  expect_false("age" %in% names(coef(a)))
})
