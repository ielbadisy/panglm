qr_sim <- function(N = 150, TT = 8, seed = 1) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
  a <- rep(rnorm(N), each = TT)
  d$x <- rnorm(N * TT) + 0.5 * a
  d$z <- runif(N * TT, 0, 2)
  d$y <- a + d$x + 0.5 * d$z + (1 + 0.5 * d$z) * rnorm(N * TT)
  d
}

test_that("Canay estimates equal rq() on the hand-adjusted outcome", {
  skip_if_missing("quantreg")
  d <- qr_sim()
  a <- panqr(y ~ x + z, d, index = c("id", "t"), tau = c(0.25, 0.9), model = "canay",
             vcov = "none")
  fe <- lm(y ~ x + z + factor(id), data = d)
  b <- coef(fe)[c("x", "z")]
  alpha <- d$y - residuals(fe) - as.vector(as.matrix(d[, c("x", "z")]) %*% b)
  d$yadj <- d$y - (alpha - mean(alpha))
  expect_equal(unname(a$within), unname(b), tolerance = 1e-8)
  for (j in 1:2) {
    r <- quantreg::rq(yadj ~ x + z, tau = c(0.25, 0.9)[j], data = d)
    expect_equal(a$coefficients[, j], coef(r), tolerance = 1e-8)
  }
  expect_null(a$vcov)
  expect_error(vcov(a), "no covariance")
})

test_that("pooled fit equals rq() and its kernel sandwich", {
  skip_if_missing("quantreg")
  d <- qr_sim(seed = 2)
  p <- panqr(y ~ x + z, d, index = c("id", "t"), tau = 0.7, model = "pooling")
  r <- quantreg::rq(y ~ x + z, tau = 0.7, data = d)
  expect_equal(coef(p), coef(r), tolerance = 1e-10)
  s <- summary(r, se = "ker", covariance = TRUE)
  X <- model.matrix(~ x + z, d)
  meat <- crossprod(rowsum(X * (0.7 - (residuals(r) < 0)), d$id))
  expect_equal(unname(vcov(p)), unname(s$Hinv %*% meat %*% s$Hinv), tolerance = 1e-8)
})

test_that("MM-QR equals a hand implementation from lm() and rq()", {
  skip_if_missing("quantreg")
  # n = 255: n * tau is not an integer, so the sample quantile is unique
  d <- qr_sim(N = 51, TT = 5, seed = 3)
  taus <- c(0.1, 0.5, 0.9)
  m <- panqr(y ~ x + z, d, index = c("id", "t"), tau = taus, model = "mmqr")
  l1 <- lm(y ~ x + z + factor(id), data = d)
  u <- residuals(l1)
  au <- 2 * u * ((u >= 0) - mean(u >= 0))
  l2 <- lm(au ~ x + z + factor(id), data = d)
  us <- u / fitted(l2)
  b <- coef(l1)[c("x", "z")]
  g <- coef(l2)[c("x", "z")]
  expect_equal(m$location$coefficients, b, tolerance = 1e-8)
  expect_equal(m$scale$coefficients, g, tolerance = 1e-8)
  for (j in seq_along(taus)) {
    q <- unname(coef(quantreg::rq(us ~ 1, tau = taus[j])))
    expect_equal(unname(m$q[j]), q, tolerance = 1e-10)
    expect_equal(m$coefficients[, j], b + g * q, tolerance = 1e-8)
  }
  expect_true(all(sapply(m$vcov, function(v) all(diag(v) > 0))))
  expect_equal(dim(coef(m)), c(2L, 3L))
  expect_identical(colnames(coef(m)), c("tau=0.1", "tau=0.5", "tau=0.9"))
})

test_that("MM-QR location and scale covariances match plm cluster-robust", {
  skip_if_missing("quantreg")
  skip_if_missing("plm")
  d <- qr_sim(N = 80, TT = 5, seed = 4)
  m <- panqr(y ~ x + z, d, index = c("id", "t"), tau = 0.5, model = "mmqr")
  w <- plm::plm(y ~ x + z, d, index = c("id", "t"), model = "within")
  expect_equal(m$location$vcov,
               unclass(plm::vcovHC(w, method = "arellano", type = "sss", cluster = "group"))[
                 c("x", "z"), c("x", "z")], tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("known beta(tau) is recovered in a location-scale panel", {
  skip_if_missing("quantreg")
  skip_on_cran()
  d <- qr_sim(N = 300, TT = 40, seed = 5)
  tau <- c(0.25, 0.75)
  truth <- rbind(x = c(1, 1), z = 0.5 + 0.5 * qnorm(tau))
  m <- panqr(y ~ x + z, d, index = c("id", "t"), tau = tau, model = "mmqr", vcov = "none")
  cn <- panqr(y ~ x + z, d, index = c("id", "t"), tau = tau, model = "canay", vcov = "none")
  expect_lt(max(abs(m$coefficients - truth)), 0.08)
  expect_lt(max(abs(cn$coefficients[c("x", "z"), ] - truth)), 0.08)
})

test_that("cluster bootstrap is reproducible and close to the analytic SE", {
  skip_if_missing("quantreg")
  d <- qr_sim(N = 100, TT = 8, seed = 6)
  set.seed(11)
  b1 <- panqr(y ~ x + z, d, index = c("id", "t"), tau = 0.5, model = "canay", boot = 40)
  set.seed(11)
  b2 <- panqr(y ~ x + z, d, index = c("id", "t"), tau = 0.5, model = "canay", boot = 40)
  expect_identical(vcov(b1), vcov(b2))
  expect_identical(b1$vcov_type, "bootstrap")
  expect_equal(dim(b1$boot), c(40L, 3L, 1L))
  an <- panqr(y ~ x + z, d, index = c("id", "t"), tau = 0.5, model = "canay",
              vcov = "analytic")
  ratio <- sqrt(diag(vcov(b1))) / sqrt(diag(vcov(an)))
  expect_true(all(ratio > 0.5 & ratio < 2))
  set.seed(12)
  mb <- panqr(y ~ x + z, d, index = c("id", "t"), tau = c(0.25, 0.75), model = "mmqr",
              vcov = "bootstrap", boot = 30)
  expect_length(vcov(mb), 2L)
  expect_equal(dim(vcov(mb, tau = 0.75)), c(2L, 2L))
})

test_that("methods, singletons, time-invariant regressors and input checks", {
  skip_if_missing("quantreg")
  d <- qr_sim(N = 60, TT = 5, seed = 7)
  d$female <- rep(rbinom(60, 1, 0.5), each = 5)
  expect_warning(
    m <- panqr(y ~ x + z + female, d, index = c("id", "t"), tau = c(0.5, 0.9),
               model = "mmqr"),
    "time-invariant")
  expect_identical(m$dropped, "female")
  expect_identical(rownames(m$coefficients), c("x", "z"))
  expect_equal(nobs(m), 300)
  ci <- confint(m, tau = 0.9)
  expect_equal(dim(ci), c(2L, 2L))
  expect_length(confint(m), 2L)
  expect_output(print(m), "MM-QR")
  expect_output(print(summary(m)), "Scale")
  expect_s3_class(summary(m), "summary.panqr")
  expect_error(vcov(m, tau = 0.3), "fitted quantiles")

  d2 <- rbind(d, data.frame(id = 999, t = 1, x = 0, z = 1, y = 2, female = 0))
  expect_message(
    m2 <- suppressWarnings(panqr(y ~ x + z, d2, index = c("id", "t"), tau = 0.5,
                                 model = "canay", vcov = "none")),
    "observed only once")
  m1 <- panqr(y ~ x + z, d, index = c("id", "t"), tau = 0.5, model = "canay", vcov = "none")
  expect_equal(coef(m2), coef(m1))
  expect_equal(m2$n_singletons, 1L)

  expect_error(panqr(y ~ x, d, index = c("id", "t"), tau = 1), "strictly between")
  expect_error(panqr(y ~ x, d, tau = 0.5), "index")
  expect_error(panqr(y ~ x - 1, d, index = c("id", "t"), model = "mmqr"), "intercept")
  expect_error(suppressWarnings(panqr(y ~ female, d, index = c("id", "t"), model = "mmqr",
                                      vcov = "none")), "no time-varying")
})

test_that("panqr() takes boot and seed; R is a deprecated alias of boot", {
  skip_if_missing("quantreg")
  d <- qr_sim(N = 60, TT = 6, seed = 6)
  set.seed(11)
  a <- panqr(y ~ x + z, d, index = c("id", "t"), model = "canay", boot = 20)
  set.seed(11)
  panglm:::deprecate_reset()
  expect_warning(
    b <- panqr(y ~ x + z, d, index = c("id", "t"), model = "canay", R = 20),
    "boot", class = "deprecatedWarning"
  )
  expect_identical(vcov(a), vcov(b))
  expect_identical(a$R, 20)
  expect_error(panqr(y ~ x + z, d, index = c("id", "t"), boot = 20, R = 20),
               "only one of")

  # seed: reproducible draws, global stream left untouched
  set.seed(5); before <- runif(1); set.seed(5)
  s1 <- panqr(y ~ x + z, d, index = c("id", "t"), model = "canay", boot = 20, seed = 3)
  expect_identical(runif(1), before)
  s2 <- panqr(y ~ x + z, d, index = c("id", "t"), model = "canay", boot = 20, seed = 3)
  expect_identical(vcov(s1), vcov(s2))
  set.seed(3)
  s3 <- panqr(y ~ x + z, d, index = c("id", "t"), model = "canay", boot = 20)
  expect_identical(vcov(s1), vcov(s3))
})

test_that("MM-QR analytic covariance is the cluster sandwich of the influence function", {
  d <- qr_sim(N = 60, TT = 5, seed = 8)
  taus <- c(0.25, 0.9)
  m <- panqr(y ~ x + z, d, index = c("id", "t"), tau = taus, model = "mmqr")
  mi <- panqr(y ~ x + z, d, index = c("id", "t"), tau = taus, model = "mmqr", vcov = "iid")
  expect_identical(m$vcov_type, "analytic")
  expect_identical(m$n_clusters, 60L)
  # hand version from lm() with individual dummies
  l1 <- lm(y ~ x + z + factor(id), data = d)
  u <- unname(residuals(l1))
  au <- 2 * u * ((u >= 0) - mean(u >= 0))
  l2 <- lm(au ~ x + z + factor(id), data = d)
  sigma <- unname(fitted(l2))
  us <- u / sigma
  g <- coef(l2)[c("x", "z")]
  Xw <- as.matrix(d[, c("x", "z")]) - apply(d[, c("x", "z")], 2, ave, d$id)
  si <- sigma - as.vector(sweep(as.matrix(d[, c("x", "z")]), 2,
                                colMeans(d[, c("x", "z")])) %*% g)
  Qinv <- solve(crossprod(Xw))
  n <- nrow(d)
  for (j in seq_along(taus)) {
    tau <- taus[j]
    # q(tau) from the hand residuals: the indicator I(U <= q) at the
    # quantile itself would flip on a rounding difference
    sq <- sort(us)
    q <- sq[ceiling(n * tau - 1e-10)]
    expect_equal(q, unname(m$q[j]), tolerance = 1e-10)
    av <- 2 * us * ((us >= 0) - mean(us >= 0)) - 1
    h <- panglm:::qr_bandwidth_hs(tau, n)
    while ((tau - h <= 0) || (tau + h >= 1)) h <- h / 2
    fq <- 2 * h / (sq[ceiling(n * (tau + h) - 1e-10)] - sq[ceiling(n * (tau - h) - 1e-10)])
    w <- (tau - (us <= q)) / fq - us - q * av
    psi <- Xw * (u + q * sigma * av)
    psi <- psi %*% Qinv + outer(si * w / sum(si), g)
    s <- rowsum(psi, d$id)
    V <- crossprod(s) * 60 / 59
    expect_equal(unname(m$vcov[[j]]), unname(V), tolerance = 1e-8)
  }
  # the iid formula stays available and differs
  expect_false(isTRUE(all.equal(m$vcov[[1]], mi$vcov[[1]])))
  # clustering on the individual column reproduces the default
  mc <- panqr(y ~ x + z, d, index = c("id", "t"), tau = taus, model = "mmqr", cluster = "id")
  expect_equal(mc$vcov, m$vcov)
  d$obs <- seq_len(n)
  mo <- panqr(y ~ x + z, d, index = c("id", "t"), tau = taus, model = "mmqr", cluster = "obs")
  expect_identical(mo$n_clusters, n)
  expect_error(panqr(y ~ x + z, d, index = c("id", "t"), model = "canay", vcov = "iid"),
               "mmqr")
  expect_error(panqr(y ~ x + z, d, index = c("id", "t"), model = "mmqr", cluster = "nope"),
               "column")
  expect_error(panqr(y ~ x + z, d, index = c("id", "t"), model = "mmqr", cluster = "id",
                     vcov = "bootstrap", boot = 5), "analytic")
  expect_output(print(summary(m)), "cluster by individual, G = 60")
})
