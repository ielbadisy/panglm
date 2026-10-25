sim_tpcre <- function(seed, rho = 0.6, N = 300, T = 5) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(N), each = T), t = rep(seq_len(T), N))
  d$x <- rnorm(N * T)
  S <- matrix(c(0.64, rho * 0.56, rho * 0.56, 0.49), 2)
  U <- MASS::mvrnorm(N, c(0, 0), S)
  d$any <- as.numeric(-0.2 + 0.5 * d$x + U[d$id, 1] + rnorm(N * T) > 0)
  d$cost <- ifelse(d$any == 1, exp(6 + 0.3 * d$x + U[d$id, 2] + rnorm(N * T, sd = 0.5)), 0)
  d
}

test_that("correlated two-part likelihood matches 2D numerical integration", {
  set.seed(3)
  T <- 4
  x <- rnorm(T); d <- c(1, 0, 1, 1); ly <- c(5.2, 6.4, 5.9)
  Xb <- cbind(1, x); Xp <- cbind(1, x[d == 1])
  theta <- c(-0.1, 0.4, 5.8, 0.2, log(0.9), log(0.6), log(0.5), atanh(0.5))
  f <- twopart_cre_contrib(Xb, d, rep(1L, T), Xp, ly, rep(1L, 3), "probit", R = 15)
  ll <- f(theta, FALSE)$ll
  s1 <- 0.9; s2 <- 0.6; se <- 0.5; rho <- 0.5
  g <- seq(-6, 6, length.out = 601); h <- diff(g[1:2])
  grid <- expand.grid(z1 = g, z2 = g)
  u1 <- s1 * grid$z1
  u2 <- s2 * (rho * grid$z1 + sqrt(1 - rho^2) * grid$z2)
  a <- as.numeric(Xb %*% theta[1:2]); b <- as.numeric(Xp %*% theta[3:4])
  lb <- rowSums(sapply(seq_len(T), function(t) pnorm((2 * d[t] - 1) * (a[t] + u1), log.p = TRUE)))
  lp <- rowSums(sapply(seq_along(ly), function(t) dnorm(ly[t], b[t] + u2, se, log = TRUE) - ly[t]))
  brute <- log(sum(exp(lb + lp) * dnorm(grid$z1) * dnorm(grid$z2)) * h^2)
  expect_equal(ll, brute, tolerance = 1e-6)
})

test_that("at rho = 0 the joint likelihood separates", {
  set.seed(4)
  T <- 3
  x <- rnorm(T); d <- c(1, 1, 0); ly <- c(5, 6)
  Xb <- cbind(1, x); Xp <- cbind(1, x[1:2])
  th <- c(0.1, 0.3, 5.5, 0.1, log(0.7), log(0.8), log(0.4), 0)
  f <- twopart_cre_contrib(Xb, d, rep(1L, T), Xp, ly, rep(1L, 2), "logit", R = 25)
  a <- as.numeric(Xb %*% th[1:2]); b <- as.numeric(Xp %*% th[3:4])
  bin <- log(stats::integrate(function(u) sapply(u, function(v)
    prod(plogis((2 * d - 1) * (a + v)))) * dnorm(u, 0, 0.7), -Inf, Inf, rel.tol = 1e-10)$value)
  S <- diag(0.4^2, 2) + 0.8^2
  pos <- -0.5 * (2 * log(2 * pi) + log(det(S)) +
    drop(t(ly - b) %*% solve(S, ly - b))) - sum(ly)
  expect_equal(f(th, FALSE)$ll, bin + pos, tolerance = 1e-7)
})

test_that("population means agree with Monte Carlo", {
  j <- list(sigma_binary = 0.8, sigma_positive = 0.7, sigma_e = 0.5, rho = 0.6)
  set.seed(1)
  x <- rnorm(500)
  m <- twopart_cre_means(-0.2 + 0.5 * x, 6 + 0.3 * x, j, "probit")
  U <- MASS::mvrnorm(4e5, c(0, 0), matrix(c(0.64, 0.336, 0.336, 0.49), 2))
  xx <- sample(x, 4e5, TRUE)
  mc_total <- mean(pnorm(-0.2 + 0.5 * xx + U[, 1]) * exp(6 + 0.3 * xx + U[, 2] + 0.125))
  mc_prob <- mean(pnorm(-0.2 + 0.5 * xx + U[, 1]))
  expect_equal(mean(m$total), mc_total, tolerance = 0.01)
  expect_equal(mean(m$prob), mc_prob, tolerance = 0.005)
})

test_that("correlated two-part fit recovers parameters and predicts the mean", {
  skip_on_cran()
  d <- sim_tpcre(11)
  tc <- twopart(cost ~ x, d, c("id", "t"), model = "random", positive = "lognormal",
                correlated = TRUE, vcov = "cluster")
  expect_true(tc$correlated)
  expect_equal(tc$joint$convergence, 0L)
  vp <- tc$joint$variance_params
  expect_lt(abs(vp[["rho"]] - 0.6), 0.25)
  expect_lt(abs(vp[["sigma_u1"]] - 0.8), 0.25)
  expect_lt(abs(vp[["sigma_e"]] - 0.5), 0.05)
  expect_lt(max(abs(coef(tc) - c(-0.2, 0.5, 6, 0.3))), 0.2)
  V <- vcov(tc)
  expect_equal(dim(V), c(4L, 4L))
  expect_true(all(is.finite(sqrt(diag(V)))))
  expect_equal(names(coef(tc)), rownames(V))
  # implied population mean at the true parameters
  j <- list(sigma_binary = 0.8, sigma_positive = 0.7, sigma_e = 0.5, rho = 0.6)
  truth <- mean(twopart_cre_means(-0.2 + 0.5 * d$x, 6 + 0.3 * d$x, j, "probit")$total)
  expect_lt(abs(mean(predict(tc)) / truth - 1), 0.15)
  p <- predict(tc, type = "prob"); pos <- predict(tc, type = "positive")
  expect_equal(predict(tc), p * pos)
  ts <- twopart(cost ~ x, d, c("id", "t"), model = "random", positive = "lognormal")
  expect_gt(as.numeric(logLik(tc)), as.numeric(logLik(ts)))
  pe <- partialeffects(tc)
  expect_true(is.finite(pe$estimate[1]) && pe$std.error[1] > 0)
  expect_output(print(summary(tc)), "rho")
})

test_that("correlated two-part argument checks", {
  d <- sim_tpcre(2, N = 40, T = 3)
  expect_error(twopart(cost ~ x, d, c("id", "t"), model = "pooling", positive = "lognormal",
                       correlated = TRUE), "random")
  expect_error(twopart(cost ~ x, d, c("id", "t"), model = "random", positive = "gamma",
                       correlated = TRUE), "lognormal")
  expect_error(twopart(cost ~ x, d, c("id", "t"), model = "random", positive = "lognormal",
                       correlated = TRUE, vcov = "HC1"))
})
