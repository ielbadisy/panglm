make_dynamic_data <- function(G = 800, Tt = 5, seed = 11) {
  set.seed(seed)
  x <- matrix(rnorm(G * (Tt + 1)), G)
  a <- rnorm(G, 0, 0.8)
  y <- matrix(0L, G, Tt + 1)
  y[, 1] <- as.integer(0.3 * x[, 1] + a + rnorm(G) > 0)
  for (t in 2:(Tt + 1)) y[, t] <- as.integer(-0.5 + 0.8 * y[, t - 1] + 0.5 * x[, t] + a + rnorm(G) > 0)
  data.frame(id = rep(seq_len(G), each = Tt + 1), wave = rep(0:Tt, G), x = c(t(x)), y = c(t(y)))
}

test_that("dynamic = TRUE equals the hand-built Wooldridge specification", {
  d <- make_dynamic_data()
  m <- panglm(y ~ x, d, index = c("id", "wave"), model = "random", family = "probit",
              dynamic = TRUE)
  expect_true(m$cre)
  expect_identical(m$dynamic$lag, "lag_y")
  d2 <- d[order(d$id, d$wave), ]
  d2$lag_y <- ave(d2$y, d2$id, FUN = function(v) c(NA, utils::head(v, -1)))
  d2$initial_y <- ave(d2$y, d2$id, FUN = function(v) v[1])
  d2 <- d2[d2$wave > 0, ]
  d2$mean_x <- ave(d2$x, d2$id)
  mm <- panglm(y ~ x + lag_y + initial_y + mean_x, d2, index = c("id", "wave"),
               model = "random", family = "probit")
  expect_equal(unname(coef(m)), unname(coef(mm)), tolerance = 1e-6)
  ape <- partialeffects(m, variables = "lag_y")
  expect_identical(ape$contrast, "1 - 0")
})

test_that("dynamic models skip gaps and handle factor outcomes", {
  d <- make_dynamic_data(G = 300)
  d <- d[!(d$id %% 3 == 0 & d$wave == 2), ]
  m <- panglm(y ~ x, d, index = c("id", "wave"), model = "pooling", family = "binomial",
              dynamic = TRUE)
  # rows of wave 0, and wave 3 for individuals missing wave 2, have no lag
  expect_equal(m$dynamic$dropped, 300 + 100)
  d$o <- factor(cut(d$x + d$y, c(-Inf, -0.5, 0.5, Inf), labels = c("a", "b", "c")),
                ordered = TRUE)
  mo <- panglm(o ~ x, d, index = c("id", "wave"), model = "pooling", family = "ologit",
               dynamic = TRUE)
  expect_true(all(c("lag_ob", "lag_oc", "initial_ob", "initial_oc") %in% names(coef(mo))))
  expect_false("mean_lag_ob" %in% names(coef(mo)))
})

test_that("dynamic input checks", {
  d <- make_dynamic_data(G = 50)
  expect_error(panglm(y ~ x, d, index = c("id", "wave"), model = "within", dynamic = TRUE),
               "Nickell")
  expect_error(panglm(y ~ x, d, index = "id", dynamic = TRUE), "index")
  d$wave <- as.character(d$wave)
  expect_error(panglm(y ~ x, d, index = c("id", "wave"), dynamic = TRUE), "numeric time")
})
