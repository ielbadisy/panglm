test_that("weights and offset may be given as a column name of data", {
  data(copd, package = "panglm")
  d <- copd
  set.seed(1)
  d$w <- stats::runif(nrow(d), 0.2, 3)
  f_vec <- panglm(exacerbations ~ treatment + crp, d, c("id", "visit"),
                  family = "poisson", weights = d$w, vcov = "cluster")
  f_nam <- panglm(exacerbations ~ treatment + crp, d, c("id", "visit"),
                  family = "poisson", weights = "w", vcov = "cluster")
  expect_identical(coef(f_nam), coef(f_vec))
  expect_identical(vcov(f_nam), vcov(f_vec))
  expect_error(panglm(fev1 ~ crp, d, c("id", "visit"), weights = "nope"),
               "weights column 'nope' not found")

  d$zero <- 0
  f_off <- panglm(exacerbations ~ treatment + crp, d, c("id", "visit"),
                  family = "poisson", offset = "zero")
  expect_identical(coef(f_off),
                   coef(panglm(exacerbations ~ treatment + crp, d, c("id", "visit"),
                               family = "poisson")))

  dc <- make_cost_data()
  dc$w <- d$w
  t_vec <- twopart(cost0 ~ treatment + crp, dc, c("id", "visit"), weights = dc$w)
  t_nam <- twopart(cost0 ~ treatment + crp, dc, c("id", "visit"), weights = "w")
  expect_identical(coef(t_nam), coef(t_vec))
})

test_that("panfmm() takes weights and offset as column names", {
  set.seed(3)
  N <- 60; TT <- 4
  d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N),
                  x = rnorm(N * TT))
  cl <- rep(rbinom(N, 1, 0.5), each = TT)
  d$y <- rpois(N * TT, exp(0.2 + ifelse(cl == 1, 1, -0.5) * d$x))
  d$w <- rep(runif(N, 0.5, 2), each = TT)
  d$off <- log(rep(c(1, 2), length.out = N * TT))
  a <- panfmm(y ~ x, d, c("id", "t"), nstart = 2, seed = 1, weights = d$w, offset = d$off)
  b <- panfmm(y ~ x, d, c("id", "t"), nstart = 2, seed = 1, weights = "w", offset = "off")
  expect_identical(coef(a), coef(b))
  expect_identical(vcov(a), vcov(b))
})
