test_that("pooled ordered logit and probit match MASS::polr", {
  skip_if_not_installed("MASS")
  d <- make_ordered_data()
  p <- panglm(y ~ x1 + z, d, c("id", "t"), model = "pooling", family = "ologit")
  m <- MASS::polr(y ~ x1 + z, d, Hess = TRUE)
  expect_equal(unname(coef(p)), unname(c(coef(m), m$zeta)), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(p)))), unname(sqrt(diag(vcov(m)))), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(p)), as.numeric(logLik(m)), tolerance = 1e-8)
  expect_s3_class(p, "panordered")

  pw <- panglm(y ~ x1 + z, d, c("id", "t"), model = "pooling",
               family = orderedfamily("probit"), weights = d$w)
  mw <- suppressWarnings(MASS::polr(y ~ x1 + z, d, weights = w, method = "probit"))
  expect_equal(unname(coef(pw)), unname(c(coef(mw), mw$zeta)), tolerance = 1e-5)
})

test_that("pooled ordered robust covariances are sandwiches of the scores", {
  d <- make_ordered_data(G = 120)
  p <- panglm(y ~ x1 + z, d, c("id", "t"), model = "pooling", family = "oprobit")
  vc <- vcov(p, type = "cluster")
  vh <- vcov(p, type = "HC1")
  expect_equal(dim(vc), c(5L, 5L))
  expect_true(all(diag(vc) > 0) && all(diag(vh) > 0))
  pc <- panglm(y ~ x1 + z, d, c("id", "t"), model = "pooling", family = "oprobit",
               vcov = "cluster")
  expect_equal(vcov(pc), vc)
  expect_equal(vcov(pc, type = "classical"), vcov(p))
})

test_that("random-effects ordered logit matches ordinal::clmm", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- make_ordered_data()
  r <- panglm(y ~ x1 + z, d, c("id", "t"), model = "random", family = "ologit", R = 40)
  cm <- ordinal::clmm(y ~ x1 + z + (1 | id), data = transform(d, id = factor(id)), nAGQ = 30)
  expect_equal(unname(coef(r)), unname(c(cm$beta, cm$alpha)), tolerance = 1e-4)
  expect_equal(r$sigma, unname(sqrt(unlist(ordinal::VarCorr(cm)))), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(r)), as.numeric(logLik(cm)), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(r)))), unname(sqrt(diag(vcov(cm)))[c(4, 5, 1:3)]),
               tolerance = 1e-3)
})

test_that("random-effects ordered models with R = 1 match the Laplace fit of ordinal::clmm", {
  skip_on_cran()
  skip_if_not_installed("ordinal")
  d <- make_ordered_data()
  dd <- transform(d, id = factor(id))
  for (lk in c("logit", "probit")) {
    r <- panglm(y ~ x1 + z, d, c("id", "t"), model = "random", family = orderedfamily(lk), R = 1)
    # clmm's inner Laplace mode search warns on a few individuals at this
    # tight tolerance; its fit is still the reference optimum
    cm <- suppressWarnings(ordinal::clmm(y ~ x1 + z + (1 | id), data = dd, link = lk, nAGQ = 1,
                                         control = ordinal::clmm.control(gradTol = 1e-10)))
    expect_equal(unname(coef(r)), unname(c(cm$beta, cm$alpha)), tolerance = 1e-4)
    expect_equal(as.numeric(logLik(r)), as.numeric(logLik(cm)), tolerance = 1e-7)
  }
})

test_that("random-effects ordered predictions, weights and cluster vcov", {
  d <- make_ordered_data(G = 150)
  r <- panglm(y ~ x1 + z, d, c("id", "t"), model = "random", family = "oprobit")
  pr <- predict(r)
  expect_equal(dim(pr), c(nrow(d), 4L))
  expect_equal(unname(rowSums(pr)), rep(1, nrow(d)), tolerance = 1e-10)
  expect_equal(unname(predict(r, newdata = d)), unname(pr), tolerance = 1e-10)
  expect_s3_class(predict(r, newdata = d[1:3, ], type = "class"), "factor")
  expect_true(all(diag(vcov(r, type = "cluster")) > 0))
  expect_error(vcov(r, type = "HC1"), "not defined")

  # integer weights equal duplicated individuals
  w <- rep(rep(1:2, length.out = 150), each = 5)
  rw <- panglm(y ~ x1 + z, d, c("id", "t"), model = "random", family = "oprobit", weights = w)
  dup <- d[rep(seq_len(nrow(d)), w), ]
  dup$id <- dup$id + 1000 * (duplicated(dup[, c("id", "t")]))
  rd <- panglm(y ~ x1 + z, dup, c("id", "t"), model = "random", family = "oprobit")
  expect_equal(coef(rw), coef(rd), tolerance = 1e-6)
  expect_error(panglm(y ~ x1, d, c("id", "t"), model = "random", family = "oprobit",
                      weights = seq_len(nrow(d))), "constant within individuals")
})

test_that("fixed-effects ordered logit is the blow-up and cluster estimator", {
  skip_if_not_installed("survival")
  d <- make_ordered_data()
  f <- panglm(y ~ x1, d, c("id", "t"), model = "within", family = "ologit")
  st <- do.call(rbind, lapply(2:4, function(k) {
    dd <- d; dd$yb <- as.numeric(as.integer(d$y) >= k); dd$strat <- paste(d$id, k); dd
  }))
  st$time <- 1
  cl <- survival::coxph(survival::Surv(time, yb) ~ x1 + survival::strata(strat), data = st,
                        method = "exact")
  expect_equal(unname(coef(f)), unname(coef(cl)), tolerance = 1e-6)
  expect_identical(f$vcov_type, "cluster")
  expect_error(panglm(y ~ x1, d, c("id", "t"), model = "within", family = "oprobit"),
               "logit link only")
  expect_error(predict(f, newdata = d), "not identified")
})

test_that("ordered partial effects sum to zero over categories", {
  d <- make_ordered_data(G = 150)
  r <- panglm(y ~ x1 + z, d, c("id", "t"), model = "random", family = "oprobit", cre = TRUE)
  ape <- partialeffects(r)
  expect_equal(nrow(ape), 8L)
  sums <- tapply(ape$estimate, ape$term, sum)
  expect_equal(as.numeric(sums), c(0, 0), tolerance = 1e-8)
  # numeric response with integer codes gives the same fit as the factor
  d$yn <- as.integer(d$y)
  p1 <- panglm(y ~ x1 + z, d, c("id", "t"), family = "ologit")
  p2 <- panglm(yn ~ x1 + z, d, c("id", "t"), family = "ologit")
  expect_equal(unname(coef(p1)), unname(coef(p2)))
})
