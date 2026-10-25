skip_if_not_installed("splines")

spline_panel <- function(seed, f, N = 200, T = 5) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(N), each = T), t = rep(seq_len(T), N))
  a <- rnorm(N)
  d$z <- rnorm(nrow(d))
  d$x <- runif(nrow(d), 0, 4) + 0.5 * a[d$id]
  d$y <- a[d$id] + f(d$x) + 0.3 * d$z + rnorm(nrow(d), sd = 0.3)
  d[sample(nrow(d)), ]
}

test_that("fixed effects with ns() recovers a nonlinear curve and predict agrees", {
  f <- function(x) sin(1.5 * x)
  d <- spline_panel(1, f)
  fit <- panglm(y ~ splines::ns(x, 5) + z, data = d, index = c("id", "t"),
                model = "within", vcov = "cluster")
  ec <- effectcurve(fit, "x", values = seq(0.5, 3.5, length.out = 20), ref = 2)
  expect_s3_class(ec, "panglm_effectcurve")
  expect_lt(max(abs(ec$estimate - (f(ec$value) - f(2)))), 0.15)
  nd <- d[rep(1, 20), ]; nd$x <- ec$value
  nd0 <- d[1, ]; nd0$x <- 2
  expect_equal(unname(predict(fit, nd) - predict(fit, nd0)), ec$estimate, tolerance = 1e-10)
  # model.matrix: in-sample equals the stored design, new data rebuilds the basis
  expect_equal(colnames(model.matrix(fit)), names(coef(fit)))
  expect_equal(nrow(model.matrix(fit)), nrow(d))
  expect_equal(unname(model.matrix(fit, d)), unname(model.matrix(fit)), tolerance = 1e-12)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(ec))
})

test_that("lineartest rejects a curve and keeps its size under linearity", {
  d <- spline_panel(2, function(x) sin(1.5 * x))
  fit <- panglm(y ~ splines::ns(x, 4) + z, data = d, index = c("id", "t"),
                model = "within", vcov = "cluster")
  lt <- lineartest(fit, "x")
  expect_equal(lt$df, 3)
  expect_lt(lt$p.value, 1e-6)
  p <- vapply(1:60, function(s) {
    d <- spline_panel(100 + s, function(x) 0.5 * x)
    lineartest(panglm(y ~ splines::ns(x, 4) + z, data = d, index = c("id", "t"),
                      model = "within"), "x")$p.value
  }, numeric(1))
  expect_lt(mean(p < 0.05), 0.15)
})

test_that("FE likelihood-ratio statistic is the Wald statistic inflated by T / (T - 1)", {
  # the ML error variance of the dummy-variable fit is RSS / (NT), not RSS / (N(T - 1))
  d <- spline_panel(3, function(x) 0.4 * x + 0.02 * (x - 2)^2, N = 400)
  f1 <- panglm(y ~ splines::ns(x, 3) + z, data = d, index = c("id", "t"), model = "within")
  f0 <- panglm(y ~ x + z, data = d, index = c("id", "t"), model = "within")
  lr <- 2 * as.numeric(logLik(f1) - logLik(f0))
  expect_equal(lineartest(f1, "x")$statistic * 5 / 4, lr, tolerance = 0.03)
})

test_that("jointtest matches a hand computation and R-matrix form", {
  d <- spline_panel(4, function(x) cos(x))
  fit <- panglm(y ~ splines::ns(x, 3) + z, data = d, index = c("id", "t"), model = "within")
  nm <- grep("ns", names(coef(fit)), value = TRUE)
  b <- coef(fit)[nm]; V <- vcov(fit)[nm, nm]
  jt <- jointtest(fit, terms = nm)
  expect_equal(jt$statistic, as.numeric(t(b) %*% solve(V, b)))
  R <- diag(length(coef(fit)))[match(nm, names(coef(fit))), ]
  expect_equal(jointtest(fit, R = R)$statistic, jt$statistic)
  expect_equal(jointtest(fit, pattern = "^splines")$statistic, jt$statistic)
  expect_false(isTRUE(all.equal(jointtest(fit, terms = nm, vcov = "cluster")$statistic,
                                jt$statistic)))
  expect_error(jointtest(fit), "give")
  expect_error(jointtest(fit, terms = "nope"), "unknown")
})

test_that("spline helpers work for Poisson FE, CRE random effects and refuse bad input", {
  d <- spline_panel(5, function(x) sin(x))
  d$cnt <- rpois(nrow(d), exp(0.5 * sin(d$x)))
  fp <- panglm(cnt ~ splines::ns(x, 3), data = d, index = c("id", "t"),
               model = "within", family = "poisson")
  expect_true(all(is.finite(effectcurve(fp, "x")$se)))
  fr <- panglm(y ~ poly(x, 3) + z, data = d, index = c("id", "t"), model = "random", cre = TRUE)
  expect_lt(lineartest(fr, "x")$p.value, 0.01)
  expect_equal(nrow(effectcurve(fr, "x", values = 1:3)), 3)
  fl <- panglm(y ~ x + z, data = d, index = c("id", "t"), model = "within")
  expect_equal(effectcurve(fl, "x", values = 3, ref = 1)$estimate, 2 * unname(coef(fl)["x"]))
  expect_error(lineartest(fl, "x"), "single column")
  fi <- panglm(y ~ splines::ns(x, 3) * z, data = d, index = c("id", "t"), model = "within")
  expect_error(effectcurve(fi, "x"), "interaction")
  expect_error(effectcurve(fl, "nope"), "numeric regressor")
})
