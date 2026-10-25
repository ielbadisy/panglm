et_panel <- function(seed, N = 240) {
  set.seed(seed)
  d <- expand.grid(id = seq_len(N), year = 2011:2018)
  coh <- sample(c(2014, 2016, NA), N, replace = TRUE)
  d$first <- coh[d$id]
  a <- rnorm(N, sd = 0.5)
  d$x <- rnorm(nrow(d))
  post <- !is.na(d$first) & d$year >= d$first
  d$y <- a[d$id] + 0.1 * (d$year - 2011) + 0.3 * d$x + 1 * post + rnorm(nrow(d))
  d$cnt <- rpois(nrow(d), exp(1 + a[d$id] + 0.05 * (d$year - 2011) + 0.3 * post))
  d$bin <- rbinom(nrow(d), 1, plogis(-0.5 + a[d$id] + 0.8 * post))
  d
}

test_that("etwfe recovers the effect and has consistent aggregations", {
  d <- et_panel(1)
  fit <- etwfe(y ~ x, d, c("id", "year"), "first")
  expect_s3_class(fit, "panglm_etwfe")
  s <- summary(fit)
  expect_equal(s$overall$estimate, 1, tolerance = 0.15)
  ev <- summary(fit, type = "event")$estimates
  expect_true(all(ev$level >= 0))
  gr <- summary(fit, type = "group")$estimates
  expect_equal(gr$level, c(2014, 2016))
  # linear model without controls: the simple ATT is the cell-size weighted
  # mean of the cell coefficients
  fit0 <- etwfe(y ~ 1, d, c("id", "year"), "first")
  s0 <- summary(fit0)
  n <- table(paste(d$first, d$year)[!is.na(d$first) & d$year >= d$first])
  cells <- fit0$cells
  w <- as.numeric(n[paste(cells$group, cells$time)])
  expect_equal(s0$overall$estimate, sum(w * cells$estimate) / sum(w), tolerance = 1e-10)
  expect_equal(nobs(fit), nrow(d))
  expect_output(print(fit), "ATT")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(fit))
})

test_that("etwfe nonlinear families and never-treated controls run", {
  d <- et_panel(2)
  p <- etwfe(cnt ~ 1, d, c("id", "year"), "first", family = "poisson")
  expect_true(summary(p)$overall$estimate > 0)
  b <- etwfe(bin ~ x, d, c("id", "year"), "first", family = "logit", control = "never")
  ev <- summary(b, type = "event")$estimates
  expect_true(any(ev$level < 0))
  expect_equal(ev$estimate[ev$level == -1], 0)
  expect_error(etwfe(y ~ 1, d, c("id", "year"), "first", family = "poisson"), "non-negative")
})

test_that("etwfe matches etwfe::etwfe and etwfe::emfx", {
  skip_if_not_installed("etwfe")
  data(mpdta, package = "did")
  mpdta$emp <- round(exp(mpdta$lemp))
  for (cg in c("notyet", "never")) {
    for (fam in c("gaussian", "poisson")) {
      f <- if (fam == "gaussian") lemp ~ lpop else emp ~ lpop
      ref <- suppressMessages(etwfe::etwfe(f, tvar = year, gvar = first.treat, data = mpdta,
                                           cgroup = cg, vcov = ~countyreal,
                                           family = if (fam == "gaussian") NULL else fam))
      fit <- etwfe(f, mpdta, c("countyreal", "year"), "first.treat", family = fam, control = cg)
      for (type in c("simple", "group")) {
        e <- etwfe::emfx(ref, type = type)
        s <- summary(fit, type = type)
        tab <- if (type == "simple") s$overall else s$estimates
        expect_equal(tab$estimate, e$estimate, tolerance = 1e-6)
        expect_equal(tab$se, e$std.error, tolerance = 1e-4)
      }
    }
  }
})

test_that("etwfe sampling weights", {
  d <- et_panel(3)
  fit <- etwfe(y ~ x, d, c("id", "year"), "first")
  d$w2 <- 3
  fw <- etwfe(y ~ x, d, c("id", "year"), "first", weights = "w2")
  expect_equal(summary(fw, type = "event")$estimates, summary(fit, type = "event")$estimates,
               tolerance = 1e-8)
  d$wv <- runif(nrow(d), 0.5, 2)
  p <- etwfe(cnt ~ 1, d, c("id", "year"), "first", family = "poisson", weights = d$wv)
  expect_true(is.finite(summary(p)$overall$estimate))
  expect_error(etwfe(y ~ 1, d, c("id", "year"), "first", weights = -d$wv), "non-negative")
  # weights varying within individuals: refused by attgt(), accepted by etwfe()
  expect_error(attgt(y ~ 1, d, c("id", "year"), "first", weights = d$wv), "constant within")
  fr <- etwfe(y ~ 1, d, c("id", "year"), "first", weights = d$wv, ivar = TRUE)
  expect_true(all(is.finite(summary(fr, type = "event")$estimates$se)))
  skip_if_not_installed("etwfe")
  data(mpdta, package = "did")
  set.seed(5)
  mpdta$w <- runif(nrow(mpdta), 0.2, 3)
  ref <- suppressMessages(etwfe::etwfe(lemp ~ lpop, tvar = year, gvar = first.treat, data = mpdta,
                                       vcov = ~countyreal, weights = ~w))
  fit <- etwfe(lemp ~ lpop, mpdta, c("countyreal", "year"), "first.treat", weights = "w")
  for (type in c("simple", "event")) {
    e <- etwfe::emfx(ref, type = type, wts = "w")
    s <- summary(fit, type = type)
    tab <- if (type == "simple") s$overall else s$estimates
    expect_equal(tab$estimate, e$estimate, tolerance = 1e-8)
    expect_equal(tab$se, e$std.error, tolerance = 1e-6)
  }
})

test_that("etwfe xvar: effects by a covariate", {
  d <- et_panel(4)
  d$g <- as.numeric(d$id %% 2 == 0)
  post <- !is.na(d$first) & d$year >= d$first
  d$y2 <- d$y + 1.5 * post * d$g
  fit <- etwfe(y2 ~ x, d, c("id", "year"), "first", xvar = "g")
  s <- summary(fit)
  expect_equal(names(s$estimates)[1], "g")
  expect_equal(s$estimates$g, c(0, 1))
  expect_equal(diff(s$estimates$estimate), 1.5, tolerance = 0.25)
  # the overall ATT is the treated-count weighted mean of the by-xvar ATTs
  n <- table(d$g[post])
  expect_equal(s$overall$estimate, sum(n * s$estimates$estimate) / sum(n), tolerance = 1e-10)
  ev <- summary(fit, type = "event")$estimates
  expect_equal(nrow(ev), 2 * length(unique(ev$level)))
  expect_null(summary(fit, byxvar = FALSE)$estimates)
  expect_equal(summary(fit, type = "group", byxvar = FALSE)$estimates$level, c(2014, 2016))
  expect_output(print(fit), "effects by g")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(fit))
  expect_error(plot(summary(fit, type = "event")), "byxvar")
  b <- etwfe(bin ~ 1, d, c("id", "year"), "first", family = "binomial", xvar = "g",
             control = "never")
  expect_equal(nrow(summary(b, type = "calendar")$estimates), 2 * 5)
  d$gl <- d$g == 1
  expect_equal(summary(etwfe(y2 ~ x, d, c("id", "year"), "first", xvar = "gl"))$estimates$estimate,
               s$estimates$estimate, tolerance = 1e-10)
  d$gc <- as.character(d$g)
  expect_error(etwfe(y2 ~ 1, d, c("id", "year"), "first", xvar = "gc"), "numeric or logical")
  expect_error(etwfe(y2 ~ 1, d, c("id", "year"), "first", xvar = "nope"), "not found")
  expect_error(etwfe(y2 ~ 1, d, c("id", "year"), "first", xvar = c("g", "x")), "one column")
  expect_error(summary(etwfe(y2 ~ 1, d, c("id", "year"), "first"), byxvar = TRUE), "xvar")
})

test_that("etwfe ivar: individual fixed effects", {
  d <- et_panel(5)
  d$z <- rnorm(240)[d$id]
  f0 <- etwfe(y ~ 1, d, c("id", "year"), "first")
  f1 <- etwfe(y ~ 1, d, c("id", "year"), "first", ivar = TRUE)
  expect_true(f1$ivar)
  # balanced panel without time-varying controls: same cell effects and
  # ATTs as with cohort effects
  for (type in c("simple", "event")) {
    a <- summary(f0, type = type)
    b <- summary(f1, type = type)
    expect_equal(b$overall$estimate, a$overall$estimate, tolerance = 1e-8)
    expect_equal(b$estimates$estimate, a$estimates$estimate, tolerance = 1e-8)
  }
  expect_equal(f1$cells$estimate, f0$cells$estimate, tolerance = 1e-8)
  expect_output(print(f1), "individual and period")
  # a time-invariant control is absorbed
  f2 <- etwfe(y ~ x + z, d, c("id", "year"), "first", ivar = TRUE)
  expect_false("z" %in% names(coef(f2)))
  expect_true(is.finite(summary(f2)$overall$se))
  # unbalanced and weighted fits run and differ from the cohort version
  du <- d[-seq(1, nrow(d), by = 7), ]
  du$w <- runif(nrow(du), 0.5, 2)
  fu <- etwfe(y ~ x, du, c("id", "year"), "first", ivar = TRUE, weights = "w")
  expect_true(is.finite(summary(fu)$overall$estimate))
  expect_error(etwfe(cnt ~ 1, d, c("id", "year"), "first", family = "poisson", ivar = TRUE),
               "gaussian")
  expect_error(etwfe(y ~ 1, d, c("id", "year"), "first", ivar = NA), "TRUE or FALSE")
})

test_that("etwfe xvar and ivar match etwfe::etwfe and etwfe::emfx", {
  skip_if_not_installed("etwfe")
  data(mpdta, package = "did")
  mpdta$big <- as.integer(mpdta$lpop > stats::median(mpdta$lpop))
  mpdta$emp <- round(exp(mpdta$lemp))
  idx <- c("countyreal", "year")
  ref <- suppressMessages(etwfe::etwfe(emp ~ 1, tvar = year, gvar = first.treat, data = mpdta,
                                       xvar = big, family = "poisson", glm.tol = 1e-12,
                                       vcov = ~countyreal))
  fit <- etwfe(emp ~ 1, mpdta, idx, "first.treat", family = "poisson", xvar = "big")
  e <- etwfe::emfx(ref, type = "event")
  s <- summary(fit, type = "event")$estimates
  expect_equal(s$level, e$event)
  expect_equal(s$big, e$big)
  expect_equal(s$estimate, e$estimate, tolerance = 1e-6)
  expect_equal(s$se, e$std.error, tolerance = 1e-4)
  ref <- suppressMessages(etwfe::etwfe(lemp ~ lpop, tvar = year, gvar = first.treat,
                                       data = mpdta, ivar = countyreal, vcov = ~countyreal))
  fit <- etwfe(lemp ~ lpop, mpdta, idx, "first.treat", ivar = TRUE)
  e <- etwfe::emfx(ref, type = "group")
  s <- summary(fit, type = "group")$estimates
  expect_equal(s$estimate, e$estimate, tolerance = 1e-8)
  expect_equal(s$se, e$std.error, tolerance = 1e-6)
})
