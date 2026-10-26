att_panel <- function(seed, N = 300, effect = function(e) 1 + 0.5 * e) {
  set.seed(seed)
  d <- expand.grid(id = seq_len(N), year = 2011:2018)
  coh <- sample(c(2013, 2015, 2017, NA), N, replace = TRUE)
  d$first <- coh[d$id]
  a <- rnorm(N)
  d$x <- rnorm(N)[d$id]
  e <- d$year - d$first
  d$y <- a[d$id] + 0.1 * (d$year - 2011) + 0.2 * d$x * (d$year - 2011) +
    ifelse(!is.na(e) & e >= 0, effect(e), 0) + rnorm(nrow(d), sd = 0.5)
  d[sample(nrow(d)), ]
}

test_that("attgt recovers dynamic effects and has a consistent structure", {
  d <- att_panel(1)
  fit <- attgt(y ~ x, data = d, index = c("id", "year"), treat = "first")
  expect_s3_class(fit, "panglm_attgt")
  est <- fit$estimates
  expect_setequal(unique(est$group), c(2013, 2015, 2017))
  expect_equal(nrow(fit$inf), 300)
  expect_equal(ncol(fit$inf), nrow(est))
  expect_equal(unname(sqrt(diag(vcov(fit)))), est$se, tolerance = 1e-10)
  dyn <- summary(fit, type = "dynamic")
  post <- dyn$estimates[dyn$estimates$level >= 0, ]
  expect_equal(post$estimate, 1 + 0.5 * post$level, tolerance = 0.3)
  pre <- dyn$estimates[dyn$estimates$level < 0, ]
  expect_true(all(abs(pre$estimate) < 4 * pre$se))
  for (type in c("simple", "group", "calendar")) {
    s <- summary(fit, type = type)
    expect_true(is.finite(s$overall$estimate) && s$overall$se > 0)
  }
  expect_output(print(fit), "Group-time")
  expect_output(print(dyn), "event")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(fit))
  expect_silent(plot(summary(fit, type = "group")))
})

test_that("attgt options run and the universal base has a zero reference cell", {
  d <- att_panel(2)
  for (m in c("dr", "ipw", "reg")) for (cg in c("never", "notyet")) {
    fit <- attgt(y ~ x, d, c("id", "year"), "first", control = cg, method = m)
    expect_true(all(is.finite(fit$estimates$estimate)))
  }
  u <- attgt(y ~ 1, d, c("id", "year"), "first", base = "universal")
  ref <- u$estimates[u$estimates$event == -1, ]
  expect_equal(ref$estimate, rep(0, nrow(ref)))
  expect_error(attgt(y ~ 1, d, c("id", "year"), "nope"))
})

test_that("attgt matches did::att_gt and did::aggte", {
  skip_if_not_installed("did")
  data(mpdta, package = "did")
  for (cg in c("never", "notyet")) for (m in c("dr", "reg")) {
    ref <- did::att_gt("lemp", tname = "year", idname = "countyreal", gname = "first.treat",
                       xformla = ~lpop, data = mpdta, est_method = m, bstrap = FALSE,
                       cband = FALSE,
                       control_group = if (cg == "never") "nevertreated" else "notyettreated")
    fit <- attgt(lemp ~ lpop, mpdta, c("countyreal", "year"), "first.treat",
                 control = cg, method = m)
    expect_equal(fit$estimates$estimate, ref$att, tolerance = 1e-8)
    expect_equal(fit$estimates$se, ref$se, tolerance = 1e-6)
    for (type in c("simple", "dynamic", "group", "calendar")) {
      ag <- did::aggte(ref, type = type, bstrap = FALSE, cband = FALSE)
      s <- summary(fit, type = type)
      expect_equal(s$overall$estimate, ag$overall.att, tolerance = 1e-8)
      expect_equal(s$overall$se, ag$overall.se, tolerance = 1e-6)
      if (type != "simple") {
        expect_equal(s$estimates$estimate, ag$att.egt, tolerance = 1e-8)
        expect_equal(s$estimates$se, ag$se.egt, tolerance = 1e-6)
      }
    }
  }
})

test_that("attgt sampling weights: constant weights change nothing, varying ones are rejected", {
  d <- att_panel(3)
  fit <- attgt(y ~ x, d, c("id", "year"), "first")
  d$w2 <- 2.5
  fw <- attgt(y ~ x, d, c("id", "year"), "first", weights = "w2")
  expect_equal(fw$estimates$estimate, fit$estimates$estimate, tolerance = 1e-10)
  expect_equal(fw$estimates$se, fit$estimates$se, tolerance = 1e-10)
  expect_equal(summary(fw, type = "dynamic")$estimates, summary(fit, type = "dynamic")$estimates,
               tolerance = 1e-10)
  d$wv <- runif(nrow(d))
  expect_error(attgt(y ~ x, d, c("id", "year"), "first", weights = "wv"), "constant within")
  expect_error(attgt(y ~ x, d, c("id", "year"), "first", weights = -d$w2), "non-negative")
  expect_error(attgt(y ~ x, d, c("id", "year"), "first", weights = "nope"), "not found")
  # unit weights of 2 for half the units move the estimates
  d$wu <- ifelse(d$id %% 2 == 0, 2, 1)
  fu <- attgt(y ~ x, d, c("id", "year"), "first", weights = d$wu)
  expect_false(isTRUE(all.equal(fu$estimates$estimate, fit$estimates$estimate)))
})

test_that("weighted attgt matches did::att_gt(weightsname = )", {
  skip_if_not_installed("did")
  data(mpdta, package = "did")
  set.seed(4)
  wc <- stats::setNames(runif(500, 0.2, 3), unique(mpdta$countyreal))
  mpdta$w <- unname(wc[as.character(mpdta$countyreal)])
  for (cg in c("never", "notyet")) for (m in c("dr", "ipw")) {
    ref <- did::att_gt("lemp", tname = "year", idname = "countyreal", gname = "first.treat",
                       xformla = ~lpop, data = mpdta, est_method = m, bstrap = FALSE,
                       cband = FALSE, weightsname = "w",
                       control_group = if (cg == "never") "nevertreated" else "notyettreated")
    fit <- attgt(lemp ~ lpop, mpdta, c("countyreal", "year"), "first.treat",
                 control = cg, method = m, weights = "w")
    expect_equal(fit$estimates$estimate, ref$att, tolerance = 1e-8)
    expect_equal(fit$estimates$se, ref$se, tolerance = 1e-6)
    for (type in c("simple", "dynamic", "group")) {
      ag <- did::aggte(ref, type = type, bstrap = FALSE, cband = FALSE)
      s <- summary(fit, type = type)
      expect_equal(s$overall$estimate, ag$overall.att, tolerance = 1e-8)
      expect_equal(s$overall$se, ag$overall.se, tolerance = 1e-6)
    }
  }
})

test_that("attgt bootstrap: reproducible with a seed, global RNG untouched, uniform bands", {
  d <- att_panel(5)
  d$region <- (d$id - 1) %/% 10
  base <- attgt(y ~ x, d, c("id", "year"), "first")
  set.seed(99)
  before <- .Random.seed
  f1 <- attgt(y ~ x, d, c("id", "year"), "first", cband = TRUE, biters = 299, seed = 4)
  expect_identical(.Random.seed, before)
  f2 <- attgt(y ~ x, d, c("id", "year"), "first", cband = TRUE, biters = 299, seed = 4)
  expect_identical(f1$estimates, f2$estimates)
  expect_identical(summary(f1, type = "dynamic")$estimates, summary(f2, type = "dynamic")$estimates)
  expect_true(f1$boot$cband)
  expect_equal(f1$estimates$estimate, base$estimates$estimate)
  expect_true(all(f1$estimates$se != base$estimates$se))
  expect_gte(f1$boot$crit, stats::qnorm(0.975))
  est <- f1$estimates
  expect_true(all(c("band_lower", "band_upper") %in% names(est)))
  expect_equal(est$band_upper - est$estimate, f1$boot$crit * est$se)
  s <- summary(f1, type = "dynamic")
  expect_gte(s$crit, stats::qnorm(0.975))
  expect_true(all(s$estimates$band_lower <= s$estimates$lower))
  expect_silent(plot(s))
  expect_silent(plot(f1))
  # no band columns without cband, and the default results are unchanged
  f3 <- attgt(y ~ x, d, c("id", "year"), "first", bstrap = TRUE, biters = 99, seed = 1,
              multiplier = "mammen")
  expect_false("band_lower" %in% names(f3$estimates))
  expect_null(base$boot)
  expect_false("band_lower" %in% names(base$estimates))
  # clustered bootstrap runs; non-constant cluster is rejected
  f4 <- attgt(y ~ x, d, c("id", "year"), "first", cluster = "region", cband = TRUE,
              biters = 99, seed = 2)
  expect_gte(f4$boot$crit, stats::qnorm(0.975))
  d$bad <- seq_len(nrow(d))
  expect_error(attgt(y ~ x, d, c("id", "year"), "first", cluster = "bad"), "constant within")
  expect_error(attgt(y ~ x, d, c("id", "year"), "first", bstrap = TRUE, biters = 1), "biters")
})

test_that("clustered attgt matches did::att_gt(clustervars = )", {
  skip_if_not_installed("did")
  data(mpdta, package = "did")
  mpdta$state <- mpdta$countyreal %/% 1000
  ref <- suppressMessages(did::att_gt("lemp", tname = "year", idname = "countyreal",
                                      gname = "first.treat", xformla = ~lpop, data = mpdta,
                                      clustervars = c("countyreal", "state"),
                                      bstrap = FALSE, cband = FALSE))
  fit <- attgt(lemp ~ lpop, mpdta, c("countyreal", "year"), "first.treat", cluster = "state")
  expect_equal(fit$estimates$se, ref$se, tolerance = 1e-6)
  ag <- did::aggte(ref, type = "dynamic", bstrap = FALSE, cband = FALSE)
  s <- summary(fit, type = "dynamic")
  expect_equal(s$estimates$se, ag$se.egt, tolerance = 1e-6)
  expect_equal(s$overall$se, ag$overall.se, tolerance = 1e-6)
})

test_that("attgt(), did() and eventstudy() ignore unused levels of a factor index", {
  set.seed(1)
  d <- expand.grid(province = 1:40, year = 2012:2019)
  first <- rep(c(2015, 2017, NA, NA), each = 10)
  d$adoption <- first[d$province]
  e <- d$year - d$adoption
  d$y <- rnorm(40)[d$province] + 0.3 * (d$year - 2012) +
    ifelse(!is.na(e) & e >= 0, 1, 0) + rnorm(nrow(d))
  d$province <- factor(d$province)
  s <- subset(d, as.integer(province) > 5)
  s2 <- droplevels(s)
  ix <- c("province", "year")
  expect_equal(coef(attgt(y ~ 1, data = s, index = ix, treat = "adoption")),
               coef(attgt(y ~ 1, data = s2, index = ix, treat = "adoption")))
  expect_equal(coef(did(y ~ 1, data = s, index = ix, treat = "adoption")),
               coef(did(y ~ 1, data = s2, index = ix, treat = "adoption")))
  expect_equal(coef(eventstudy(y ~ 1, data = s, index = ix, treat = "adoption")),
               coef(eventstudy(y ~ 1, data = s2, index = ix, treat = "adoption")))
})
