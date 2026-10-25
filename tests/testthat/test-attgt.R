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
