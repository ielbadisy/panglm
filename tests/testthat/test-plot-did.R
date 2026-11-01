did_plot_panel <- function(seed = 1, N = 150) {
  set.seed(seed)
  d <- expand.grid(id = seq_len(N), t = 1:8)
  d$first <- rep(c(4, 6, NA), length.out = N)[d$id]
  d$y <- rnorm(N)[d$id] + 0.2 * d$t +
    ifelse(!is.na(d$first) & d$t >= d$first, 1, 0) + rnorm(nrow(d))
  d
}

test_that("attgt plots run for every type and return invisibly", {
  d <- did_plot_panel()
  fit <- attgt(y ~ 1, d, index = c("id", "t"), treat = "first", cband = TRUE,
               biters = 99, seed = 1)
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off())
  for (type in c("dynamic", "group", "calendar", "gt")) {
    expect_no_error(expect_invisible(plot(fit, type = type)))
    expect_no_error(plot(fit, type = type, legend = TRUE, main = "m", ylim = c(-2, 3)))
  }
  s <- expect_invisible(plot(fit, min_e = -2, max_e = 2))
  expect_s3_class(s, "summary.panglm_attgt")
  expect_equal(s$estimates$level, -2:2)
  gt <- plot(fit, type = "gt", ncol = 1)
  expect_setequal(unique(gt$group), c(4, 6))
  expect_no_error(plot(summary(fit, type = "dynamic"), min_e = -1, legend = TRUE))
  expect_error(plot(fit, min_e = 3, max_e = 1), "min_e")
  expect_error(plot(fit, legend = NA), "legend")
  expect_error(plot(fit, type = "gt", ncol = 0), "ncol")
})

test_that("uniform bands are drawn only when the fit has them", {
  d <- did_plot_panel(2)
  fb <- attgt(y ~ 1, d, index = c("id", "t"), treat = "first", cband = TRUE,
              biters = 99, seed = 1)
  f0 <- attgt(y ~ 1, d, index = c("id", "t"), treat = "first")
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off())
  s <- plot(fb)
  expect_true(all(c("band_lower", "band_upper") %in% names(s$estimates)))
  tab <- panglm:::did_drawes(panglm:::did_plottab(s$estimates))
  expect_true(attr(tab, "band"))
  expect_true(all(tab$band_lower <= tab$lower))
  s0 <- plot(f0)
  expect_false("band_lower" %in% names(s0$estimates))
  expect_false(attr(panglm:::did_drawes(panglm:::did_plottab(s0$estimates)), "band"))
})

test_that("etwfe and eventstudy plots share the drawer and take an event window", {
  d <- did_plot_panel(3)
  ew <- etwfe(y ~ 1, d, index = c("id", "t"), treat = "first", control = "never")
  es <- eventstudy(y ~ 1, d, index = c("id", "t"), treat = "first", leads = 3, lags = 3)
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off())
  s <- expect_invisible(plot(ew, min_e = -2, max_e = 2, legend = TRUE))
  expect_s3_class(s, "summary.panglm_etwfe")
  ptab <- panglm:::etwfe_plottab(s, -2, 2)
  expect_equal(ptab$level, -2:2)
  expect_equal(ptab$level[ptab$ref], -1)
  for (type in c("group", "calendar")) expect_no_error(plot(ew, type = type, legend = TRUE))
  expect_no_error(expect_invisible(plot(es)))
  expect_identical(plot(es, min_e = -2, max_e = 1, legend = TRUE, main = "es"), es)
  ptab <- panglm:::es_plottab(es, -2, 1)
  expect_equal(ptab$level, -2:1)
  expect_equal(ptab$level[ptab$ref], es$ref)
  expect_error(plot(es, min_e = 5), "window")
})

test_that("treatmentplot shows cohorts and scales to many units", {
  d <- did_plot_panel(5, N = 30)
  d$first[is.na(d$first)] <- 0
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off())
  out <- expect_invisible(treatmentplot(d, index = c("id", "t"), treat = "first"))
  expect_equal(out$cohort, c(4, 6, NA))
  expect_equal(out$n, c(10, 10, 10))
  expect_no_error(treatmentplot(d[-(1:5), ], c("id", "t"), "first", legend = FALSE,
                                labels = FALSE, main = "unbalanced", ylab = "Unit"))
  set.seed(1)
  big <- expand.grid(id = 1:2000, t = 1:10)
  big$first <- sample(c(3, 5, 8, NA), 2000, TRUE)[big$id]
  expect_no_error(treatmentplot(big, c("id", "t"), "first"))
  expect_error(treatmentplot(d, c("id", "t"), "first", col = "red"), "three")
  d$first <- as.character(d$first)
  expect_error(treatmentplot(d, c("id", "t"), "first"), "numeric")
})
