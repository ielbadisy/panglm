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

test_that("eventplot compares the three estimators", {
  d <- did_plot_panel(4)
  es <- eventstudy(y ~ 1, d, index = c("id", "t"), treat = "first", leads = 3, lags = 3)
  cs <- attgt(y ~ 1, d, index = c("id", "t"), treat = "first")
  ew <- etwfe(y ~ 1, d, index = c("id", "t"), treat = "first", control = "never")
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off())
  out <- expect_invisible(eventplot(TWFE = es, CS = cs, ETWFE = ew, min_e = -3, max_e = 3))
  expect_setequal(out$model, c("TWFE", "CS", "ETWFE"))
  expect_true(all(out$event >= -3 & out$event <= 3))
  out <- eventplot(es, summary(cs, type = "dynamic"), summary(ew, type = "event"),
                   legend = FALSE)
  expect_setequal(out$model, c("TWFE event study", "Callaway-Sant'Anna", "ETWFE"))
  expect_equal(unique(eventplot(cs, cs)$model), c("Callaway-Sant'Anna",
                                                  "Callaway-Sant'Anna 1"))
  expect_error(eventplot(), "at least one")
  expect_error(eventplot(lm(y ~ t, d)), "eventplot")
  expect_error(eventplot(es, cs, labels = "a"), "one label")
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

test_that("top legends wrap into columns in a narrow panel", {
  lab <- c("TWFE event study", "Callaway-Sant'Anna", "Extended TWFE")
  pdf(NULL, width = 9, height = 4.5); on.exit(dev.off())
  wide <- panglm:::did_leglayout(lab, 0.8)
  expect_equal(c(wide$ncol, wide$nrow), c(3, 1))
  op <- par(mfrow = c(1, 3))
  narrow <- panglm:::did_leglayout(lab, 0.8)
  expect_lt(narrow$ncol, 3)
  expect_gt(narrow$nrow, 1)
  expect_gt(panglm:::did_legroom(c(0, 1), narrow), panglm:::did_legroom(c(0, 1), wide))
  d <- did_plot_panel(6)
  es <- eventstudy(y ~ 1, d, index = c("id", "t"), treat = "first", leads = 3, lags = 3)
  cs <- attgt(y ~ 1, d, index = c("id", "t"), treat = "first")
  expect_no_error(eventplot("TWFE event study" = es, "Callaway-Sant'Anna" = cs))
  expect_no_error(plot(cs, legend = TRUE))
  par(op)
})

test_that("event-study styles run for every plot method", {
  d <- did_plot_panel(7)
  cs <- attgt(y ~ 1, d, index = c("id", "t"), treat = "first", cband = TRUE,
              biters = 99, seed = 1)
  ew <- etwfe(y ~ 1, d, index = c("id", "t"), treat = "first", control = "never")
  es <- eventstudy(y ~ 1, d, index = c("id", "t"), treat = "first", leads = 3, lags = 3)
  pdf(NULL); on.exit(dev.off())
  for (style in c("classic", "points", "line")) {
    for (type in c("dynamic", "group", "calendar")) {
      expect_no_error(plot(cs, type = type, style = style, legend = TRUE))
    }
    expect_no_error(plot(summary(cs, type = "dynamic"), style = style, min_e = -2))
    for (type in c("event", "group", "calendar")) {
      expect_no_error(plot(ew, type = type, style = style, legend = TRUE))
    }
    expect_no_error(plot(summary(ew, type = "event"), style = style))
    expect_identical(plot(es, style = style, legend = TRUE), es)
    out <- eventplot(TWFE = es, CS = cs, ETWFE = ew, style = style)
    expect_setequal(out$model, c("TWFE", "CS", "ETWFE"))
  }
  expect_no_error(plot(es, style = "classic", onset = -0.5, onset_label = NULL,
                       caps = FALSE, main = "m"))
  expect_no_error(plot(cs, style = "points", onset_label = "Policy", caps = TRUE))
  expect_error(plot(es, style = "bars"), "should be one of")
  expect_error(plot(es, onset = NA), "onset")
  expect_error(plot(es, onset_label = 1), "onset_label")
  expect_error(plot(es, caps = NA), "caps")
  expect_error(plot(es, ref = "a"), "ref")
  st <- panglm:::did_style("classic")
  expect_equal(st[c("onset", "onset_label", "caps")],
               list(onset = 0, onset_label = "Intervention", caps = TRUE))
  st <- panglm:::did_style("line")
  expect_equal(st$onset, -0.5)
  expect_null(st$onset_label)
  expect_false(st$caps)
})

test_that("the reference period is marked only where the estimate is exactly zero", {
  d <- did_plot_panel(8)
  pdf(NULL); on.exit(dev.off())
  es <- eventstudy(y ~ 1, d, index = c("id", "t"), treat = "first", leads = 3, lags = 3)
  tab <- panglm:::did_refmark(panglm:::es_plottab(es), es$ref)
  expect_equal(tab$level[tab$ref], es$ref)
  expect_equal(tab$estimate[tab$ref], 0)
  # varying base period: no omitted period, nothing marked
  cv <- attgt(y ~ 1, d, index = c("id", "t"), treat = "first")
  sv <- summary(cv, type = "dynamic")
  tv <- panglm:::did_refmark(panglm:::did_plottab(sv$estimates), -1)
  expect_false(any(tv$ref))
  # universal base period: event time -1 is normalized to zero
  cu <- attgt(y ~ 1, d, index = c("id", "t"), treat = "first", base = "universal")
  su <- summary(cu, type = "dynamic")
  tu <- panglm:::did_refmark(panglm:::did_plottab(su$estimates), -1)
  expect_equal(tu$level[tu$ref], -1)
  expect_identical(plot(su), su)
  expect_false(any(panglm:::did_refmark(tu, NULL)$ref))
  ew <- etwfe(y ~ 1, d, index = c("id", "t"), treat = "first", control = "never")
  te <- panglm:::etwfe_plottab(summary(ew, type = "event"))
  expect_equal(te$level[te$ref], -1)
})

test_that("attgt heat maps show the ATT(g, t) table and flag pre-treatment cells", {
  d <- did_plot_panel(3)
  fit <- attgt(y ~ 1, d, index = c("id", "t"), treat = "first")
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off())
  hm <- expect_invisible(plot(fit, type = "heatmap"))
  est <- fit$estimates
  expect_equal(hm$estimate, est$estimate[is.finite(est$estimate)])
  expect_identical(hm$pre, hm$time < hm$group)
  expect_true(any(hm$pre) && any(!hm$pre))
  hw <- plot(fit, type = "heatmap", cells = TRUE, digits = 1, min_e = 3, max_e = 6,
             ylim = c(-1, 1), main = "m")
  expect_true(all(hw$time >= 3 & hw$time <= 6))
  expect_error(plot(fit, type = "heatmap", cells = NA), "cells")
  expect_error(plot(fit, type = "heatmap", min_e = 5, max_e = 2), "min_e")
})

test_that("etwfe heat maps aggregate back to the event-study estimates", {
  d <- did_plot_panel(4)
  fit <- etwfe(y ~ 1, d, index = c("id", "t"), treat = "first", control = "never")
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off())
  hm <- expect_invisible(plot(fit, type = "heatmap", cells = TRUE))
  # one row per cohort and period, the reference cell g - 1 at zero
  expect_equal(nrow(hm), 2 * 8)
  expect_equal(hm$estimate[hm$time == hm$group - 1], c(0, 0))
  expect_identical(hm$pre, hm$time < hm$group)
  # cell-size weighted means of ATT(g, t) reproduce the event aggregation
  rows <- fit$rows
  n <- table(paste(rows$group, rows$time)[is.finite(rows$group)])
  hm$n <- as.numeric(n[paste(hm$group, hm$time)])
  ev <- summary(fit, type = "event")$estimates
  agg <- sapply(split(hm, hm$event), function(z) sum(z$n * z$estimate) / sum(z$n))
  expect_equal(unname(agg[as.character(ev$level)]), ev$estimate, tolerance = 1e-10)
  # with the not-yet-treated control there are no pre-treatment cells
  fn <- etwfe(y ~ 1, d, index = c("id", "t"), treat = "first", control = "notyet")
  hn <- plot(fn, type = "heatmap")
  expect_false(any(hn$pre))
  expect_equal(nrow(hn), (8 - 4 + 1) + (8 - 6 + 1))
})
