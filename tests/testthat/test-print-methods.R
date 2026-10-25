# The print and summary-print methods share one layout: a header naming the
# estimator, an "N = ..." line (with "vcov: ..." in summaries), a
# printCoefmat() table with significance stars, and an invisible return.

expect_print_layout <- function(obj, vcov_line = FALSE) {
  out <- utils::capture.output(res <- withVisible(print(obj)))
  expect_false(res$visible)
  expect_identical(res$value, obj)
  expect_true(any(grepl("\\bN = ", out)))
  if (vcov_line) expect_true(any(grepl("N = .*vcov:", out)))
  out
}

test_that("print methods of the model classes share the layout", {
  for (nm in c("twopart", "pangmm", "paniv", "panqr", "panselect", "panfmm")) {
    expect_print_layout(methods_fit(nm))
  }
  es <- expect_print_layout(methods_fit("eventstudy"), vcov_line = TRUE)
  copd <- methods_data("copd")
  h <- suppressMessages(panglm_hurdle(exacerbations ~ crp, data = copd, index = c("id", "visit")))
  out <- expect_print_layout(h)
  expect_true(any(grepl("Zero part", out)))
})

test_that("summary print methods show the vcov line and significance stars", {
  for (nm in c("twopart", "pangmm", "paniv", "panqr", "panselect", "panfmm", "eventstudy")) {
    s <- summary(methods_fit(nm))
    out <- expect_print_layout(s)
    expect_true(any(grepl("vcov", out)), info = nm)
    expect_true(any(grepl("Signif. codes", out)), info = nm)
  }
  out <- utils::capture.output(print(summary(methods_fit("twopart"))))
  expect_true(any(grepl("^Call:", out)))
})

test_that("summary of the hurdle model agrees with coef, vcov and confint", {
  copd <- methods_data("copd")
  h <- suppressMessages(panglm_hurdle(exacerbations ~ crp, data = copd, index = c("id", "visit"),
                                      count_family = "negbin"))
  s <- summary(h)
  expect_s3_class(s, "summary.panglm_hurdle")
  tab <- rbind(s$zero, s$count)
  expect_equal(unname(tab[, 1]), unname(coef(h)))
  expect_equal(unname(tab[, 2]), unname(sqrt(diag(vcov(h)))))
  expect_equal(unname(tab[, 4]), 2 * pnorm(-abs(unname(tab[, 3]))))
  ci <- confint(h)
  expect_equal(unname(ci[, 1]), unname(tab[, 1] - qnorm(0.975) * tab[, 2]))
  out <- expect_print_layout(s, vcov_line = TRUE)
  expect_true(any(grepl("Theta", out)))
})
