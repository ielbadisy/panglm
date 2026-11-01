with_pdf <- function(code) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  force(code)
}

ix <- c("id", "visit")

test_that("the correlation heatmap draws cov2cor(vcov()) and restores par", {
  data(copd)
  fit <- panglm(fev1 ~ crp + age + treatment, data = copd, index = ix, model = "random")
  with_pdf({
    old <- graphics::par(c("mai", "pty"))
    expect_invisible(plot(fit, which = "correlation"))
    expect_identical(plot(fit, which = "correlation", cells = TRUE, parm = 1:3), fit)
    expect_equal(graphics::par(c("mai", "pty")), old)
    R <- panglm:::.panglm_corrplot(vcov(fit))
    expect_equal(R, stats::cov2cor(vcov(fit))[-1, -1])
    R <- panglm:::.panglm_corrplot(vcov(fit), parm = c("age", "crp"))
    expect_identical(rownames(R), c("age", "crp"))
    expect_error(plot(fit, which = "correlation", parm = "zz"), "parm")
    expect_error(plot(fit, which = "correlation", parm = "crp"), "two terms")
  })
})

test_that("degenerate variances are drawn as missing cells", {
  V <- diag(c(1, 0, 4))
  V[1, 3] <- V[3, 1] <- 1
  dimnames(V) <- list(c("a", "b", "c"), c("a", "b", "c"))
  with_pdf(R <- panglm:::.panglm_corrplot(V))
  expect_true(all(is.na(R["b", ])))
  expect_equal(R["a", "c"], 0.5)
})

test_that("multi-part models expose the correlation view", {
  data(copd)
  set.seed(1)
  copd$cost <- ifelse(copd$exacerbations > 0, stats::rgamma(nrow(copd), 2, 0.01), 0)
  tp <- twopart(cost ~ crp + treatment, data = copd, index = ix)
  hu <- panhurdle(exacerbations ~ crp, data = copd, index = ix)
  with_pdf({
    expect_identical(plot(tp, which = "correlation", cells = TRUE), tp)
    expect_identical(plot(hu, which = "correlation"), hu)
    R <- panglm:::.panglm_corrplot(vcov(hu))
    expect_equal(R["zero:crp", "count:crp"], 0)
  })
})
