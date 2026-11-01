test_that("single-panel plots fill the next panel of a par(mfrow) grid", {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  data(copd)
  ix <- c("id", "visit")
  fit <- panglm(fev1 ~ crp + age + treatment, data = copd, index = ix,
    model = "random")
  draws <- list(
    function() plot(fit, which = "density"),
    function() plot(fit, which = "correlation"),
    function() treatmentplot(copd, index = ix, treat = "adopt"),
    function() panelplot(copd, index = ix, value = "fev1")
  )
  for (draw in draws) {
    graphics::par(mfrow = c(1, 2))
    draw()
    expect_equal(graphics::par("mfg")[1:2], c(1L, 1L))
    expect_equal(graphics::par("mfrow"), c(1L, 2L))
  }
})
