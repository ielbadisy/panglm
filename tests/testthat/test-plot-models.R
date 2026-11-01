with_pdf <- function(code) {
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  force(code)
}

test_that("plot.twopart draws either or both parts", {
  data(copd)
  set.seed(1)
  copd$cost <- ifelse(copd$exacerbations > 0,
                      exp(6 + 0.1 * copd$crp + rnorm(nrow(copd), 0, 0.5)), 0)
  tp <- twopart(cost ~ treatment + age + crp, data = copd, index = c("id", "visit"))
  with_pdf({
    expect_identical(plot(tp), tp)
    expect_invisible(plot(tp, which = "binary"))
    expect_invisible(plot(tp, which = "positive", parm = "crp"))
  })
})

test_that("plot.panfmm draws class coefficients and posterior quality", {
  data(copd)
  fit <- panfmm(exacerbations ~ crp + treatment, data = copd, index = c("id", "visit"),
                k = 2, family = "poisson", nstart = 2, seed = 1)
  with_pdf({
    expect_identical(plot(fit), fit)
    expect_invisible(plot(fit, which = "posterior"))
    expect_invisible(plot(fit, parm = "crp"))
  })
})

test_that("GMM, IV, selection and hurdle fits have coefficient plots", {
  data(copd)
  gm <- pangmm(fev1 ~ lag(fev1, 1) + crp | lag(fev1, 2:99), data = copd,
               index = c("id", "visit"))
  iv <- paniv(fev1 ~ crp + treatment | age + treatment, data = copd,
              index = c("id", "visit"), model = "pooling")
  set.seed(1)
  d <- data.frame(id = rep(1:150, each = 4), wave = rep(1:4, 150))
  d$treat <- rnorm(600)
  d$reminder <- rnorm(600)
  v <- rnorm(600)
  d$respond <- as.numeric(0.3 + 0.5 * d$treat + d$reminder + v > 0)
  d$qol <- ifelse(d$respond == 1, 1 + 0.5 * d$treat + 0.7 * v + rnorm(600), NA)
  sel <- panselect(respond ~ treat + reminder, qol ~ treat, data = d,
                   index = c("id", "wave"), model = "pooling", boot = 0)
  hu <- suppressMessages(panhurdle(exacerbations ~ crp, data = copd, index = c("id", "visit")))
  with_pdf({
    expect_identical(plot(gm), gm)
    expect_invisible(plot(gm, parm = "crp"))
    expect_identical(plot(iv), iv)
    expect_identical(plot(sel), sel)
    expect_invisible(plot(sel, which = "selection"))
    expect_invisible(plot(sel, which = "outcome", parm = "treat"))
    expect_invisible(plot(hu, parm = "crp"))
  })
})
