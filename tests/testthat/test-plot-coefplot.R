test_that("coefplot compares several models and returns them", {
  data(copd)
  f <- exacerbations ~ crp + treatment + age
  pool <- panglm(f, data = copd, index = c("id", "visit"), model = "pooling",
                 family = "poisson")
  fe <- panglm(f, data = copd, index = c("id", "visit"), model = "random",
               family = "poisson")
  g <- glm(f, data = copd, family = poisson)
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  old <- graphics::par("mai")
  out <- coefplot(Pooled = pool, FE = fe, glm = g)
  expect_named(out, c("Pooled", "FE", "glm"))
  expect_equal(graphics::par("mai"), old)
  expect_invisible(coefplot(pool, fe, parm = "crp", level = 0.9))
  expect_named(coefplot(pool, fe), c("pool", "fe"))
  expect_named(coefplot(list(a = pool, b = fe)), c("a", "b"))
  expect_invisible(coefplot(pool, fe, labels = c("one", "two"), ref = NULL))
  expect_error(coefplot(pool, labels = c("a", "b")), "labels")
  expect_error(coefplot(pool, parm = "zz"), "parm")
  expect_error(coefplot(), "at least one")

  tab <- panglm:::.panglm_coef_matrices(
    list(a = list(est = coef(pool), se = sqrt(diag(vcov(pool)))),
         b = list(est = coef(fe), se = sqrt(diag(vcov(fe))))),
    parm = NULL, level = 0.95)
  expect_equal(rownames(tab$est), c("crp", "treatment", "age"))
  expect_equal(unname(tab$lower[, "a"]), unname(confint.default(g)[-1, 1]), tolerance = 1e-3)
})
