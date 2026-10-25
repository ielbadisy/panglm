# plm::pgmm() evaluates a call to plm() in its parent frame
predict_pgmm <- function(formula, data, ...) {
  plm <- plm::plm
  plm::pgmm(formula, data = data, ...)
}

key_match <- function(a, b) match(paste(a$id, a$t), paste(b$id, b$t))

test_that("predict.paniv agrees with plm and AER::ivreg, in sample and on newdata", {
  skip_if_missing("plm")
  d <- methods_data("iv")
  set.seed(9)
  d <- d[sample(nrow(d)), ]
  pd <- plm::pdata.frame(d, index = c("id", "t"))
  f <- y ~ x + w | z + z2 + w
  for (m in c("pooling", "within", "random")) {
    a <- paniv(f, data = d, index = c("id", "t"), model = m)
    p <- plm::plm(f, data = pd, model = m)
    ref <- as.numeric(predict(p))[key_match(d, pd)]
    expect_equal(predict(a), ref, tolerance = 1e-10)
    expect_equal(predict(a, newdata = d), predict(a), tolerance = 1e-10)
  }
  a <- paniv(f, data = d, index = c("id", "t"), model = "within", effect = "twoways")
  p <- plm::plm(f, data = pd, model = "within", effect = "twoways")
  expect_equal(predict(a), as.numeric(predict(p))[key_match(d, pd)], tolerance = 1e-10)
  expect_equal(predict(a, newdata = d), predict(a), tolerance = 1e-10)
  skip_if_missing("AER")
  a <- paniv(f, data = d, index = c("id", "t"), model = "pooling")
  iv <- AER::ivreg(f, data = d)
  expect_equal(predict(a, newdata = d), unname(predict(iv, newdata = d)), tolerance = 1e-10)
})

test_that("predict.paniv handles missing rows and unseen individuals", {
  d <- methods_data("iv")
  d$x[5] <- NA
  a <- paniv(y ~ x + w | z + z2 + w, data = d, index = c("id", "t"), model = "within")
  expect_length(predict(a), nrow(d) - 1L)
  nd <- d[1:4, ]
  nd$id <- 999
  expect_error(predict(a, newdata = nd), "not present")
  expect_equal(predict(a, newdata = nd, allow.new.levels = TRUE),
               as.numeric(cbind(nd$x, nd$w) %*% coef(a)))
})

test_that("predict.pangmm reproduces the transformed equations of plm::pgmm", {
  skip_if_missing("plm")
  g <- methods_data("gmm")
  fg <- methods_fit("pangmm")
  pg <- predict_pgmm(y ~ lag(y, 1) + x | lag(y, 2:99), data = plm::pdata.frame(g, index = c("id", "t")),
                  effect = "twoways", model = "twosteps", transformation = "d")
  expect_equal(coef(fg), coef(pg)[names(coef(fg))], tolerance = 1e-8)
  dif <- predict(fg, newdata = g, type = "difference")
  ref <- unlist(lapply(pg$model, function(m) m[, -1, drop = FALSE] %*% coef(pg)))
  mine <- unlist(lapply(split(dif, g$id), function(v) v[3:6]))
  expect_equal(unname(mine), unname(ref), tolerance = 1e-8)
  y <- unlist(lapply(pg$model, function(m) m[, 1]))
  expect_equal(unname(unlist(residuals(fg))), unname(y - mine), tolerance = 1e-8)

  sy <- pangmm(y ~ lag(y, 1) + x | lag(y, 2:99), data = g, index = c("id", "t"),
               model = "twostep", transformation = "ld")
  ps <- predict_pgmm(y ~ lag(y, 1) + x | lag(y, 2:99), data = plm::pdata.frame(g, index = c("id", "t")),
                  effect = "twoways", model = "twosteps", transformation = "ld")
  lv <- predict(sy, newdata = g)
  dif <- predict(sy, newdata = g, type = "difference")
  mine <- unlist(lapply(split(seq_len(nrow(g)), g$id), function(r) c(dif[r][3:6], lv[r][2:6])))
  ref <- unlist(lapply(ps$model, function(m) m[, -1, drop = FALSE] %*% coef(ps)))
  expect_equal(unname(mine), unname(ref), tolerance = 1e-8)
  expect_true(all(is.na(lv[g$t == 1])))
})

test_that("predict.pangmm checks its inputs", {
  fg <- methods_fit("pangmm")
  expect_error(predict(fg), "newdata")
  g <- methods_data("gmm")
  expect_error(predict(fg, newdata = g[, c("y", "x")]), "index")
  g7 <- rbind(g, transform(g[g$t == 6, ], t = 7))
  expect_warning(p <- predict(fg, newdata = g7), "no time effect")
  expect_true(all(is.na(p[g7$t == 7])))
})

test_that("predict.panselect agrees with sampleSelection::heckit on one period", {
  skip_if_missing("sampleSelection")
  d <- methods_data("sel")
  d1 <- d[d$t == 1, ]
  ps <- panselect(s ~ x + z, y ~ x, data = d1, index = c("id", "t"), model = "pooling",
                  mundlak = FALSE, lambda = "common", boot = 0)
  h <- sampleSelection::selection(s ~ x + z, y ~ x, data = d1, method = "2step")
  expect_equal(predict(ps, newdata = d1),
               unname(predict(h, newdata = d1, part = "outcome", type = "unconditional")), tolerance = 1e-6)
  expect_equal(predict(ps, newdata = d1, type = "conditional"),
               unname(predict(h, newdata = d1, part = "outcome", type = "conditional")[, "E[yo|ys=1]"]),
               tolerance = 1e-6)
  expect_equal(predict(ps, newdata = d1, type = "selection"),
               unname(predict(h, newdata = d1, part = "selection", type = "response")), tolerance = 1e-6)
})

test_that("predict.panselect matches the fit on the selected rows", {
  d <- methods_data("sel")
  pp <- panselect(s ~ x + z, y ~ x, data = d, index = c("id", "t"), model = "pooling", boot = 0)
  ok <- d$s == 1
  cond <- predict(pp, newdata = d, type = "conditional")
  expect_equal(cond[ok], d$y[ok] - pp$residuals, tolerance = 1e-10)
  pr <- predict(pp, newdata = d, type = "selection")
  expect_true(all(pr > 0 & pr < 1))
  expect_error(predict(pp), "newdata")
  w <- methods_fit("panselect")
  expect_false("(Intercept)" %in% names(coef(w)))
  expect_length(predict(w, newdata = d), nrow(d))
})
