test_that("cigar has the documented structure", {
  data(cigar, package = "panglm", envir = environment())
  expect_identical(dim(cigar), c(1380L, 13L))
  expect_identical(names(cigar),
                   c("state", "year", "price", "pop", "pop16", "cpi", "ndi",
                     "sales", "pimin", "lnc", "lnp", "lnpn", "lny"))
  expect_false(anyDuplicated(cigar[c("state", "year")]) > 0)
  expect_true(is.factor(cigar$state))
  expect_identical(nlevels(cigar$state), 46L)
  expect_false(any(c("AK", "CO", "HI", "NC", "OR") %in% levels(cigar$state)))
  expect_type(cigar$year, "integer")
  expect_identical(range(cigar$year), c(1963L, 1992L))
  expect_true(all(table(cigar$state) == 30L))
  expect_false(anyNA(cigar))
})

test_that("cigar values match the source and the derived logs", {
  data(cigar, package = "panglm", envir = environment())
  ut <- cigar[cigar$state == "UT" & cigar$year == 1988L, ]
  expect_equal(ut$sales, 55)                 # Baltagi (2005), p. 157
  expect_equal(round(cigar$sales[cigar$state == "NV" & cigar$year == 1988L]), 142)
  expect_equal(cigar$lnc, log(cigar$sales * cigar$pop / cigar$pop16))
  expect_equal(cigar$lnp, log(cigar$price / cigar$cpi * 100))
  expect_equal(cigar$lnpn, log(cigar$pimin / cigar$cpi * 100))
  expect_equal(cigar$lny, log(cigar$ndi / cigar$cpi * 100))
  expect_equal(mean(cigar$sales), 123.9509, tolerance = 1e-6)
  expect_equal(mean(cigar$price), 68.6999, tolerance = 1e-6)
})

test_that("cigar reproduces Baltagi (2005) Table 8.1, OLS and within", {
  data(cigar, package = "panglm", envir = environment())
  f <- lnc ~ lag(lnc, 1) + lnp + lnpn + lny
  ols <- panglm(f, data = cigar, index = c("state", "year"), model = "pooling")
  fe <- panglm(f, data = cigar, index = c("state", "year"), model = "within",
               effect = "twoways")
  expect_equal(round(unname(coef(ols)[-1]), c(2, 3, 3, 2)),
               c(0.97, -0.090, 0.024, -0.03))
  expect_equal(round(unname(coef(fe)), c(2, 3, 3, 2)),
               c(0.83, -0.299, 0.034, 0.10))
})
