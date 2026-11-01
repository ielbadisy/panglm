test_that("healthspend has the documented structure", {
  data(healthspend, package = "panglm", envir = environment())
  expect_identical(dim(healthspend), c(4439L, 15L))
  expect_identical(names(healthspend),
                   c("iso3c", "year", "country", "region", "income",
                     "chepc", "chereal", "chegdp", "oop", "gghe",
                     "gdppc", "lifeexp", "u5mr", "pop", "pop65"))
  expect_false(anyDuplicated(healthspend[c("iso3c", "year")]) > 0)
  expect_identical(nlevels(healthspend$iso3c), 193L)
  expect_identical(range(healthspend$year), c(2000L, 2022L))
  expect_true(all(table(healthspend$iso3c) == 23L))
  expect_true(is.factor(healthspend$iso3c))
  expect_true(is.factor(healthspend$region))
  expect_identical(levels(healthspend$income),
                   c("low", "lowermiddle", "uppermiddle", "high"))
  expect_type(healthspend$year, "integer")
  expect_type(healthspend$pop, "integer")
  expect_type(healthspend$country, "character")
  expect_identical(nlevels(healthspend$region), 7L)
})

test_that("healthspend values are in range and match known values", {
  data(healthspend, package = "panglm", envir = environment())
  h <- healthspend
  for (v in c("oop", "gghe")) {
    expect_true(all(h[[v]] >= 0 & h[[v]] <= 100, na.rm = TRUE))
  }
  expect_true(all(h$oop + h$gghe <= 100.01, na.rm = TRUE))
  expect_true(all(h$chepc > 0 & h$chereal > 0 & h$gdppc > 0, na.rm = TRUE))
  expect_true(all(h$lifeexp > 10 & h$lifeexp < 90))
  expect_identical(sum(is.na(h$lifeexp)), 0L)
  expect_identical(sum(is.na(h$chereal)), 217L)
  ## time-invariant region and income group
  expect_true(all(tapply(h$region, h$iso3c, function(x) length(unique(x))) == 1))
  expect_true(all(tapply(h$income, h$iso3c, function(x) length(unique(x))) == 1))
  r <- h[h$iso3c == "MAR" & h$year == 2021L, ]
  expect_equal(r$chepc, 477.35)
  expect_equal(r$oop, 41.795)
  expect_equal(r$lifeexp, 73.385)
  expect_equal(h$chepc[h$iso3c == "USA" & h$year == 2019L], 10550.22)
  expect_equal(r$chereal, round(r$chegdp / 100 * r$gdppc, 2))
  expect_equal(mean(h$lifeexp), 69.9566, tolerance = 1e-5)
  expect_equal(mean(h$chereal, na.rm = TRUE), 1479.2011, tolerance = 1e-6)
})
