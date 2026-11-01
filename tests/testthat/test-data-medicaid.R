test_that("medicaid has the documented structure", {
  data(medicaid, package = "panglm", envir = environment())
  expect_identical(dim(medicaid), c(561L, 20L))
  expect_identical(names(medicaid),
                   c("state", "year", "abb", "fips", "region", "pop",
                     "uninsured", "adultpop", "uninsadult", "medicaidadult",
                     "povpop", "povbase", "uninspov", "medicaidpov", "lowpop",
                     "uninslow", "adopt", "expdate", "expanded", "prior"))
  expect_false(anyDuplicated(medicaid[c("state", "year")]) > 0)
  expect_identical(nlevels(medicaid$state), 51L)
  expect_true(all(table(medicaid$state) == 11L))
  expect_identical(range(medicaid$year), c(2009L, 2019L))
  expect_true(all(is.na(medicaid$uninslow) == (medicaid$year <= 2010L)))
  expect_false(anyNA(medicaid[setdiff(names(medicaid),
                                      c("lowpop", "uninslow", "adopt",
                                        "expdate"))]))
  rd <- tryCatch(tools::Rd_db("panglm")[["medicaid.Rd"]], error = function(e) NULL)
  if (!is.null(rd)) {
    expect_match(paste(as.character(rd), collapse = ""), "and 20\\s+columns")
  }
})

test_that("medicaid matches published ACS rates and KFF expansion dates", {
  data(medicaid, package = "panglm", envir = environment())
  # Census P60-253, Appendix Table A-1 (percent uninsured, all ages)
  pub <- c(KY_2013 = 14.3, KY_2014 = 8.5, TX_2013 = 22.1, TX_2014 = 19.1,
           CA_2013 = 17.2, CA_2014 = 12.4, MA_2014 = 3.3)
  key <- paste(medicaid$abb, medicaid$year, sep = "_")
  got <- medicaid$uninsured[match(names(pub), key)]
  expect_true(all(abs(got - pub) <= 0.055))
  st <- unique(medicaid[c("abb", "adopt", "expdate")])
  expect_identical(sum(!is.na(st$adopt)), 41L)
  expect_identical(st$adopt[st$abb == "MI"], 2014L)   # April 2014
  expect_identical(st$adopt[st$abb == "NH"], 2015L)   # August 2014
  expect_identical(st$adopt[st$abb == "LA"], 2017L)   # July 2016
  expect_identical(st$expdate[st$abb == "VA"], as.Date("2019-01-01"))
  expect_true(all(is.na(st$adopt[st$abb %in% c("TX", "FL", "WI")])))
  expect_identical(medicaid$expanded,
                   as.integer(!is.na(medicaid$adopt) &
                                medicaid$year >= medicaid$adopt))
})

test_that("medicaid expansion lowers low-income uninsurance", {
  data(medicaid, package = "panglm", envir = environment())
  cs <- attgt(uninspov ~ 1, data = medicaid, index = c("state", "year"),
              treat = "adopt", control = "never")
  ov <- summary(cs)$overall
  expect_true(ov$estimate < 0 && ov$estimate + 1.96 * ov$se < 0)
  expect_equal(ov$estimate, -6.54, tolerance = 0.01)
})
