test_that("meps has the documented structure", {
  data(meps, package = "panglm", envir = environment())
  expect_s3_class(meps, "data.frame")
  expect_identical(dim(meps), c(6000L, 20L))
  expect_identical(names(meps), c(
    "id", "year", "totexp", "oopexp", "anyexp", "officevis", "ervis",
    "ipstays", "insurance", "age", "sex", "race", "region", "povcat",
    "povlev", "health", "nchronic", "longwt", "varstr", "varpsu"))
  expect_false(anyDuplicated(meps[c("id", "year")]) > 0)
  expect_identical(as.vector(table(table(meps$id))), 3000L)
  expect_true(all(table(meps$id) == 2L))
  for (v in c("id", "year", "anyexp", "officevis", "ervis", "ipstays", "age",
              "nchronic", "varstr", "varpsu")) {
    expect_type(meps[[v]], "integer")
  }
  for (v in c("insurance", "sex", "race", "region", "povcat")) {
    expect_true(is.factor(meps[[v]]))
  }
  expect_true(is.ordered(meps$health))
  expect_identical(levels(meps$insurance), c("private", "public", "uninsured"))
})

test_that("meps values are in range and match the build", {
  data(meps, package = "panglm", envir = environment())
  expect_setequal(unique(meps$year), c(2018L, 2019L))
  expect_true(all(meps$totexp >= 0) && all(meps$oopexp >= 0))
  expect_true(all(meps$oopexp <= meps$totexp))
  expect_identical(meps$anyexp, as.integer(meps$totexp > 0))
  expect_true(all(meps$age >= 18 & meps$age <= 85))
  expect_true(all(meps$longwt > 0))
  expect_true(all(meps$nchronic >= 0 & meps$nchronic <= 8))
  # longitudinal weight and design variables are constant within person
  expect_true(all(tapply(meps$longwt, meps$id, function(x) diff(range(x))) == 0))
  expect_equal(sum(is.na(meps$health)), 4L)
  expect_equal(mean(meps$totexp), 7698.3327, tolerance = 1e-5)
  expect_equal(mean(meps$officevis), 7.6312, tolerance = 1e-5)
  expect_equal(mean(meps$anyexp), 0.8585, tolerance = 1e-5)
  expect_equal(sum(meps$insurance == "uninsured"), 487L)
  first <- meps[meps$id == 1L, ]
  expect_identical(first$totexp, c(0, 54))
  expect_identical(first$age, c(41L, 42L))
})
