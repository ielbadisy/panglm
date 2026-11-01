test_that("docvisits has the documented structure", {
  data(docvisits, package = "panglm", envir = environment())
  expect_s3_class(docvisits, "data.frame")
  expect_identical(dim(docvisits), c(27326L, 21L))
  expect_identical(names(docvisits),
                   c("id", "year", "docvis", "hospvis", "doctor", "hospital",
                     "hsat", "female", "age", "educ", "married", "hhkids",
                     "hhninc", "working", "bluecollar", "whitecollar", "self",
                     "civil", "handicap", "public", "addon"))
  expect_false(anyNA(docvisits))
  expect_false(anyDuplicated(docvisits[, c("id", "year")]) > 0)
  expect_identical(length(unique(docvisits$id)), 7293L)
  expect_identical(sort(unique(docvisits$year)),
                   c(1984L, 1985L, 1986L, 1987L, 1988L, 1991L, 1994L))
  ints <- setdiff(names(docvisits), c("educ", "hhninc"))
  for (v in ints) expect_type(docvisits[[v]], "integer")
  expect_type(docvisits$educ, "double")
  expect_type(docvisits$hhninc, "double")
})

test_that("docvisits values are in range and match the source", {
  data(docvisits, package = "panglm", envir = environment())
  bin <- c("doctor", "hospital", "female", "married", "hhkids", "working",
           "bluecollar", "whitecollar", "self", "civil", "handicap",
           "public", "addon")
  for (v in bin) expect_true(all(docvisits[[v]] %in% 0:1), info = v)
  expect_true(all(docvisits$hsat %in% 0:10))
  expect_identical(range(docvisits$age), c(25L, 64L))
  expect_identical(range(docvisits$docvis), c(0L, 121L))
  expect_identical(range(docvisits$hospvis), c(0L, 51L))
  expect_identical(docvisits$doctor, as.integer(docvisits$docvis > 0))
  expect_identical(docvisits$hospital, as.integer(docvisits$hospvis > 0))
  expect_equal(range(docvisits$educ), c(7, 18))
  expect_equal(mean(docvisits$docvis), 3.183525, tolerance = 1e-6)
  expect_equal(mean(docvisits$hhninc), 3.520836, tolerance = 1e-6)
  expect_equal(mean(docvisits$female), 0.4787748, tolerance = 1e-6)
  expect_identical(sum(docvisits$addon), 514L)
  expect_identical(nrow(docvisits[docvisits$year == 1988, ]), 4483L)
  # first person-year of the file: id 1 in 1984
  r1 <- docvisits[docvisits$id == 1L & docvisits$year == 1984L, ]
  expect_identical(c(r1$docvis, r1$hsat, r1$age, r1$female), c(1L, 8L, 54L, 0L))
  expect_equal(r1$hhninc, 3.05)
})
