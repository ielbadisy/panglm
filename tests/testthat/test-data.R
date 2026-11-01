test_that("copd has the 16 columns its help page documents", {
  data(copd, package = "panglm", envir = environment())
  expect_identical(names(copd), c("id", "visit", "treatment", "age", "smoker", "bmi",
                                  "crp", "fev1", "exacerbations", "hospitalized",
                                  "clinic", "adopt", "programme", "cost", "adherence",
                                  "symptoms"))
  expect_identical(dim(copd), c(300L, 16L))
  db <- tryCatch(tools::Rd_db("panglm"), error = function(e) NULL)
  skip_if(is.null(db) || is.null(db[["copd.Rd"]]), "installed help pages not available")
  rd <- paste(as.character(db[["copd.Rd"]]), collapse = "")
  expect_match(rd, "and 16 columns", fixed = TRUE)
  for (v in names(copd)) expect_match(rd, paste0("\\item{", v, "}"), fixed = TRUE)
})

test_that("copd clinical columns are unchanged from earlier releases", {
  data(copd, package = "panglm", envir = environment())
  expect_false(anyDuplicated(copd[c("id", "visit")]) > 0)
  expect_identical(sum(copd$exacerbations), 640L)
  expect_identical(sum(copd$hospitalized), 79L)
  expect_equal(mean(copd$fev1), 2.947533, tolerance = 1e-6)
  expect_equal(sum(copd$crp), 1548.3, tolerance = 1e-8)
})

test_that("copd economic columns follow their documented design", {
  data(copd, package = "panglm", envir = environment())
  expect_true(all(tapply(copd$clinic, copd$id, function(x) length(unique(x))) == 1))
  expect_identical(as.integer(table(copd$clinic[copd$visit == 1])), rep(10L, 6))
  adopt <- tapply(copd$adopt, copd$clinic, unique)
  expect_identical(as.integer(adopt), c(3L, 3L, 4L, 4L, 0L, 0L))
  expect_identical(copd$programme,
                   as.integer(copd$adopt > 0 & copd$visit >= copd$adopt))
  zero <- mean(copd$cost == 0)
  expect_true(zero > 0.25 && zero < 0.35)
  expect_true(all(copd$cost >= 0))
  expect_true(all(copd$adherence > 0 & copd$adherence < 1))
  expect_true(is.ordered(copd$symptoms))
  expect_identical(levels(copd$symptoms), c("none", "mild", "moderate", "severe"))
})

test_that("estimators recover the copd generating values within sampling error", {
  data(copd, package = "panglm", envir = environment())
  idx <- c("id", "visit")
  near <- function(est, se, truth) expect_lt(abs(est - truth), 3 * se)

  tp <- twopart(cost ~ exacerbations + hospitalized + smoker + treatment + crp +
                  visit + programme + factor(clinic),
                binary_formula = ~ exacerbations + hospitalized + smoker + treatment + crp,
                data = copd, index = idx, positive = "lognormal")
  b <- coef(tp$positive)
  s <- sqrt(diag(vcov(tp$positive)))
  near(b[["programme"]], s[["programme"]], -0.50)
  near(b[["hospitalized"]], s[["hospitalized"]], 0.70)
  near(b[["exacerbations"]], s[["exacerbations"]], 0.10)

  sy <- panglm(symptoms ~ crp, data = copd, index = idx, model = "within",
               family = "ologit")
  near(coef(sy)[["crp"]], sqrt(vcov(sy)["crp", "crp"]), 0.25)

  att <- summary(attgt(cost ~ 1, data = copd, index = idx, treat = "adopt",
                       control = "notyet"))$overall
  near(att$estimate, att$se, -406.0)  # sample ATT documented in ?copd
})
