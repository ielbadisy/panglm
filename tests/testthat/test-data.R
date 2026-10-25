test_that("copd has the 10 columns its help page documents", {
  data(copd, package = "panglm", envir = environment())
  expect_identical(names(copd), c("id", "visit", "treatment", "age", "smoker", "bmi",
                                  "crp", "fev1", "exacerbations", "hospitalized"))
  expect_identical(dim(copd), c(300L, 10L))
  db <- tryCatch(tools::Rd_db("panglm"), error = function(e) NULL)
  skip_if(is.null(db) || is.null(db[["copd.Rd"]]), "installed help pages not available")
  rd <- paste(as.character(db[["copd.Rd"]]), collapse = "")
  expect_match(rd, "and 10 columns", fixed = TRUE)
  for (v in names(copd)) expect_match(rd, paste0("\\item{", v, "}"), fixed = TRUE)
})
