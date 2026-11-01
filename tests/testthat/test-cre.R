test_that("Mundlak pooled and random fits reproduce the within slope (gaussian)", {
  data(copd, package = "panglm")
  fe <- panglm(fev1 ~ crp, copd, c("id", "visit"), model = "within")
  po <- panglm(fev1 ~ treatment + crp, copd, c("id", "visit"), model = "pooling", cre = TRUE)
  re <- panglm(fev1 ~ treatment + crp, copd, c("id", "visit"), model = "random", cre = TRUE)
  expect_equal(unname(coef(po)["crp"]), unname(coef(fe)), tolerance = 1e-8)
  expect_equal(unname(coef(re)["crp"]), unname(coef(fe)), tolerance = 1e-6)
  expect_true("mean_crp" %in% names(coef(po)))
  expect_false("mean_treatment" %in% names(coef(po)))
})

test_that("CRE works for nonlinear families and reports a Mundlak test", {
  data(copd, package = "panglm")
  f <- panglm(hospitalized ~ treatment + age + crp, copd, c("id", "visit"),
              model = "random", family = "probit", cre = TRUE)
  s <- summary(f)
  expect_false(is.null(s$cre_test))
  expect_true(s$cre_test$p.value >= 0 && s$cre_test$p.value <= 1)
  expect_output(print(s), "Mundlak")
  expect_equal(unname(predict(f, copd, type = "response")[order(copd$id, copd$visit)]),
               unname(fitted(f)), tolerance = 1e-8)
  expect_error(panglm(fev1 ~ crp, copd, c("id", "visit"), model = "within", cre = TRUE))
})

test_that("CRE drops a mean that is collinear up to rounding (age rising one per year)", {
  # age = birth cohort + year with a rounding error of 2e-6 (as in
  # pglm::HealthIns, stored in single precision): its within part equals the
  # year trend up to rounding, so the Mundlak means of the year dummies and
  # of age are dependent up to rounding
  set.seed(11)
  G <- 80
  d <- do.call(rbind, lapply(seq_len(G), function(i) {
    years <- sort(sample(1:5, sample(2:5, 1)))
    data.frame(id = i, year = years)
  }))
  d$age <- rep(runif(G, 30, 60), table(d$id)) + d$year + runif(nrow(d), -2e-6, 2e-6)
  d$x <- rnorm(nrow(d))
  d$y <- 1 + 0.5 * d$x + 0.02 * d$age + rep(rnorm(G), table(d$id)) + rnorm(nrow(d))
  d$fy <- factor(d$year)
  expect_warning(
    fit <- panglm(y ~ age + fy + x, d, c("id", "year"), model = "random", cre = TRUE),
    "dropping unidentified column\\(s\\): mean_fy5$")
  expect_identical(fit$aliased, "mean_fy5")
  expect_true(all(abs(coef(fit)) < 10))
  # pooled CRE is least squares on the reduced Mundlak design
  expect_warning(
    po <- panglm(y ~ age + fy + x, d, c("id", "year"), model = "pooling", cre = TRUE),
    "mean_fy5")
  X <- model.matrix(~ age + fy + x, d)
  M <- apply(X[, -1], 2, ave, d$id)
  colnames(M) <- paste0("mean_", colnames(M))
  ref <- stats::lm.fit(cbind(X, M[, colnames(M) != "mean_fy5"]), d$y)$coefficients
  expect_equal(unname(coef(po)), unname(ref), tolerance = 1e-6)
})

test_that("exact collinearity drops the same columns as lm()", {
  set.seed(3)
  d <- data.frame(id = rep(1:30, each = 4), year = rep(1:4, 30))
  d$age <- rep(sample(30:60, 30, TRUE), each = 4) + d$year
  d$x <- rnorm(120)
  d$y <- d$x + rnorm(120)
  X <- model.matrix(~ age + factor(year) + x, d)
  M <- apply(X[, -1], 2, ave, d$id)
  colnames(M) <- paste0("mean_", colnames(M))
  na_lm <- names(which(is.na(stats::lm.fit(cbind(X, M), d$y)$coefficients)))
  fit <- suppressWarnings(panglm(y ~ age + factor(year) + x, d, c("id", "year"),
                                 model = "random", cre = TRUE))
  expect_setequal(fit$aliased, na_lm)
})
