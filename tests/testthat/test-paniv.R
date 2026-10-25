iv_crime_formula <- function(over = FALSE) {
  exog <- "lprbconv + lprbpris + lavgsen + ldensity + lwcon + lpctymle"
  if (over) {
    stats::as.formula(paste("lcrmrte ~ lpolpc + lprbarr +", exog,
                            "| ltaxpc + lmix + lprbarr +", exog))
  } else {
    stats::as.formula(paste("lcrmrte ~ lprbarr + lpolpc +", exog, "| ltaxpc + lmix +", exog))
  }
}

iv_sim <- function(N = 300, TT = 5, seed = 1) {
  set.seed(seed)
  d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
  alpha <- rep(rnorm(N), each = TT)
  u <- rnorm(N * TT)
  d$z <- rnorm(N * TT)
  d$w <- rnorm(N * TT)
  d$x <- d$z + 0.8 * u + 0.5 * alpha + rnorm(N * TT)
  d$y <- 1 + 0.5 * d$x + 0.3 * d$w + alpha + u
  d
}

test_that("FE-2SLS (one-way and two-way) matches plm", {
  skip_if_missing("plm")
  data("Crime", package = "plm")
  for (over in c(FALSE, TRUE)) {
    f <- iv_crime_formula(over)
    for (eff in c("individual", "twoways")) {
      a <- paniv(f, Crime, index = c("county", "year"), model = "within", effect = eff)
      b <- plm::plm(f, Crime, index = c("county", "year"), model = "within", effect = eff)
      nm <- names(coef(a))
      expect_equal(coef(a), coef(b)[nm], tolerance = 1e-8)
      expect_equal(sqrt(diag(vcov(a))), sqrt(diag(vcov(b)))[nm], tolerance = 1e-8)
      vb <- plm::vcovHC(b, method = "arellano", type = "sss", cluster = "group")
      expect_equal(sqrt(diag(vcov(a, type = "cluster"))), sqrt(diag(vb))[nm], tolerance = 1e-8)
    }
  }
})

test_that("pooled 2SLS matches plm and AER::ivreg, with first-stage F and Sargan", {
  skip_if_missing("plm")
  skip_if_missing("AER")
  data("Crime", package = "plm")
  f <- iv_crime_formula(over = TRUE)
  a <- paniv(f, Crime, index = c("county", "year"), model = "pooling")
  b <- plm::plm(f, Crime, index = c("county", "year"), model = "pooling")
  iv <- AER::ivreg(f, data = Crime)
  expect_equal(coef(a), coef(b), tolerance = 1e-8)
  expect_equal(unname(coef(a)), unname(coef(iv)), tolerance = 1e-8)
  expect_equal(unname(sqrt(diag(vcov(a)))), unname(sqrt(diag(vcov(iv)))), tolerance = 1e-8)
  dg <- summary(iv, diagnostics = TRUE)$diagnostics
  expect_equal(a$first_stage$F, unname(dg[grepl("^Weak", rownames(dg)), "statistic"]),
               tolerance = 1e-8)
  expect_equal(a$overid$statistic, unname(dg["Sargan", "statistic"]), tolerance = 1e-8)
  expect_equal(a$overid$df, 1)
  expect_identical(a$endogenous, "lpolpc")
})

test_that("random-effects G2SLS and EC2SLS match plm, balanced and unbalanced", {
  skip_if_missing("plm")
  data("Crime", package = "plm")
  set.seed(5)
  CrimeU <- Crime[-sample(nrow(Crime), 60), ]
  f <- iv_crime_formula()
  for (d in list(Crime, CrimeU)) {
    for (im in c("bvk", "baltagi")) {
      a <- paniv(f, d, index = c("county", "year"), model = "random", inst.method = im)
      b <- plm::plm(f, d, index = c("county", "year"), model = "random", inst.method = im)
      expect_equal(coef(a), coef(b), tolerance = 1e-8)
      expect_equal(sqrt(diag(vcov(a))), sqrt(diag(vcov(b))), tolerance = 1e-8)
      expect_equal(unname(a$ercomp$sigma2), unname(b$ercomp$sigma2[c("idios", "id")]),
                   tolerance = 1e-8)
    }
  }
  a <- paniv(f, Crime, index = c("county", "year"), model = "random")
  b <- plm::plm(f, Crime, index = c("county", "year"), model = "random")
  vb <- plm::vcovHC(b, method = "arellano", type = "sss", cluster = "group")
  expect_equal(sqrt(diag(vcov(a, type = "cluster"))), sqrt(diag(vb)), tolerance = 1e-8)
})

test_that("Hausman-Taylor matches plm::pht, balanced and unbalanced", {
  skip_if_missing("plm")
  data("Wages", package = "plm")
  Wages$id <- rep(1:595, each = 7)
  Wages$t <- rep(1:7, 595)
  f <- lwage ~ wks + south + smsa + married + exp + I(exp^2) + bluecol + ind + union +
    sex + black + ed | sex + black + bluecol + south + smsa + ind
  # pht() re-evaluates its call as plm(...) in the calling frame
  plm <- plm::plm
  set.seed(7)
  for (d in list(Wages, Wages[-sample(nrow(Wages), 600), ])) {
    a <- paniv(f, d, index = c("id", "t"), model = "hausmantaylor")
    b <- suppressWarnings(plm::pht(f, data = d, index = c("id", "t")))
    nm <- names(coef(a))
    expect_equal(coef(a), coef(b)[nm], tolerance = 1e-8)
    expect_equal(sqrt(diag(vcov(a))), sqrt(diag(vcov(b)))[nm], tolerance = 1e-8)
    expect_equal(unname(a$ercomp$sigma2[c("idios", "id")]),
                 unname(unlist(b$ercomp$sigma2)[c("idios", "id")]), tolerance = 1e-8)
  }
  expect_identical(a$varlist$nc, "ed")
  expect_identical(a$varlist$xc, c("sexfemale", "blackyes"))
  expect_false(is.null(a$overid))
})

test_that("FE-2SLS removes the endogeneity bias of the within estimator", {
  est <- vapply(1:20, function(r) {
    d <- iv_sim(N = 200, TT = 4, seed = r)
    fe <- panglm(y ~ x + w, data = d, index = c("id", "t"), model = "within")
    iv <- paniv(y ~ x + w | z + w, data = d, index = c("id", "t"), model = "within")
    c(fe = unname(coef(fe)["x"]), iv = unname(coef(iv)["x"]))
  }, numeric(2))
  bias <- rowMeans(est) - 0.5
  expect_gt(bias[["fe"]], 0.15)
  expect_lt(abs(bias[["iv"]]), 0.03)
})

test_that("diagnostics, methods and the cluster variants are consistent", {
  d <- iv_sim()
  d$z2 <- d$z + rnorm(nrow(d))
  fit <- paniv(y ~ x + w | z + z2 + w, data = d, index = c("id", "t"), model = "within",
               vcov = "cluster")
  expect_s3_class(fit, "paniv")
  expect_identical(fit$vcov_type, "cluster")
  expect_equal(vcov(fit), vcov(fit, type = "cluster"))
  expect_false(isTRUE(all.equal(vcov(fit), vcov(fit, type = "classical"))))
  expect_identical(fit$overid$method, "Hansen J")
  expect_equal(fit$overid$df, 1)
  expect_true(all(c("partial.R2", "F", "F.robust") %in% names(fit$first_stage)))
  expect_gt(fit$first_stage$F, 100)
  expect_equal(nobs(fit), nrow(d))
  ci <- confint(fit)
  expect_equal(unname(ci[, 1]), unname(coef(fit) - qnorm(0.975) * sqrt(diag(vcov(fit)))))
  expect_output(print(fit), "fixed-effects 2SLS")
  expect_output(print(summary(fit)), "Hansen J")
  # a user cluster equal to the individual reproduces the default
  expect_equal(vcov(fit, type = "cluster", cluster = d$id), vcov(fit, type = "cluster"))
  # the first stage is unaffected by row order
  shuffled <- d[sample(nrow(d)), ]
  fit2 <- paniv(y ~ x + w | z + z2 + w, data = shuffled, index = c("id", "t"),
                model = "within", vcov = "cluster")
  expect_equal(coef(fit2), coef(fit))
  expect_equal(vcov(fit2), vcov(fit))
  classical <- paniv(y ~ x + w | z + z2 + w, data = d, index = c("id", "t"), model = "within")
  expect_identical(classical$overid$method, "Sargan")
  expect_gt(classical$overid$p.value, 0.001)
})

test_that("missing values are dropped and an exactly identified model has no overid test", {
  d <- iv_sim(N = 50)
  d$x[3] <- NA
  d$z[10] <- NA
  fit <- paniv(y ~ x + w | z + w, data = d, index = c("id", "t"), model = "pooling")
  expect_equal(nobs(fit), nrow(d) - 2)
  expect_null(fit$overid)
  expect_length(fit$na.action, 2)
})

test_that("paniv() reports invalid input clearly", {
  d <- iv_sim(N = 40)
  idx <- c("id", "t")
  expect_error(paniv(y ~ x + w, data = d, index = idx), "two parts")
  expect_error(paniv(y ~ x | z | w, data = d, index = idx), "exactly two parts")
  expect_error(paniv(y ~ x | z, data = d), "'index' is required")
  expect_error(paniv(y ~ x | z, data = d, index = "nope"), "not found")
  expect_error(paniv(y ~ x + w | z, data = d, index = idx, model = "pooling"),
               "insufficient number of instruments")
  expect_error(paniv(y ~ x | z, data = d, index = idx, model = "random", effect = "twoways"),
               "only available for model = 'within'")
  d$g <- rep(rnorm(40), each = 5)
  expect_error(paniv(y ~ x | g, data = d, index = idx, model = "within"),
               "insufficient number of instruments")
  expect_warning(paniv(y ~ x + w | x + w, data = d, index = idx, model = "pooling"),
                 "no endogenous regressor")
  # Hausman-Taylor identification and formula conventions
  d$s <- rep(rnorm(40), each = 5)
  expect_error(paniv(y ~ x + s | 1, data = d, index = idx, model = "hausmantaylor"),
               "not identified")
  expect_error(paniv(y ~ x + s | z, data = d, index = idx, model = "hausmantaylor"),
               "subset of the first part")
  expect_error(vcov(paniv(y ~ x | z, data = d, index = idx), type = "HC1"))
})

test_that("random-effects 2SLS keeps time-invariant regressors without a spurious warning", {
  d <- iv_sim(N = 200, TT = 4, seed = 7)
  d$female <- rep(rbinom(200, 1, 0.5), each = 4)
  d$y <- d$y + 0.4 * d$female
  f <- y ~ x + w + female | z + w + female
  for (im in c("bvk", "baltagi")) {
    a <- expect_no_warning(paniv(f, d, index = c("id", "t"), model = "random",
                                 inst.method = im))
    expect_true("female" %in% names(coef(a)))
    skip_if_missing("plm")
    b <- plm::plm(f, d, index = c("id", "t"), model = "random", inst.method = im)
    expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-8)
    expect_equal(unname(a$ercomp$sigma2), unname(b$ercomp$sigma2[c("idios", "id")]),
                 tolerance = 1e-8)
  }
  # a regressor that is truly unidentified in the final fit still warns
  d$female2 <- 2 * d$female
  expect_warning(paniv(y ~ x + w + female + female2 | z + w + female + female2, d,
                       index = c("id", "t"), model = "random"),
                 "dropping unidentified column")
})
