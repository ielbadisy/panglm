tidy_cols <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high")

expect_tidy_matches <- function(td, b, V, ci) {
  expect_true(all(tidy_cols %in% names(td)))
  expect_equal(td$term, names(b))
  expect_equal(td$estimate, unname(b), tolerance = 1e-12)
  expect_equal(td$std.error, unname(sqrt(diag(V))), tolerance = 1e-12)
  expect_equal(td$statistic, td$estimate / td$std.error, tolerance = 1e-12)
  expect_equal(td$p.value, 2 * pnorm(-abs(td$statistic)), tolerance = 1e-12)
  expect_equal(td$conf.low, unname(ci[, 1]), tolerance = 1e-10)
  expect_equal(td$conf.high, unname(ci[, 2]), tolerance = 1e-10)
}

test_that("tidy agrees with coef, vcov and confint for the coefficient classes", {
  skip_if_missing("generics")
  for (nm in c("pangmm", "paniv", "panselect", "panfmm", "twopart")) {
    fit <- methods_fit(nm)
    td <- generics::tidy(fit, conf.int = TRUE, conf.level = 0.9)
    expect_tidy_matches(td, coef(fit), vcov(fit), confint(fit, level = 0.9))
    expect_false(any(c("conf.low", "conf.high") %in% names(generics::tidy(fit))))
    gl <- generics::glance(fit)
    expect_s3_class(gl, "data.frame")
    expect_equal(nrow(gl), 1L)
    expect_equal(gl$nobs, nobs(fit))
  }
})

test_that("tidy adds the block columns of twopart and panfmm", {
  skip_if_missing("generics")
  td <- generics::tidy(methods_fit("twopart"))
  expect_setequal(td$component, c("binary", "positive"))
  expect_true(all(startsWith(td$term, paste0(td$component, ":"))))
  gl <- generics::glance(methods_fit("twopart"))
  expect_equal(gl$logLik, as.numeric(logLik(methods_fit("twopart"))))
  expect_equal(gl$AIC, AIC(methods_fit("twopart")))
  expect_equal(gl$BIC, BIC(methods_fit("twopart")))

  fm <- methods_fit("panfmm")
  td <- generics::tidy(fm)
  expect_equal(td$class, c(1L, 1L, 2L, 2L, 2L))
  expect_equal(td$component, c("class", "class", "class", "class", "prior"))
  gl <- generics::glance(fm)
  expect_equal(gl$logLik, as.numeric(logLik(fm)))
  expect_equal(gl$AIC, AIC(fm))
  expect_equal(gl$k, 2L)
})

test_that("glance of pangmm, paniv and panselect reports the specification tests", {
  skip_if_missing("generics")
  g <- methods_fit("pangmm")
  gl <- generics::glance(g)
  expect_equal(gl$hansen, g$hansen$statistic)
  expect_equal(gl$sargan.p.value, g$sargan$p.value)
  expect_equal(gl$ar2.p.value, g$artest$p.value[2])
  expect_equal(gl$n.instruments, g$n_instruments)
  iv <- methods_fit("paniv")
  expect_equal(generics::glance(iv)$overid, iv$overid$statistic)
  ps <- methods_fit("panselect")
  expect_equal(generics::glance(ps)$selection.p.value, ps$selection_test$cluster$p.value)
})

test_that("tidy.panqr has one row per term and quantile", {
  skip_if_missing("generics")
  fit <- methods_fit("panqr")
  td <- generics::tidy(fit, conf.int = TRUE)
  expect_equal(nrow(td), length(coef(fit)))
  expect_equal(td$tau, rep(fit$tau, each = nrow(fit$coefficients)))
  for (j in seq_along(fit$tau)) {
    r <- td[td$tau == fit$tau[j], ]
    expect_equal(r$estimate, unname(fit$coefficients[, j]))
    expect_equal(r$std.error, unname(sqrt(diag(fit$vcov[[j]]))))
    ci <- confint(fit, tau = fit$tau[j])
    expect_equal(r$conf.low, unname(ci[, 1]))
    expect_equal(r$conf.high, unname(ci[, 2]))
  }
  gl <- generics::glance(fit)
  expect_equal(gl$n.tau, 3L)
  expect_equal(gl$nobs, nobs(fit))
})

test_that("tidy.panglm_attgt gives the group-time effects and the aggregations", {
  skip_if_missing("generics")
  fit <- methods_fit("attgt")
  td <- generics::tidy(fit, conf.int = TRUE)
  expect_equal(td$term, names(coef(fit)))
  expect_equal(td$estimate, unname(coef(fit)))
  expect_equal(td$std.error, fit$estimates$se)
  expect_equal(td$conf.low, fit$estimates$lower, tolerance = 1e-10)
  expect_equal(td$event, td$time - td$group)
  s <- summary(fit, type = "dynamic")
  tdd <- generics::tidy(fit, type = "dynamic", conf.int = TRUE)
  expect_equal(tdd$event, s$estimates$level)
  expect_equal(tdd$estimate, s$estimates$estimate)
  expect_equal(tdd$conf.high, s$estimates$upper, tolerance = 1e-10)
  tds <- generics::tidy(fit, type = "simple")
  expect_equal(tds$term, "ATT")
  expect_equal(tds$estimate, summary(fit)$overall$estimate)
  expect_equal(generics::tidy(fit, type = "group")$group, summary(fit, type = "group")$estimates$level)
  gl <- generics::glance(fit)
  expect_equal(gl$nobs, nobs(fit))
  expect_equal(gl$pretrend.p.value, fit$pretrend$p.value)
})

test_that("tidy.panglm_etwfe gives the ATT aggregations and the coefficients", {
  skip_if_missing("generics")
  fit <- methods_fit("etwfe")
  s <- summary(fit)
  td <- generics::tidy(fit, conf.int = TRUE)
  expect_equal(td$term, "ATT")
  expect_equal(td$estimate, s$overall$estimate)
  expect_equal(td$p.value, s$overall$p.value)
  expect_equal(td$conf.low, s$overall$lower)
  se <- summary(fit, type = "event")$estimates
  tde <- generics::tidy(fit, type = "event", conf.int = TRUE)
  expect_equal(tde$event, se$level)
  expect_equal(tde$statistic, se$statistic)
  expect_equal(tde$conf.high, se$upper)
  tdc <- generics::tidy(fit, type = "coefficients", conf.int = TRUE)
  expect_equal(tdc$estimate, unname(coef(fit)))
  expect_equal(tdc$std.error, unname(sqrt(diag(vcov(fit)))))
  expect_equal(generics::glance(fit)$nobs, nobs(fit))
})

test_that("tidy.panglm_eventstudy keeps the reference period", {
  skip_if_missing("generics")
  fit <- methods_fit("eventstudy")
  td <- generics::tidy(fit, conf.int = TRUE)
  expect_equal(td$event, fit$estimates$event)
  expect_equal(td$estimate, fit$estimates$estimate)
  expect_equal(td$conf.low, fit$estimates$lower, tolerance = 1e-10)
  ref <- td$event == fit$ref
  expect_true(is.na(td$statistic[ref]))
  expect_equal(td$estimate[ref], 0)
  gl <- generics::glance(fit)
  expect_equal(gl$pretrend.statistic, fit$pretrend$statistic)
})

test_that("tidy and glance of the test objects are one-row htest-like frames", {
  skip_if_missing("generics")
  copd <- methods_data("copd")
  fe <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"), model = "within")
  re <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"), model = "random")
  h <- panhausman(fe, re)
  td <- generics::tidy(h)
  expect_equal(nrow(td), 1L)
  expect_equal(td$statistic, h$statistic)
  expect_equal(td$p.value, h$p.value)
  expect_equal(td$parameter, as.numeric(h$parameter))
  expect_identical(generics::glance(h), td)
  hm <- panhausman(fe, method = "mundlak", data = copd)
  tdm <- generics::tidy(hm)
  expect_equal(c(tdm$parameter, tdm$den.df), as.numeric(hm$parameter))
  expect_equal(tdm$distribution, "F")

  jt <- jointtest(re, terms = "crp")
  expect_equal(generics::tidy(jt)$statistic, jt$statistic)
  expect_equal(generics::tidy(jt)$parameter, jt$df)

  pois <- panglm(exacerbations ~ crp, data = copd, index = c("id", "visit"),
                 model = "pooling", family = "poisson")
  dt <- pandispersiontest(pois)
  tdd <- generics::tidy(dt)
  expect_equal(tdd$statistic, dt$statistic)
  expect_equal(tdd$ratio, dt$ratio)
  expect_equal(tdd$p.value, dt$p.value)
})

test_that("modelsummary tabulates the main classes", {
  skip_if_missing("generics")
  skip_if_missing("modelsummary")
  copd <- methods_data("copd")
  models <- list(
    panglm = panglm(fev1 ~ treatment + age + crp, data = copd, index = c("id", "visit"),
                    model = "random"),
    twopart = methods_fit("twopart"),
    pangmm = methods_fit("pangmm"),
    paniv = methods_fit("paniv"),
    panselect = methods_fit("panselect"),
    panfmm = methods_fit("panfmm"),
    etwfe = methods_fit("etwfe"),
    eventstudy = methods_fit("eventstudy")
  )
  tab <- modelsummary::modelsummary(models, output = "data.frame")
  expect_s3_class(tab, "data.frame")
  expect_true(all(names(models) %in% names(tab)))
  tq <- modelsummary::modelsummary(methods_fit("panqr"), output = "data.frame",
                                   shape = term + tau ~ model)
  expect_s3_class(tq, "data.frame")
})

test_that("glance() adjusts R-squared for the absorbed fixed effects", {
  data(copd, package = "panglm", envir = environment())
  pooled <- panglm(fev1 ~ crp + age, data = copd, index = c("id", "visit"),
                   model = "pooling", family = "gaussian")
  ref <- summary(stats::lm(fev1 ~ crp + age, data = copd))
  expect_equal(generics::glance(pooled)$r.squared, ref$r.squared, tolerance = 1e-10)
  expect_equal(generics::glance(pooled)$adj.r.squared, ref$adj.r.squared, tolerance = 1e-10)
  within <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
                   model = "within", family = "gaussian")
  g <- generics::glance(within)
  expect_lt(g$adj.r.squared, g$r.squared)
  dummies <- summary(stats::lm(fev1 ~ crp + factor(id), data = copd))
  expect_equal(g$r.squared, dummies$r.squared, tolerance = 1e-8)
  expect_equal(g$adj.r.squared, dummies$adj.r.squared, tolerance = 1e-8)
  twoway <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
                   model = "within", effect = "twoways", family = "gaussian")
  dummies2 <- summary(stats::lm(fev1 ~ crp + factor(id) + factor(visit), data = copd))
  expect_equal(generics::glance(twoway)$r.squared, dummies2$r.squared, tolerance = 1e-6)
  expect_equal(generics::glance(twoway)$adj.r.squared, dummies2$adj.r.squared,
               tolerance = 1e-6)
  skip_if_not_installed("fixest")
  fx <- fixest::feols(fev1 ~ crp | id, data = copd)
  expect_equal(g$adj.r.squared, unname(fixest::r2(fx, "ar2")), tolerance = 1e-8)
})
