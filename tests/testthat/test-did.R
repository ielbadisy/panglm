did_panel <- function(seed, N = 240, dyn = function(e) 1) {
  set.seed(seed)
  d <- expand.grid(id = seq_len(N), year = 2011:2020)
  coh <- sample(c(2014, 2016, 2018, NA), N, replace = TRUE)
  d$first <- coh[d$id]
  a <- rnorm(N)
  d$x <- rnorm(nrow(d))
  e <- d$year - d$first
  d$post <- as.numeric(!is.na(e) & e >= 0)
  d$y <- a[d$id] + 0.1 * (d$year - 2011) + 0.3 * d$x +
    ifelse(d$post == 1, dyn(e), 0) + rnorm(nrow(d))
  d$cnt <- rpois(nrow(d), exp(0.2 * a[d$id] + 0.2 * d$x + 0.3 * d$post))
  d$nb <- rnbinom(nrow(d), mu = exp(0.5 + 0.5 * a[d$id] + 0.2 * d$x + 0.3 * d$post), size = 2)
  d[sample(nrow(d)), ]
}

test_that("did() is the TWFE regression on the treatment indicator", {
  d <- did_panel(1)
  fit <- did(y ~ x, data = d, index = c("id", "year"), treat = "first")
  ref <- panglm(y ~ post + x, data = d, index = c("id", "year"), model = "within",
                effect = "twoways", vcov = "cluster")
  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-10)
  expect_equal(unname(vcov(fit)), unname(vcov(ref)), tolerance = 1e-10)
  expect_equal(fit$did$cohorts, c(2014, 2016, 2018))
  expect_equal(fit$did$treated_units + fit$did$control_units, 240)
  expect_true(abs(coef(fit)[["treated"]] - 1) < 0.3)
  # y ~ 1 works, and 0 / Inf code never-treated units
  d2 <- d; d2$first[is.na(d2$first)] <- Inf
  expect_equal(coef(did(y ~ 1, d2, c("id", "year"), "first")),
               coef(did(y ~ 1, d, c("id", "year"), "first")))
})

test_that("did() and eventstudy() agree with fixest", {
  skip_if_not_installed("fixest")
  d <- did_panel(2, dyn = function(e) 0.5 + 0.2 * e)
  d$treated2 <- d$post
  fx <- fixest::feols(y ~ treated2 + x | id + year, d, cluster = ~id)
  fit <- did(y ~ x, data = d, index = c("id", "year"), treat = "first")
  expect_equal(unname(coef(fit)), unname(coef(fx)), tolerance = 1e-8)
  # cluster SE equal up to fixest's small-sample factor
  expect_equal(unname(sqrt(diag(vcov(fit)))), as.numeric(fixest::se(fx)), tolerance = 0.01)

  es <- eventstudy(y ~ x, d, c("id", "year"), "first", leads = 3, lags = 3)
  d$rel <- pmin(pmax(d$year - d$first, -3), 3)
  d$rel[is.na(d$rel)] <- -1000
  ex <- fixest::feols(y ~ i(rel, ref = c(-1, -1000)) + x | id + year, d, cluster = ~id)
  tab <- es$estimates[es$estimates$event != -1, ]
  expect_equal(tab$estimate, unname(coef(ex)[1:6]), tolerance = 1e-8)
  expect_equal(tab$se, as.numeric(fixest::se(ex)[1:6]), tolerance = 0.01)
  w <- fixest::wald(ex, keep = "rel::-[23]", print = FALSE)
  expect_equal(es$pretrend$statistic, 2 * w$stat, tolerance = 0.01)
  expect_equal(es$pretrend$df, 2)
  expect_equal(coef(es)[["-1"]], 0)
})

test_that("eventstudy() detects a violated pre-trend and supports bin = FALSE", {
  set.seed(3)
  N <- 300
  d <- expand.grid(id = seq_len(N), year = 1:10)
  d$first <- ifelse(d$id <= 150, 7, NA)
  a <- rnorm(N)
  # treated units trend upwards before treatment: parallel trends fails
  d$y <- a[d$id] + ifelse(d$id <= 150, 0.3 * d$year, 0) + rnorm(nrow(d), sd = 0.5)
  es <- eventstudy(y ~ 1, d, c("id", "year"), "first", leads = 4, lags = 2)
  expect_lt(es$pretrend$p.value, 1e-6)
  es2 <- eventstudy(y ~ 1, d, c("id", "year"), "first", leads = 2, lags = 2, bin = FALSE)
  expect_equal(nobs(es2$fit), nrow(d) - 150 * 5)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(es))
  expect_output(print(es), "Pre-trend test")
})

test_that("two-way Poisson and negbin robust vcov profile out the time effects", {
  d <- did_panel(4)
  skip_if_not_installed("fixest")
  fp <- panglm(cnt ~ post + x, d, c("id", "year"), model = "within", family = "poisson",
               effect = "twoways", vcov = "cluster")
  fx <- fixest::fepois(cnt ~ post + x | id + year, d, cluster = ~id)
  expect_equal(unname(coef(fp)), unname(coef(fx)), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(fp)))), as.numeric(fixest::se(fx)), tolerance = 0.01)
  skip_if_not_installed("sandwich")
  skip_if_not_installed("MASS")
  fn <- panglm(nb ~ post + x, d, c("id", "year"), model = "within", family = "negbin",
               effect = "twoways", vcov = "cluster")
  g <- MASS::glm.nb(nb ~ post + x + factor(id) + factor(year), d)
  expect_equal(unname(coef(fn)), unname(coef(g)[2:3]), tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(vcov(fn)))),
               unname(sqrt(diag(sandwich::vcovCL(g, cluster = ~id, type = "HC0")))[2:3]),
               tolerance = 0.01)
  expect_true(all(is.finite(vcov(fn, type = "HC1"))))
})

test_that("did() and eventstudy() with a binary outcome are the conditional logit", {
  set.seed(6)
  N <- 300
  d <- expand.grid(id = seq_len(N), year = 2011:2018)
  d$first <- sample(c(2014, 2016, NA), N, replace = TRUE)[d$id]
  a <- rnorm(N)
  d$x <- rnorm(nrow(d))
  d$z <- rnorm(N)[d$id]
  d$post <- as.numeric(!is.na(d$first) & d$year >= d$first)
  d$yb <- rbinom(nrow(d), 1, plogis(a[d$id] - 0.2 + 0.1 * (d$year - 2011) + 0.4 * d$x +
                                     0.7 * d$post))
  fit <- suppressMessages(did(yb ~ x, d, c("id", "year"), "first", family = "binomial"))
  ref <- suppressMessages(panglm(yb ~ post + x, d, c("id", "year"), model = "within",
                                 family = "binomial", effect = "twoways", vcov = "cluster"))
  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-10)
  expect_equal(unname(vcov(fit)), unname(vcov(ref)), tolerance = 1e-10)
  expect_identical(fit$family$family, "binomial")
  expect_lt(abs(coef(fit)[["treated"]] - 0.7), 0.35)
  # a time-invariant control is absorbed by the individual effects
  expect_warning(f2 <- suppressMessages(did(yb ~ x + z, d, c("id", "year"), "first",
                                            family = "binomial")), "z")
  expect_equal(coef(f2), coef(fit), tolerance = 1e-10)
  skip_if_not_installed("survival")
  suppressMessages(library(survival))
  cl <- clogit(yb ~ post + x + factor(year) + strata(id), d, method = "exact")
  fc <- suppressMessages(did(yb ~ x, d, c("id", "year"), "first", family = "binomial",
                             vcov = "classical"))
  expect_equal(unname(coef(fc)), unname(coef(cl)[1:2]), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(fc)))), unname(sqrt(diag(vcov(cl)))[1:2]),
               tolerance = 1e-6)
  es <- suppressMessages(eventstudy(yb ~ x, d, c("id", "year"), "first", leads = 2, lags = 2,
                                    family = "binomial", vcov = "classical"))
  d$rel <- pmin(pmax(d$year - d$first, -2), 2)
  d$rel[is.na(d$rel)] <- -1
  cle <- clogit(yb ~ relevel(factor(rel), "-1") + x + factor(year) + strata(id), d,
                method = "exact")
  tab <- es$estimates[es$estimates$event != -1, ]
  expect_equal(tab$estimate, unname(coef(cle)[1:4]), tolerance = 1e-6)
  expect_equal(tab$se, unname(sqrt(diag(vcov(cle)))[1:4]), tolerance = 1e-6)
  expect_equal(es$pretrend$df, 1)
  expect_output(print(es), "binomial")
})

test_that("did() and eventstudy() reject bad input", {
  d <- did_panel(5)
  expect_error(did(y ~ x, d, "id", "first"), "c\\(individual, time\\)")
  expect_error(did(y ~ x, d, c("id", "year"), "nope"), "not found")
  d2 <- d; d2$first[1] <- 1999
  expect_error(did(y ~ x, d2, c("id", "year"), "first"), "constant within")
  d3 <- d; d3$year <- as.character(d3$year)
  expect_error(did(y ~ x, d3, c("id", "year"), "first"), "numeric")
  expect_error(eventstudy(y ~ x, d, c("id", "year"), "first", leads = 2, ref = -5), "window")
})

test_that("did() and eventstudy() summaries show the user's call, not the internal one", {
  d <- did_panel(3)
  d$bin <- as.numeric(d$y > median(d$y))
  fit <- did(y ~ x, data = d, index = c("id", "year"), treat = "first")
  expect_identical(fit$call[[1]], as.name("did"))
  out <- paste(capture.output(summary(fit)), collapse = "\n")
  expect_match(out, "did(formula = y ~ x", fixed = TRUE)
  expect_no_match(out, "formula = f", fixed = TRUE)
  expect_match(out, "Residual df:", fixed = TRUE)
  # conditional logit: no residual df, so the line is not printed
  fitb <- did(bin ~ x, data = d, index = c("id", "year"), treat = "first",
              family = "binomial")
  outb <- paste(capture.output(summary(fitb)), collapse = "\n")
  expect_match(outb, "did(formula = bin ~ x", fixed = TRUE)
  expect_no_match(outb, "Residual df: NA", fixed = TRUE)
  expect_no_match(outb, "NA", fixed = TRUE)
  es <- eventstudy(y ~ x, data = d, index = c("id", "year"), treat = "first",
                   leads = 2, lags = 2)
  expect_identical(es$call[[1]], as.name("eventstudy"))
  expect_identical(summary(es)$call[[1]], as.name("eventstudy"))
  expect_match(paste(capture.output(summary(es$fit)), collapse = "\n"),
               "eventstudy(formula = y ~ x", fixed = TRUE)
})
