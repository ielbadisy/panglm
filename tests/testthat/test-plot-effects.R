effects_data <- function() {
  data(copd, package = "panglm", envir = environment())
  copd$smk <- factor(copd$smoker, labels = c("no", "yes"))
  copd
}

with_pdf <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

test_that("response, slope, factor and by curves match marginaleffects on pooled fits", {
  skip_if_missing("marginaleffects")
  skip_if_missing("MASS")
  old <- options(marginaleffects_numDeriv = list(method = "Richardson"))
  on.exit(options(old), add = TRUE)
  d <- effects_data()
  vals <- c(1, 3, 5, 7)
  for (fam in c("binomial", "probit", "poisson", "gamma", "negbin")) {
    y <- switch(fam, binomial = , probit = "hospitalized", poisson = , negbin = "exacerbations",
                gamma = "fev1")
    form <- stats::as.formula(paste(y, "~ crp + I(crp^2) + age + smk"))
    f <- panglm(form, d, c("id", "visit"), family = fam)
    g <- switch(fam,
      binomial = glm(form, binomial, d),
      # glm's default convergence leaves probit coefficients about 1e-5 off
      probit = glm(form, binomial("probit"), d,
                   control = glm.control(epsilon = 1e-14, maxit = 100)),
      poisson = glm(form, poisson, d), gamma = glm(form, Gamma("log"), d),
      negbin = MASS::glm.nb(form, d)
    )
    grid <- function(...) marginaleffects::datagrid(..., model = g, grid_type = "counterfactual")
    # panglm's probit optimiser stops about 1e-6 (relative) from glm's
    tol <- if (fam == "probit") 1e-5 else 1e-6
    check <- function(a, m) {
      m <- as.data.frame(m)
      expect_equal(a$estimate, m$estimate, tolerance = tol)
      expect_equal(a$se, m$std.error, tolerance = 1e-4)
    }
    check(effectcurve(f, "crp", values = vals, scale = "response"),
          marginaleffects::avg_predictions(g, newdata = grid(crp = vals), by = "crp"))
    check(effectcurve(f, "crp", values = vals, scale = "slope"),
          marginaleffects::avg_slopes(g, variables = "crp", newdata = grid(crp = vals),
                                      by = "crp", numderiv = "fdcenter"))
    check(effectcurve(f, "smk", scale = "response"),
          marginaleffects::avg_predictions(g, variables = "smk", by = "smk"))
    check(effectcurve(f, "crp", values = vals, scale = "response", by = "smk"),
          marginaleffects::avg_predictions(g, newdata = grid(crp = vals, smk = c("no", "yes")),
                                           by = c("smk", "crp")))
    a <- partialeffects(f)
    a$se <- a$std.error
    m <- as.data.frame(marginaleffects::avg_slopes(g, numderiv = "fdcenter"))
    check(a, m[match(a$term, m$term), ])
  }
})

test_that("spline response curves equal averaged glm counterfactual predictions", {
  d <- effects_data()
  form <- exacerbations ~ splines::ns(crp, 3) + age + smk
  f <- panglm(form, d, c("id", "visit"), family = "poisson")
  g <- glm(form, poisson, d)
  vals <- c(0.5, 2, 4, 6)
  manual <- vapply(vals, function(v) {
    nd <- d
    nd$crp <- v
    mean(predict(g, nd, type = "response"))
  }, numeric(1))
  ec <- effectcurve(f, "crp", values = vals, scale = "response")
  expect_equal(ec$estimate, manual, tolerance = 1e-8)
  expect_equal(attr(ec, "scale"), "response")
  expect_true(length(attr(ec, "rug")) <= 101L)
  # the slope curve is the derivative of the response curve
  sl <- effectcurve(f, "crp", values = 3, scale = "slope")
  rs <- effectcurve(f, "crp", values = 3 + c(-1, 1) * 1e-4, scale = "response")
  expect_equal(sl$estimate, diff(rs$estimate) / 2e-4, tolerance = 1e-5)
})

test_that("random-effects and CRE curves use the partialeffects() mean", {
  d <- effects_data()
  pr <- panglm(hospitalized ~ crp + age + smk, d, c("id", "visit"), model = "random",
               family = "probit", cre = TRUE, vcov = "cluster")
  po <- panglm(exacerbations ~ splines::ns(crp, 3) + smk, d, c("id", "visit"),
               model = "random", family = "poisson")
  for (f in list(pr, po)) {
    ec <- effectcurve(f, "smk", scale = "response")
    ape <- partialeffects(f, variables = "smk")
    expect_equal(diff(ec$estimate), ape$estimate, tolerance = 1e-10)
    expect_true(all(is.finite(ec$se) & ec$se > 0))
  }
  lin <- panglm(fev1 ~ crp + age, d, c("id", "visit"), model = "random")
  sl <- effectcurve(lin, "crp", values = c(1, 4), scale = "slope")
  expect_equal(sl$estimate, rep(unname(coef(lin)["crp"]), 2), tolerance = 1e-8)
  expect_equal(sl$se, rep(unname(sqrt(vcov(lin)["crp", "crp"])), 2), tolerance = 1e-6)
})

test_that("within fits refuse the response scale but give slopes", {
  d <- effects_data()
  fe <- panglm(exacerbations ~ splines::ns(crp, 3), d, c("id", "visit"),
               model = "within", family = "poisson")
  expect_error(effectcurve(fe, "crp", scale = "response"), "scale = \"link\"")
  sl <- effectcurve(fe, "crp", scale = "slope")
  expect_true(all(is.finite(sl$estimate)))
  cl <- panglm(hospitalized ~ crp, d, c("id", "visit"), model = "within", family = "binomial")
  expect_error(effectcurve(cl, "crp", scale = "slope"), "not identified")
})

test_that("effectcurve() argument errors are informative", {
  d <- effects_data()
  f <- panglm(exacerbations ~ crp + age + smk, d, c("id", "visit"), family = "poisson")
  expect_error(effectcurve(f, "smk", scale = "slope"), "continuous")
  expect_error(effectcurve(f, "crp", scale = "response", by = "bmi"), "regressor")
  expect_error(effectcurve(f, "crp", by = "smk"), "by")
  expect_error(effectcurve(f, "crp", scale = "response", by = "age"), "at most 10")
  expect_error(effectcurve(f, "crp", scale = "response", by = "crp"), "differ")
  expect_error(effectcurve(list(), "crp"), "panglm fit")
})

test_that("every effect plot path draws and returns invisibly", {
  d <- make_cost_data()
  d$smk <- factor(d$smoker, labels = c("no", "yes"))
  f <- panglm(exacerbations ~ splines::ns(crp, 3) + age + smk + treatment, d, c("id", "visit"),
              model = "random", family = "poisson", cre = TRUE)
  fe <- panglm(fev1 ~ splines::ns(crp, 3), d, c("id", "visit"), model = "within")
  tp <- twopart(cost0 ~ treatment + crp, d, c("id", "visit"), positive = "gamma")
  curves <- list(
    effectcurve(fe, "crp"),
    effectcurve(f, "crp", values = seq(1, 6, length.out = 10), scale = "response"),
    effectcurve(f, "crp", values = seq(1, 6, length.out = 10), scale = "response", by = "smk"),
    effectcurve(f, "crp", values = seq(1, 6, length.out = 10), scale = "slope"),
    effectcurve(f, "crp", values = seq(1, 6, length.out = 10), scale = "slope", by = "smk"),
    effectcurve(f, "smk", scale = "response"),
    effectcurve(f, "smk", scale = "response", by = "treatment")
  )
  with_pdf({
    for (ec in curves) {
      expect_no_error(expect_invisible(plot(ec)))
    }
    expect_invisible(plot(curves[[3]], legend = FALSE, col = c("red", "blue"), rug = FALSE))
    expect_invisible(plot(curves[[3]], legend = "top"))
    a <- partialeffects(f)
    expect_no_error(expect_invisible(plot(a)))
    expect_invisible(plot(a, labels = letters[seq_len(nrow(a))], main = "APE"))
    expect_error(plot(a, labels = "x"), "one entry per row")
    for (part in c("total", "binary", "positive")) {
      expect_no_error(expect_invisible(plot(partialeffects(tp, part = part))))
    }
    mai <- graphics::par("mai")
    plot(a)
    expect_equal(graphics::par("mai"), mai)
  })
})
