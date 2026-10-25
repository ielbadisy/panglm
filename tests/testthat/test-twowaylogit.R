twowaylogit_data <- function(G = 150, T = 5, seed = 7) {
  set.seed(seed)
  d <- expand.grid(wave = seq_len(T), id = seq_len(G))[, c("id", "wave")]
  a <- rnorm(G)
  d$price <- rnorm(nrow(d)) + 0.2 * d$wave
  d$income <- rnorm(nrow(d), 0, 0.5) + 0.3 * a[d$id]
  shock <- c(0, 0.3, -0.4, -0.7, -0.2, 0.1)[seq_len(T)]
  d$smoking <- rbinom(nrow(d), 1, plogis(a[d$id] - 0.5 * d$price + 0.4 * d$income +
                                             shock[d$wave]))
  d <- d[!(d$wave > 3 & d$id %% 4 == 0) & !(d$wave == 1 & d$id %% 7 == 0), ]
  rownames(d) <- NULL
  d
}

fit_twowaylogit <- function(d, ...) {
  suppressMessages(panglm(smoking ~ price + income, d, c("id", "wave"),
                          model = "within", family = "binomial", effect = "twoways", ...))
}

test_that("two-way conditional logit matches survival::clogit with time dummies", {
  skip_if_missing("survival")
  suppressMessages(library(survival))
  d <- twowaylogit_data()
  f <- fit_twowaylogit(d)
  cl <- suppressWarnings(clogit(
    smoking ~ price + income + factor(wave) + strata(id), data = d,
    method = "exact", control = coxph.control(eps = 1e-12, toler.chol = 1e-14)
  ))
  expect_equal(unname(coef(f)), unname(coef(cl)[1:2]), tolerance = 1e-7)
  expect_equal(unname(f$time_effects[-1]), unname(coef(cl)[-(1:2)]), tolerance = 1e-7)
  expect_equal(unname(f$time_effects[1]), 0)
  expect_equal(names(f$time_effects), as.character(1:5))
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(vcov(cl)))[1:2]),
               tolerance = 1e-7)
  expect_equal(unname(f$time_effects_se[-1]), unname(sqrt(diag(vcov(cl)))[-(1:2)]),
               tolerance = 1e-7)
  expect_equal(as.numeric(logLik(f)), cl$loglik[2], tolerance = 1e-8)
  expect_equal(attr(logLik(f), "df"), length(coef(cl)))
  expect_named(coef(f), c("price", "income"))
})

test_that("two-way conditional logit robust covariances match a T = 2 differenced logit", {
  skip_if_missing("sandwich")
  d <- twowaylogit_data(G = 300, T = 2)
  d <- d[d$id %in% names(which(table(d$id) == 2)), ]
  f <- fit_twowaylogit(d)
  w1 <- d[d$wave == 1, ]; w2 <- d[d$wave == 2, ]
  fd <- data.frame(y = w2$smoking, dprice = w2$price - w1$price,
                   dincome = w2$income - w1$income)[w1$smoking != w2$smoking, ]
  g <- glm(y ~ dprice + dincome, binomial, fd, control = glm.control(epsilon = 1e-14))
  b <- c("dprice", "dincome")
  expect_equal(unname(coef(f)), unname(coef(g)[b]), tolerance = 1e-7)
  expect_equal(unname(f$time_effects[2]), unname(coef(g)[1]), tolerance = 1e-7)
  expect_equal(unname(sqrt(diag(vcov(f, type = "HC1")))),
               unname(sqrt(diag(sandwich::vcovHC(g, type = "HC1")))[b]), tolerance = 1e-6)
  expect_equal(unname(sqrt(diag(vcov(f, type = "cluster")))),
               unname(sqrt(diag(sandwich::vcovHC(g, type = "HC1")))[b]), tolerance = 1e-6)
})

test_that("two-way conditional logit vcov argument, summary, print and predict", {
  d <- twowaylogit_data()
  f <- fit_twowaylogit(d, vcov = "cluster")
  expect_identical(f$vcov_type, "cluster")
  expect_equal(vcov(f), vcov(f, type = "cluster"))
  expect_false(isTRUE(all.equal(vcov(f), vcov(f, type = "classical"))))
  expect_equal(dim(vcov(f, type = "HC1")), c(2L, 2L))
  coarse <- (d$id - 1) %/% 5
  expect_equal(dim(vcov(f, type = "cluster", cluster = coarse)), c(2L, 2L))
  expect_error(vcov(f, type = "cluster", cluster = seq_len(nrow(d))), "must not split")

  s <- summary(f)
  expect_identical(s$effect, "twoways")
  expect_equal(rownames(s$coefficients), c("price", "income"))
  expect_output(print(s), "Effect: twoways")

  fitted_p <- fitted(f)
  expect_length(fitted_p, nrow(d))
  expect_true(all(fitted_p >= 0 & fitted_p <= 1))
  # conditional probabilities add up to each individual's success total
  expect_equal(unname(tapply(fitted_p, d$id, sum)), unname(tapply(d$smoking, d$id, sum)),
               tolerance = 1e-8)
  expect_equal(predict(f, type = "response"), fitted_p)
  expect_error(predict(f, newdata = d[1:3, ]), "not defined for exact conditional")
})

test_that("two-way conditional logit is invariant to row order and drops a time trend", {
  d <- twowaylogit_data()
  f <- fit_twowaylogit(d)
  f2 <- fit_twowaylogit(d[sample(nrow(d)), ])
  expect_equal(coef(f), coef(f2), tolerance = 1e-8)
  expect_equal(f$time_effects, f2$time_effects, tolerance = 1e-8)

  # a slope collinear with the time effects (a linear trend) is dropped,
  # not traded against a time dummy
  d$trend <- d$wave
  expect_warning(
    f3 <- suppressMessages(panglm(smoking ~ price + trend, d, c("id", "wave"),
                                  model = "within", family = "binomial",
                                  effect = "twoways")),
    "trend"
  )
  expect_named(coef(f3), "price")
  expect_true(all(is.finite(f3$time_effects)))
})

test_that("two-way conditional logit drops a wave seen only in uninformative strata", {
  d <- twowaylogit_data()
  # wave 6 observed only for two individuals who never smoke
  extra <- d[d$id %in% c(1, 2) & d$wave == 1, ]
  extra$wave <- 6
  never <- d$id %in% c(1, 2)
  d$smoking[never] <- 0
  extra$smoking <- 0
  d <- rbind(d, extra)
  f <- fit_twowaylogit(d)
  expect_true(is.na(f$time_effects[["6"]]))
  expect_true(all(is.finite(f$time_effects[as.character(1:5)])))
  expect_equal(attr(logLik(f), "df"), 2L + 4L)
  expect_true(all(is.finite(sqrt(diag(vcov(f, type = "cluster"))))))
})

test_that("effect = 'twoways' is still rejected for unsupported binomial models", {
  d <- twowaylogit_data()
  expect_error(panglm(smoking ~ price, d, c("id", "wave"), model = "pooling",
                      family = "binomial", effect = "twoways"), "only implemented")
  expect_error(panglm(smoking ~ price, d, c("id", "wave"), model = "within",
                      family = "probit", effect = "twoways"), "conditional logit")
})
