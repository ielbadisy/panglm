skip_if_missing("sandwich")
skip_if_missing("MASS")
data("birthwt", package = "MASS", envir = environment())
bw <- transform(birthwt, id = seq_along(bwt), race = factor(race))
hc3 <- function(fit) sandwich::vcovHC(fit, type = "HC3")
rel <- function(a, b) max(abs(a - b)) / max(abs(b))

test_that("HC3 equals sandwich::vcovHC for pooled lm and glm families", {
  f <- bwt ~ smoke + age + lwt + race
  m <- panglm(f, bw, index = "id", model = "pooling", vcov = "HC3")
  expect_lt(rel(vcov(m), hc3(lm(f, bw))), 1e-8)
  g <- panglm(low ~ smoke + age + lwt + race, bw, index = "id", model = "pooling", family = "binomial")
  expect_lt(rel(vcov(g, type = "HC3"), hc3(glm(low ~ smoke + age + lwt + race, binomial, bw))), 1e-5)
  p <- panglm(ptl ~ smoke + age, bw, index = "id", model = "pooling", family = "poisson")
  expect_lt(rel(vcov(p, type = "HC3"), hc3(glm(ptl ~ smoke + age, poisson, bw))), 1e-5)
  ga <- panglm(bwt ~ smoke + age, bw, index = "id", model = "pooling", family = gammafamily("log"))
  expect_lt(rel(vcov(ga, type = "HC3"), hc3(glm(bwt ~ smoke + age, Gamma("log"), bw))), 1e-4)
  nb <- panglm(ftv ~ smoke + age, bw, index = "id", model = "pooling", family = "negbin")
  expect_lt(rel(vcov(nb, type = "HC3"), hc3(MASS::glm.nb(ftv ~ smoke + age, bw))), 1e-3)
  # the default variance is unchanged
  expect_equal(vcov(panglm(f, bw, index = "id", model = "pooling")), vcov(lm(f, bw)), tolerance = 1e-8)
})

test_that("HC3 with prior weights equals sandwich::vcovHC", {
  set.seed(2); bw$w <- runif(nrow(bw), 0.5, 2)
  g <- panglm(low ~ smoke + age + lwt, bw, index = "id", model = "pooling", family = "binomial", weights = "w")
  r <- suppressWarnings(glm(low ~ smoke + age + lwt, binomial, bw, weights = w))
  expect_lt(rel(vcov(g, type = "HC3"), hc3(r)), 1e-5)
})

test_that("HC3 for the within model equals the dummy-variable lm on an unbalanced, unsorted panel", {
  set.seed(3)
  d <- data.frame(id = rep(1:40, times = sample(2:6, 40, TRUE)))
  d$t <- ave(d$id, d$id, FUN = seq_along)
  d$x <- rnorm(nrow(d)) + d$id / 20; d$z <- rnorm(nrow(d))
  d$y <- d$id / 10 + d$x - 0.5 * d$z + rnorm(nrow(d), sd = 1 + abs(d$x))
  d <- d[sample(nrow(d)), ]
  w <- panglm(y ~ x + z, d, index = c("id", "t"), model = "within", vcov = "HC3")
  v <- hc3(lm(y ~ x + z + factor(id), d))[c("x", "z"), c("x", "z")]
  expect_lt(rel(vcov(w), v), 1e-8)
})

test_that("HC3 stops for models it does not cover", {
  d <- data.frame(id = rep(1:30, each = 4), t = 1:4, x = rnorm(120))
  d$y <- rnorm(30)[d$id] + d$x + rnorm(120)
  r <- panglm(y ~ x, d, index = c("id", "t"), model = "random")
  expect_error(vcov(r, type = "HC3"), "HC3")
})
