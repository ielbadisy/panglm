# Difference and system GMM: pangmm() against plm::pgmm() and its
# vcovHC(), mtest(), sargan() and pwaldtest() methods.
# Data: plm::EmplUK (140 UK firms, 1976-1984, unbalanced), with the
# employment equations of Arellano and Bond (1991) and Blundell and Bond
# (1998); and a simulated BMI panel with habit persistence (lagged BMI),
# a predetermined smoking indicator, exercise and gaps (dropped waves).
#
# Conventions aligned before comparing:
# * plm's one-step non-robust covariance B1 = (X'Z A1 Z'X)^{-1} with
#   A1 = N (sum_i Z_i' H Z_i)^{-1} lacks the error variance: pangmm's
#   classical one-step covariance is s2 N B1, with
#   s2 = sum(de^2) / (2 (n_d - k)) from the differenced one-step residuals.
#   s2 is recomputed below from plm's own residuals.
# * plm's sargan(weights = "onestep") is z' A1 z = N z'(sum Z'HZ)^{-1} z;
#   pangmm's Sargan statistic divides by s2 instead of multiplying by N.
# * pangmm's Hansen J is plm's sargan() (two-step weights) of the
#   two-step fit, whatever the model of the pangmm fit.
# * The AR tests and Wald test use the robust covariance (vcovHC), as in
#   summary(pgmm()).
# * System GMM: plm instruments the level equation with the difference
#   lagged once whatever the first GMM lag; pangmm uses the difference
#   lagged (a - 1) times, which coincides for the lag(., 2:99) instruments
#   used here. With effect = "individual" pangmm keeps a constant in the
#   level equation and plm does not, so only effect = "twoways" is
#   compared for system GMM.
source("validation/common.R")
suppressPackageStartupMessages(library(plm))
data("EmplUK", package = "plm")
se <- function(v) sqrt(diag(v))
out <- list()
add <- function(...) out[[length(out) + 1L]] <<- vcompare(...)

# simulated BMI panel: 400 people over 8 waves after a 50-wave burn-in
# (mean-stationary start, the system-GMM condition), 10% of person-waves
# missing at random (internal gaps). Waves are numbered 1..8:
# plm::pgmm() lags by the numeric value of the time index, so waves coded
# as years two apart (1996, 1998, ...) would have no lag at all in plm,
# while pangmm lags on the sorted distinct periods.
set.seed(20261025)
N <- 400; TT <- 8; burn <- 50
bmi <- data.frame(id = rep(seq_len(N), each = TT), wave = rep(seq_len(TT), N),
                  exercise = NA_real_, smoke = NA_real_, bmi = NA_real_)
mu <- rnorm(N, 0, 1.5)
for (i in seq_len(N)) {
  r <- (i - 1) * TT + seq_len(TT)
  y <- 22 + mu[i] / 0.4
  e_prev <- 0
  for (t in seq_len(burn + TT)) {
    s <- as.numeric(0.3 * mu[i] + 0.8 * e_prev + rnorm(1) > 0.5)  # predetermined
    x <- rnorm(1, 3, 1)  # weekly hours of exercise
    e <- rnorm(1)
    y <- 10 + 0.6 * y - 0.5 * s - 0.3 * x + mu[i] + e
    if (t > burn) bmi[r[t - burn], c("exercise", "smoke", "bmi")] <- c(x, s, y)
    e_prev <- e
  }
}
bmi <- bmi[-sample(nrow(bmi), round(0.1 * nrow(bmi))), ]

fits <- list(
  list(data = "EmplUK", label = "difference GMM",
       f = log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) + log(capital) +
         lag(log(output), 0:1) | lag(log(emp), 2:99),
       d = EmplUK, idx = c("firm", "year"), args = list()),
  list(data = "EmplUK", label = "difference GMM, collapsed",
       f = log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) + log(capital) +
         lag(log(output), 0:1) | lag(log(emp), 2:99),
       d = EmplUK, idx = c("firm", "year"), args = list(collapse = TRUE)),
  list(data = "EmplUK", label = "difference GMM, individual effect",
       f = log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) + log(capital) +
         lag(log(output), 0:1) | lag(log(emp), 2:99),
       d = EmplUK, idx = c("firm", "year"), args = list(effect = "individual")),
  list(data = "EmplUK", label = "difference GMM, maxlag = 4",
       f = log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) + log(capital) +
         lag(log(output), 0:1) | lag(log(emp), 2:99),
       fplm = log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) + log(capital) +
         lag(log(output), 0:1) | lag(log(emp), 2:4),
       d = EmplUK, idx = c("firm", "year"), args = list(maxlag = 4)),
  list(data = "EmplUK", label = "difference GMM, standard instruments",
       f = log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) + log(capital) |
         lag(log(emp), 2:99) | lag(log(wage), 1:2) + log(capital),
       d = EmplUK, idx = c("firm", "year"), args = list()),
  list(data = "EmplUK", label = "system GMM",
       f = log(emp) ~ lag(log(emp), 1) + lag(log(wage), 0:1) + lag(log(capital), 0:1) |
         lag(log(emp), 2:99) + lag(log(wage), 2:99) + lag(log(capital), 2:99),
       d = EmplUK, idx = c("firm", "year"), args = list(transformation = "ld")),
  list(data = "EmplUK", label = "system GMM, collapsed",
       f = log(emp) ~ lag(log(emp), 1) + lag(log(wage), 0:1) + lag(log(capital), 0:1) |
         lag(log(emp), 2:99) + lag(log(wage), 2:99) + lag(log(capital), 2:99),
       d = EmplUK, idx = c("firm", "year"),
       args = list(transformation = "ld", collapse = TRUE)),
  list(data = "simulated BMI (gaps)", label = "difference GMM",
       f = bmi ~ lag(bmi, 1) + smoke + exercise | lag(bmi, 2:99) + lag(smoke, 1:99),
       d = bmi, idx = c("id", "wave"), args = list()),
  list(data = "simulated BMI (gaps)", label = "difference GMM, collapsed",
       f = bmi ~ lag(bmi, 1) + smoke + exercise | lag(bmi, 2:99) + lag(smoke, 1:99),
       d = bmi, idx = c("id", "wave"), args = list(collapse = TRUE)),
  list(data = "simulated BMI (gaps)", label = "system GMM",
       f = bmi ~ lag(bmi, 1) + smoke + exercise | lag(bmi, 2:99) + lag(smoke, 2:99),
       d = bmi, idx = c("id", "wave"), args = list(transformation = "ld")),
  list(data = "simulated BMI (gaps)", label = "system GMM, collapsed, maxlag = 3",
       f = bmi ~ lag(bmi, 1) + smoke + exercise | lag(bmi, 2:99) + lag(smoke, 2:99),
       fplm = bmi ~ lag(bmi, 1) + smoke + exercise | lag(bmi, 2:3) + lag(smoke, 2:3),
       d = bmi, idx = c("id", "wave"),
       args = list(transformation = "ld", collapse = TRUE, maxlag = 3))
)

# differenced one-step residuals of a plm fit and the implied s2
plm_sigma2 <- function(b) {
  res <- b$residuals
  if (b$args$transformation == "ld") {
    m <- length(b$args$namest) - 1L
    res <- lapply(res, function(x) x[seq_len(m)])
  }
  e <- unlist(res, use.names = FALSE)
  k <- length(if (is.list(b$coefficients)) b$coefficients[[1L]] else b$coefficients)
  sum(e^2) / (2 * (sum(e != 0) - k))
}

for (s in fits) {
  fplm <- if (is.null(s$fplm)) s$f else s$fplm
  plm_args <- s$args[intersect(names(s$args), c("transformation", "collapse", "effect"))]
  b1 <- do.call(plm::pgmm, c(list(fplm, data = s$d, index = s$idx, model = "onestep"), plm_args))
  b2 <- do.call(plm::pgmm, c(list(fplm, data = s$d, index = s$idx, model = "twosteps"), plm_args))
  a1 <- do.call(pangmm, c(list(s$f, data = s$d, index = s$idx, model = "onestep"), s$args))
  a2 <- do.call(pangmm, c(list(s$f, data = s$d, index = s$idx, model = "twostep"), s$args))
  Np <- length(b1$W)
  vb1 <- plm::vcovHC(b1)
  vb2 <- plm::vcovHC(b2)
  m1 <- paste("pangmm one-step", s$label)
  m2 <- paste("pangmm two-step", s$label)
  ref <- function(x) list(plm = unname(x))

  add(m1, s$data, "coefficients", coef(a1), ref(coef(b1)), 1e-6)
  add(m1, s$data, "robust SE", se(vcov(a1, type = "robust")), ref(se(vb1)), 1e-6)
  s2 <- plm_sigma2(b1)
  add(m1, s$data, "error variance s2", a1$sigma2, ref(s2), 1e-6)
  add(m1, s$data, "classical SE (s2 N B1)", se(vcov(a1, type = "classical")),
      ref(se(vcov(b1) * Np * s2)), 1e-6)
  add(m1, s$data, "AR(1), AR(2) tests", a1$artest$statistic,
      ref(c(plm::mtest(b1, 1, vb1)$statistic, plm::mtest(b1, 2, vb1)$statistic)), 1e-6)
  add(m1, s$data, "Sargan (z'A1z / (N s2))", a1$sargan$statistic,
      ref(plm::sargan(b1, weights = "onestep")$statistic / (Np * s2)), 1e-6)
  add(m1, s$data, "Sargan/Hansen df", c(a1$sargan$df, a1$hansen$df),
      ref(rep(plm::sargan(b1)$parameter, 2)), 1e-12)
  add(m1, s$data, "Wald test (robust)", a1$wald$statistic,
      ref(plm::pwaldtest(b1, vcov = vb1)$statistic), 1e-6)

  add(m2, s$data, "coefficients", coef(a2), ref(coef(b2)), 1e-6)
  add(m2, s$data, "one-step coefficients", a2$coefficients_onestep, ref(coef(b1)), 1e-6)
  add(m2, s$data, "Windmeijer SE", se(vcov(a2, type = "robust")), ref(se(vb2)), 1e-6)
  add(m2, s$data, "conventional two-step SE", se(vcov(a2, type = "classical")),
      ref(se(vcov(b2))), 1e-6)
  add(m2, s$data, "AR(1), AR(2) tests", a2$artest$statistic,
      ref(c(plm::mtest(b2, 1, vb2)$statistic, plm::mtest(b2, 2, vb2)$statistic)), 1e-6)
  add(m2, s$data, "Hansen J", a2$hansen$statistic, ref(plm::sargan(b2)$statistic), 1e-6)
  add(m1, s$data, "Hansen J (two-step)", a1$hansen$statistic, ref(plm::sargan(b2)$statistic),
      1e-6)
  add(m2, s$data, "Wald test (robust)", a2$wald$statistic,
      ref(plm::pwaldtest(b2, vcov = vb2)$statistic), 1e-6)
  add(m2, s$data, "number of instruments", a2$n_instruments, ref(ncol(b2$W[[1L]])), 1e-12)
}
do.call(vsave, c(out, name = "gmm"))
