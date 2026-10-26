# Panel quantile regression: panqr() with model = "canay", "mmqr" and
# "pooling" against hand implementations built from plm::plm(), lm() and
# quantreg::rq().
# * Canay (2011): the individual effects from the within regression are
#   subtracted by hand and quantreg::rq() is run on the adjusted outcome.
# * MM-QR (Machado and Santos Silva 2019): location and scale within
#   regressions by plm (and by lm() with individual dummies on a small
#   panel), q(tau) by rq() on a constant, and the asymptotic covariance
#   transcribed from the Stata command xtqreg 1.5, both as published
#   (centered = FALSE) and with the centered scale influence of the paper.
# * Pooled: rq() itself, and its kernel sandwich clustered by individual
#   built from the Hinv matrix of quantreg::summary.rq(se = "ker").
# * Simulated location-scale panels with known beta(tau): recovery of the
#   true quantile coefficients, and Monte Carlo checks that the analytic
#   standard errors match the sampling standard deviation (ratio near 1).
# Data: RAND Health Insurance Experiment panel pglm::HealthIns (5,908
# individuals, 1 to 5 years, unbalanced), outcome log(1 + medical
# expenses).
source("validation/common.R")
suppressPackageStartupMessages({
  library(plm)
  library(quantreg)
})
data("HealthIns", package = "pglm")
hi <- HealthIns
hi$lmed <- log1p(hi$med)
hi$kid <- as.numeric(hi$child == "yes")
hi <- hi[order(hi$id, hi$year), ]
# individuals observed once carry no within information; panqr() drops them
# for the fixed-effects estimators, the hand versions use the same sample
hi_all <- hi
hi <- hi[stats::ave(hi$id, hi$id, FUN = length) > 1, ]
taus <- c(0.25, 0.5, 0.75, 0.9)
f <- lmed ~ age + size + kid
idx <- c("id", "year")
xv <- c("age", "size", "kid")
out <- list()
add <- function(...) out[[length(out) + 1L]] <<- vcompare(...)
tl <- function(tau) paste0("tau=", tau)

# ---- Canay two-step -------------------------------------------------------
a <- panqr(f, hi, index = idx, tau = taus, model = "canay", vcov = "none")
w <- plm::plm(f, hi, index = idx, model = "within")
add("panqr canay", "HealthIns", "within slopes", a$within, list(plm = coef(w)), 1e-8)
a_all <- suppressMessages(panqr(f, hi_all, index = idx, tau = taus, model = "canay",
                                vcov = "none"))
add("panqr canay", "HealthIns", "coefficients, singletons dropped by panqr",
    as.vector(a_all$coefficients), list(panglm = as.vector(a$coefficients)), 1e-12)
# individual effects normalized as in Stata's xtreg, fe (mean zero over obs)
lev <- plm::fixef(w, type = "level")
alpha <- as.numeric(lev[as.character(hi$id)])
alpha <- alpha - mean(alpha)
hi$yadj <- hi$lmed - alpha
for (t in taus) {
  r <- quantreg::rq(yadj ~ age + size + kid, tau = t, data = hi)
  add("panqr canay", "HealthIns", paste("coefficients", tl(t)),
      a$coefficients[, tl(t)], list(quantreg = coef(r)), 1e-8)
}

# ---- Pooled quantile regression ------------------------------------------
p <- panqr(f, hi, index = idx, tau = taus, model = "pooling")
for (j in seq_along(taus)) {
  t <- taus[j]
  r <- quantreg::rq(f, tau = t, data = hi)
  add("panqr pooling", "HealthIns", paste("coefficients", tl(t)),
      p$coefficients[, j], list(quantreg = coef(r)), 1e-10)
  s <- summary(r, se = "ker", covariance = TRUE)
  X <- model.matrix(f, hi)
  e <- residuals(r)
  meat <- crossprod(rowsum(X * (t - (e < 0)), hi$id))
  vh <- s$Hinv %*% meat %*% s$Hinv
  add("panqr pooling", "HealthIns", paste("cluster kernel SE", tl(t)),
      sqrt(diag(p$vcov[[j]])), list(quantreg = sqrt(diag(vh))), 1e-8)
}

# ---- MM-QR ---------------------------------------------------------------
# hand implementation: within regressions (plm or lm with dummies), rq()
# for q(tau) and for the Hall-Sheather difference quotient of the density,
# and the xtqreg 1.5 covariance written out step by step
mmqr_hand <- function(y, X, id, taus, loc_fit) {
  pieces <- loc_fit(y, X, id)
  u <- pieces$resid
  up <- (u >= 0) - mean(u >= 0)
  au <- 2 * u * up
  sp <- loc_fit(au, X, id)
  sigma <- au - sp$resid
  us <- u / sigma
  # within-demeaned regressors and the scale at the overall regressor means
  Xd <- X - apply(X, 2, function(v) ave(v, id))
  xbar_i <- apply(X, 2, function(v) ave(v, id))
  si <- ave(au, id) - as.vector(sweep(xbar_i, 2, colMeans(X)) %*% sp$coef)
  n <- length(y); k <- ncol(X)
  Qxx <- crossprod(Xd)
  Pxx <- t(Xd) %*% (sigma^2 * Xd)
  Px <- colSums(sigma^2 * Xd)
  res <- list(location = pieces$coef, scale = sp$coef)
  for (t in taus) {
    rq1 <- function(p) unname(coef(suppressWarnings(quantreg::rq(us ~ 1, tau = p))))
    q <- rq1(t)
    h <- quantreg::bandwidth.rq(t, n, hs = TRUE)
    fq <- 2 * h / (rq1(t + h) - rq1(t - h))
    V <- list()
    for (centered in c(FALSE, TRUE)) {
      vp <- (us >= 0) - mean(us >= 0)
      av <- 2 * us * vp - centered
      wq <- (t - (us <= q)) / fq - us - q * av
      O <- rbind(cbind(mean(us^2) * Pxx, mean(us * av) * Pxx, mean(us * wq) * Px),
                 cbind(mean(us * av) * Pxx, mean(av^2) * Pxx, mean(av * wq) * Px),
                 c(mean(us * wq) * Px, mean(av * wq) * Px, sum(wq^2) * mean(si^2)))
      XI <- cbind(solve(Qxx), solve(Qxx) * q, sp$coef / sum(si))
      V[[if (centered) "paper" else "xtqreg"]] <- XI %*% O %*% t(XI)
    }
    res[[tl(t)]] <- list(q = q, coef = pieces$coef + sp$coef * q, V = V)
  }
  res
}
within_plm <- function(y, X, id) {
  d <- data.frame(y = y, X, id = id, tt = stats::ave(id, id, FUN = seq_along))
  fit <- plm::plm(stats::reformulate(colnames(X), "y"), d, index = c("id", "tt"),
                  model = "within")
  list(coef = coef(fit), resid = as.numeric(residuals(fit)))
}
within_lm <- function(y, X, id) {
  fit <- lm(y ~ X + factor(id))
  b <- coef(fit)[paste0("X", colnames(X))]
  names(b) <- colnames(X)
  # within residual = OLS residual of the dummy-variable regression
  list(coef = b, resid = as.numeric(residuals(fit)))
}

check_mmqr <- function(dat, dname, f, xv, taus, loc_fit, ref) {
  m <- panqr(f, dat, index = idx, tau = taus, model = "mmqr", vcov = "iid")
  X <- as.matrix(dat[, xv])
  hand <- mmqr_hand(dat[[all.vars(f)[1]]], X, dat$id, taus, loc_fit)
  rl <- function(v) stats::setNames(list(v), ref)
  add("panqr mmqr", dname, "location coefficients", m$location$coefficients,
      rl(hand$location), 1e-8)
  add("panqr mmqr", dname, "scale coefficients", m$scale$coefficients, rl(hand$scale), 1e-8)
  for (j in seq_along(taus)) {
    hj <- hand[[tl(taus[j])]]
    add("panqr mmqr", dname, paste("q(tau)", tl(taus[j])), m$q[j], rl(hj$q), 1e-10)
    add("panqr mmqr", dname, paste("coefficients", tl(taus[j])), m$coefficients[, j],
        rl(hj$coef), 1e-8)
    add("panqr mmqr", dname, paste("analytic SE (paper)", tl(taus[j])),
        sqrt(diag(m$vcov[[j]])), rl(sqrt(diag(hj$V$paper))), 1e-8)
    # xtqreg 1.5 as published, through the internal switch
    Xs <- X
    loc <- panglm:::qr_within_fit(dat[[all.vars(f)[1]]], Xs, match(dat$id, unique(dat$id)))
    g <- m$scale$coefficients
    grp <- match(dat$id, unique(dat$id))
    u <- loc$residuals
    au <- 2 * u * ((u >= 0) - mean(u >= 0))
    au_i <- ave(au, grp)
    sigma <- au_i + as.vector(loc$Xw %*% g)
    si <- au_i - as.vector(sweep(loc$xbar_i, 2, colMeans(Xs))[grp, , drop = FALSE] %*% g)
    vx <- panglm:::mmqr_vcov(loc$Xw, u / sigma, sigma, si, g, m$q[j], taus[j], centered = FALSE)
    add("panqr mmqr", dname, paste("analytic SE (xtqreg 1.5)", tl(taus[j])),
        sqrt(diag(vx)), rl(sqrt(diag(hj$V$xtqreg))), 1e-8)
  }
  m
}

m <- suppressWarnings(check_mmqr(hi, "HealthIns", f, xv, taus, within_plm, "plm+quantreg"))
# location and scale covariances: cluster-robust within, as plm::vcovHC(sss)
wl <- plm::plm(f, hi, index = idx, model = "within")
add("panqr mmqr", "HealthIns", "location cluster SE", sqrt(diag(m$location$vcov)),
    list(plm = sqrt(diag(plm::vcovHC(wl, method = "arellano", type = "sss", cluster = "group")))),
    1e-8)
hi$au <- {
  u <- as.numeric(residuals(wl))
  2 * u * ((u >= 0) - mean(u >= 0))
}
ws <- plm::plm(au ~ age + size + kid, hi, index = idx, model = "within")
add("panqr mmqr", "HealthIns", "scale cluster SE", sqrt(diag(m$scale$vcov)),
    list(plm = sqrt(diag(plm::vcovHC(ws, method = "arellano", type = "sss", cluster = "group")))),
    1e-8)

# small simulated panel, hand version from lm() with individual dummies
set.seed(20261025)
N <- 61; TT <- 7
sm <- data.frame(id = rep(seq_len(N), each = TT), year = rep(seq_len(TT), N))
al <- rep(rnorm(N), each = TT)
sm$x1 <- rnorm(N * TT) + al
sm$x2 <- runif(N * TT)
sm$y <- al + sm$x1 + 0.5 * sm$x2 + (1 + 0.5 * sm$x2) * rnorm(N * TT)
invisible(check_mmqr(sm, "simulated N=61 T=7", y ~ x1 + x2, c("x1", "x2"), c(0.1, 0.5, 0.9),
           within_lm, "lm+quantreg"))

# ---- Recovery of known beta(tau) in a location-scale panel ---------------
# y = alpha_i + x + 0.5 z + (1 + 0.5 z) e, e ~ N(0, 1), x correlated with
# alpha_i: beta_x(tau) = 1 and beta_z(tau) = 0.5 + 0.5 qnorm(tau)
sim_ls <- function(N, TT) {
  d <- data.frame(id = rep(seq_len(N), each = TT), year = rep(seq_len(TT), N))
  a <- rep(rnorm(N), each = TT)
  d$x <- rnorm(N * TT) + 0.5 * a
  d$z <- runif(N * TT, 0, 2)
  d$y <- a + d$x + 0.5 * d$z + (1 + 0.5 * d$z) * rnorm(N * TT)
  d
}
truth <- function(t) c(x = 1, z = 0.5 + 0.5 * qnorm(t))
set.seed(1)
big <- sim_ls(500, 60)
rt <- c(0.1, 0.25, 0.5, 0.75, 0.9)
cb <- panqr(y ~ x + z, big, index = idx, tau = rt, model = "canay", vcov = "none")
mb <- panqr(y ~ x + z, big, index = idx, tau = rt, model = "mmqr", vcov = "none")
pb <- panqr(y ~ x + z, big, index = idx, tau = rt, model = "pooling", vcov = "none")
for (j in seq_along(rt)) {
  add("panqr canay", "location-scale N=500 T=60", paste("recovery", tl(rt[j])),
      cb$coefficients[c("x", "z"), j], list(truth = truth(rt[j])), 0.05)
  add("panqr mmqr", "location-scale N=500 T=60", paste("recovery", tl(rt[j])),
      mb$coefficients[c("x", "z"), j], list(truth = truth(rt[j])), 0.05)
}
# the pooled fit ignores alpha_i, which is correlated with x: its slope on x
# is biased (about 1.4), which the fixed-effects estimators remove
stopifnot(abs(pb$coefficients["x", "tau=0.5"] - 1) > 0.2)

# ---- Monte Carlo: analytic SE versus sampling SD -------------------------
mc <- lapply(seq_len(300), function(s) {
  set.seed(1000 + s)
  d <- sim_ls(200, 30)
  cn <- panqr(y ~ x + z, d, index = idx, tau = c(0.5, 0.9), model = "canay", vcov = "analytic")
  mm <- panqr(y ~ x + z, d, index = idx, tau = c(0.5, 0.9), model = "mmqr")
  c(cn = as.vector(cn$coefficients[c("x", "z"), ]),
    cnse = as.vector(sapply(cn$vcov, function(v) sqrt(diag(v))[c("x", "z")])),
    mm = as.vector(mm$coefficients),
    mmse = as.vector(sapply(mm$vcov, function(v) sqrt(diag(v)))))
})
mc <- do.call(rbind, mc)
lab <- c("x tau=0.5", "z tau=0.5", "x tau=0.9", "z tau=0.9")
ratio <- function(est, se) {
  stats::setNames(colMeans(mc[, se]) / apply(mc[, est], 2, stats::sd), lab)
}
add("panqr canay", "Monte Carlo N=200 T=30 (300 reps)", "analytic SE / MC SD",
    ratio(paste0("cn", 1:4), paste0("cnse", 1:4)), list(one = rep(1, 4)), 0.15)
add("panqr mmqr", "Monte Carlo N=200 T=30 (300 reps)", "analytic SE / MC SD",
    ratio(paste0("mm", 1:4), paste0("mmse", 1:4)), list(one = rep(1, 4)), 0.15)

# cluster bootstrap against the analytic Canay SE on one data set
set.seed(7)
d1 <- sim_ls(200, 30)
cbt <- panqr(y ~ x + z, d1, index = idx, tau = 0.9, model = "canay", vcov = "bootstrap", boot = 300)
can <- panqr(y ~ x + z, d1, index = idx, tau = 0.9, model = "canay", vcov = "analytic")
# slopes only: the bootstrap intercept also reflects the sampling variation of the
# mean individual effect, which the analytic (fixed-effects) SE conditions on
sl <- c("x", "z")
add("panqr canay", "simulated N=200 T=30", "bootstrap SE / analytic SE (tau=0.9)",
    sqrt(diag(vcov(cbt))[sl]) / sqrt(diag(vcov(can))[sl]), list(one = rep(1, 2)), 0.2)

do.call(vsave, c(out, name = "qr"))
