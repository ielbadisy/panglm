# Validation of the two-way fixed-effects conditional (Chamberlain) logit,
# panglm(model = "within", family = "binomial", effect = "twoways").
#
# Reference: survival::clogit(y ~ x + factor(time) + strata(id),
# method = "exact"), the same estimator (individual effects conditioned
# out, time effects as dummies inside the conditional likelihood).
# Compared: slopes, time effects, classical standard errors of both, and
# the conditional log-likelihood.
#
# Robust covariances: clogit cannot give a robust covariance with
# method = "exact" (no dfbeta residuals), so two independent references
# are used.
#   (a) Brute force: each stratum's score is computed in R by enumerating
#       every 0/1 sequence with the observed number of successes (combn),
#       not by panglm's recursion, and combined with clogit's bread. The
#       sandwich is over the full parameter vector (slopes + time dummies)
#       and the slope block is compared.
#   (b) T = 2: the conditional logit is the logit of y_i2 on x_i2 - x_i1
#       among individuals whose outcome changes; the time effect is its
#       intercept. glm + sandwich (HC1, and vcovCL with one cluster per
#       individual) give the reference.
# Small-sample factors (m informative strata, K = slopes + estimated time
# dummies, G clusters): HC1 m / (m - K); cluster G / (G - 1) (m - 1) /
# (m - K). With one cluster per stratum both equal m / (m - K), which is
# sandwich's HC1 on the T = 2 differenced logit (K counts its intercept).
#
# Data: a simulated unbalanced smoking panel (survey waves with year
# shocks), the copd package data (hospitalization by visit), and
# pglm::UnionWage (union membership 1980-1987, made unbalanced).
#
# Run from the package root: Rscript validation/validate-twowaylogit.R
source("validation/common.R")
suppressPackageStartupMessages({
  library(survival)
  library(sandwich)
})
se <- function(v) sqrt(diag(v))
tol <- 1e-6
out <- list()

# ---- data -------------------------------------------------------------------
# smoking (0/1) over survey waves, with a common year shock (e.g. a tobacco
# tax or an anti-smoking campaign), a cigarette price index and log income
# that vary within individuals, and individual propensities correlated with
# income. Unbalanced: later waves missing for some respondents, the first
# wave missing for others.
sim_smoking <- function(G = 400, T = 7, seed = 11) {
  set.seed(seed)
  d <- expand.grid(wave = seq_len(T), id = seq_len(G))[, c("id", "wave")]
  a <- rnorm(G, 0, 1.2)
  d$income <- rnorm(nrow(d), 0, 0.5) + 0.4 * a[d$id]
  d$price <- rnorm(nrow(d), 0, 1) + 0.2 * d$wave
  d$stress <- rbinom(nrow(d), 1, 0.3)
  shock <- c(0, 0.2, -0.3, -0.6, -0.4, -0.9, -1.1)[seq_len(T)]
  d$smoking <- rbinom(nrow(d), 1, plogis(a[d$id] - 0.6 * d$price + 0.3 * d$income +
                                             0.5 * d$stress + shock[d$wave]))
  d <- d[!(d$wave > 5 & d$id %% 4 == 0) & !(d$wave == 1 & d$id %% 9 == 0), ]
  rownames(d) <- NULL
  d
}
smk <- sim_smoking()
sname <- sprintf("simulated smoking, unbalanced (%d ids, %d obs)",
                 length(unique(smk$id)), nrow(smk))

data(copd, package = "panglm")
cname <- "copd (60 ids, 5 visits)"

data("UnionWage", package = "pglm")
uw <- UnionWage
uw$union01 <- as.numeric(uw$union == "yes")
uw$health01 <- as.numeric(uw$health == "yes")
uw$rural01 <- as.numeric(uw$rural == "yes")
uw <- uw[!(uw$year >= 1986 & uw$id %% 3 == 0) & !(uw$year == 1980 & uw$id %% 5 == 0), ]
uname <- sprintf("UnionWage, unbalanced (%d obs)", nrow(uw))

# ---- brute-force stratum scores (independent of panglm's recursion) --------
brute_scores <- function(Z, y, id, theta) {
  eta <- as.numeric(Z %*% theta)
  groups <- split(seq_along(y), id)
  t(vapply(groups, function(rows) {
    S <- sum(y[rows]); Tn <- length(rows)
    if (S == 0 || S == Tn) return(rep(NA_real_, ncol(Z)))
    sets <- utils::combn(Tn, S, simplify = FALSE)
    lw <- vapply(sets, function(s) sum(eta[rows][s]), 0)
    w <- exp(lw - max(lw)); w <- w / sum(w)
    ex <- Reduce(`+`, Map(function(s, wi) wi * colSums(Z[rows[s], , drop = FALSE]), sets, w))
    colSums(Z[rows[y[rows] == 1], , drop = FALSE]) - ex
  }, numeric(ncol(Z))))
}

# robust references: HC1 (one stratum = one unit) and cluster (by id, or
# by a coarser cluster made of whole strata)
brute_robust <- function(cl, Z, y, id, cluster, k) {
  theta <- coef(cl)
  U <- brute_scores(Z, y, id, theta)
  used <- stats::complete.cases(U)
  U <- U[used, , drop = FALSE]
  m <- nrow(U); K <- ncol(U)
  B <- vcov(cl)
  hc1 <- B %*% crossprod(U) %*% B * m / (m - K)
  cl_of_id <- tapply(cluster, id, function(v) as.character(v[1]))[used]
  Uc <- rowsum(U, cl_of_id)
  Gc <- nrow(Uc)
  clu <- B %*% crossprod(Uc) %*% B * (Gc / (Gc - 1)) * ((m - 1) / (m - K))
  list(HC1 = se(hc1)[seq_len(k)], cluster = se(clu)[seq_len(k)])
}

check <- function(label, dname, formula, data, idx, xnames, cluster_coarse = NULL) {
  y <- all.vars(formula)[1]
  p <- suppressMessages(panglm(formula, data, idx, model = "within",
                               family = "binomial", effect = "twoways"))
  cf <- stats::update(formula, stats::as.formula(
    sprintf(". ~ . + factor(%s) + strata(%s)", idx[2], idx[1])))
  cl <- suppressWarnings(clogit(cf, data, method = "exact",
                                control = coxph.control(eps = 1e-12, toler.chol = 1e-14,
                                                        iter.max = 100)))
  k <- length(xnames)
  tcl <- coef(cl)[-seq_len(k)]
  tp <- p$time_effects[-1]
  names(tp) <- names(tcl)

  # design in clogit's column order, rows in data order
  Z <- model.matrix(stats::update(formula, stats::as.formula(
    sprintf(". ~ . + factor(%s)", idx[2]))), data)[, -1, drop = FALSE]
  Z <- Z[, names(coef(cl)), drop = FALSE]
  yv <- data[[y]]; if (!is.numeric(yv)) yv <- as.numeric(yv == levels(yv)[2])
  rob <- brute_robust(cl, Z, yv, data[[idx[1]]], data[[idx[1]]], k)
  res <- rbind(
    vcompare(label, dname, "coef", coef(p), list(survival = coef(cl)[xnames]), tol),
    vcompare(label, dname, "time_effects", tp, list(survival = tcl), tol),
    vcompare(label, dname, "se_classical", se(vcov(p)),
             list(survival = se(vcov(cl))[xnames]), tol),
    vcompare(label, dname, "se_time_effects", p$time_effects_se[-1],
             list(survival = se(vcov(cl))[-seq_len(k)]), tol),
    vcompare(label, dname, "logLik", as.numeric(logLik(p)),
             list(survival = cl$loglik[2]), tol),
    vcompare(label, dname, "se_HC1", se(vcov(p, type = "HC1")),
             list(brute_force = rob$HC1), tol),
    vcompare(label, dname, "se_cluster", se(vcov(p, type = "cluster")),
             list(brute_force = rob$cluster), tol)
  )
  if (!is.null(cluster_coarse)) {
    rob2 <- brute_robust(cl, Z, yv, data[[idx[1]]], cluster_coarse, k)
    res <- rbind(res, vcompare(label, dname, "se_cluster_coarse",
                               se(vcov(p, type = "cluster", cluster = cluster_coarse)),
                               list(brute_force = rob2$cluster), tol))
  }
  stopifnot(attr(logLik(p), "df") == length(coef(cl)))
  res
}

lab <- "FE logit two-way (conditional)"
out$smoking <- check(lab, sname, smoking ~ price + income + stress, smk, c("id", "wave"),
                     c("price", "income", "stress"),
                     cluster_coarse = (smk$id - 1) %/% 5)
out$copd <- check(lab, cname, hospitalized ~ crp, copd, c("id", "visit"), "crp")
out$union <- check(lab, uname, union01 ~ wage + married + rural01 + health01, uw,
                   c("id", "year"), c("wage", "married", "rural01", "health01"))

# ---- T = 2: differenced logit with intercept --------------------------------
s2 <- smk[smk$wave %in% 2:3, ]
s2 <- s2[s2$id %in% names(which(table(s2$id) == 2)), ]
p2 <- suppressMessages(panglm(smoking ~ price + income + stress, s2, c("id", "wave"),
                              model = "within", family = "binomial", effect = "twoways"))
w1 <- s2[s2$wave == 2, ]; w2 <- s2[s2$wave == 3, ]
stopifnot(identical(w1$id, w2$id))
fd <- data.frame(y = w2$smoking, dprice = w2$price - w1$price,
                 dincome = w2$income - w1$income, dstress = w2$stress - w1$stress)
fd <- fd[w1$smoking != w2$smoking, ]
gd <- glm(y ~ dprice + dincome + dstress, binomial, fd,
          control = glm.control(epsilon = 1e-14, maxit = 100))
b <- c("dprice", "dincome", "dstress")
lab2 <- "FE logit two-way (conditional), T = 2"
d2 <- "simulated smoking, waves 2-3"
out$t2 <- rbind(
  vcompare(lab2, d2, "coef", coef(p2), list(glm_fd = coef(gd)[b]), tol),
  vcompare(lab2, d2, "time_effects", p2$time_effects[2], list(glm_fd = coef(gd)[1]), tol),
  vcompare(lab2, d2, "se_classical", se(vcov(p2)), list(glm_fd = se(vcov(gd))[b]), tol),
  vcompare(lab2, d2, "se_HC1", se(vcov(p2, type = "HC1")),
           list(sandwich = se(vcovHC(gd, type = "HC1"))[b]), tol),
  vcompare(lab2, d2, "se_cluster", se(vcov(p2, type = "cluster")),
           list(sandwich = se(vcovCL(gd, cluster = seq_len(nrow(fd)), type = "HC1"))[b]), tol)
)

do.call(vsave, c(unname(out), name = "twowaylogit"))
