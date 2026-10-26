# MM-QR analytic covariance: cluster-robust (the default, vcov = "analytic")
# against the iid formula of Machado and Santos Silva (2019) (vcov = "iid").
# 1. Monte Carlo (mcstatsim::runsim): location-scale panel whose error U_it
#    is AR(1) within the individual, with regressors that also have an
#    individual component, N in {200, 1000}, T = 5, tau in {0.25, 0.5, 0.9},
#    500 replications (PANGLM_MC_REPS overrides). Coverage of 95% Wald
#    intervals and the ratio of the mean SE to the Monte Carlo SD, for the
#    cluster and iid covariances.
# 2. One cluster per observation on pglm::HealthIns (log(1 + medical
#    expenses)): the cluster sandwich keeps the observation-level products
#    that the iid formula replaces by sample means, so the two differ only
#    through heteroskedasticity or dependence between U and the regressors.
#    The ratios are reported; they are not expected to be exactly one.
# 3. Cluster SE against the cluster bootstrap over individuals (boot = 500)
#    on HealthIns: ratio near one.
source("validation/common.R")
suppressPackageStartupMessages({
  library(mcstatsim)
})
options(panglm.threads = 6)
MC_REPS <- as.integer(Sys.getenv("PANGLM_MC_REPS", "500"))
data("HealthIns", package = "pglm")
hi <- HealthIns
hi$lmed <- log1p(hi$med)
hi$kid <- as.numeric(hi$child == "yes")
hi <- hi[order(hi$id, hi$year), ]
hi <- hi[stats::ave(hi$id, hi$id, FUN = length) > 1, ]
hi$obs <- seq_len(nrow(hi))
taus <- c(0.25, 0.5, 0.75, 0.9)
f <- lmed ~ age + size + kid
idx <- c("id", "year")
tl <- function(tau) paste0("tau=", tau)
se <- function(m) sapply(m$vcov, function(v) sqrt(diag(v)))
out <- list()
add <- function(...) out[[length(out) + 1L]] <<- vcompare(...)

# ---- 1. Monte Carlo coverage -----------------------------------------------
# y_it = alpha_i + x_it + 0.5 z_it + (1 + 0.5 z_it) U_it, with
# U_it = 0.6 U_i,t-1 + sqrt(1 - 0.36) e_it stationary N(0, 1), x_it with an
# individual component correlated with alpha_i, z_it in [0, 2] with an
# individual component. True slopes: x 1, z 0.5 + 0.5 qnorm(tau).
mc_taus <- c(0.25, 0.5, 0.9)
truth <- function(t) c(x = 1, z = 0.5 + 0.5 * stats::qnorm(t))
sim_ar <- function(N, TT, rho = 0.6) {
  d <- data.frame(id = rep(seq_len(N), each = TT), year = rep(seq_len(TT), N))
  a <- rnorm(N)
  U <- matrix(0, N, TT)
  U[, 1] <- rnorm(N)
  for (t in 2:TT) U[, t] <- rho * U[, t - 1] + sqrt(1 - rho^2) * rnorm(N)
  d$x <- rep(0.5 * a + rnorm(N), each = TT) + rnorm(N * TT)
  d$z <- rep(runif(N), each = TT) + runif(N * TT)
  d$y <- rep(a, each = TT) + d$x + 0.5 * d$z + (1 + 0.5 * d$z) * as.vector(t(U))
  d
}
one_rep <- function(N) {
  d <- sim_ar(N, 5)
  m_c <- suppressWarnings(panqr(y ~ x + z, d, c("id", "year"), tau = mc_taus, model = "mmqr"))
  m_i <- suppressWarnings(panqr(y ~ x + z, d, c("id", "year"), tau = mc_taus, model = "mmqr",
                                vcov = "iid"))
  do.call(rbind, lapply(seq_along(mc_taus), function(j) {
    data.frame(tau = mc_taus[j], term = c("x", "z"), truth = truth(mc_taus[j]),
               est = m_c$coefficients[, j],
               se_cluster = sqrt(diag(m_c$vcov[[j]])),
               se_iid = sqrt(diag(m_i$vcov[[j]])))
  }))
}
# Run before any multithreaded fit in this process: forked workers deadlock
# if the parent has already started the parallel thread pool. Each worker
# fits single-threaded (6 workers in all).
options(panglm.threads = 1)
reps <- runsim(MC_REPS, data.frame(N = c(200, 1000)), one_rep, show_progress = FALSE,
               num_cores = 6, seed = 20261026)
options(panglm.threads = 6)
cover <- function(est, se, truth) mean(abs(est - truth) <= stats::qnorm(0.975) * se)
summ <- do.call(rbind, lapply(split(reps, list(reps$N, reps$tau, reps$term)), function(r) {
  data.frame(N = r$N[1], tau = r$tau[1], term = r$term[1],
             bias = mean(r$est - r$truth), mc_sd = stats::sd(r$est),
             se_cluster_ratio = mean(r$se_cluster) / stats::sd(r$est),
             se_iid_ratio = mean(r$se_iid) / stats::sd(r$est),
             cover_cluster = cover(r$est, r$se_cluster, r$truth),
             cover_iid = cover(r$est, r$se_iid, r$truth),
             cover_cluster_mean = cover(r$est, r$se_cluster, mean(r$est)))
}))
summ <- summ[order(summ$N, summ$tau, summ$term), ]
rownames(summ) <- NULL
cat("Monte Carlo, AR(1) errors, T = 5,", MC_REPS, "replications\n")
print(summ, digits = 3)
lab <- paste0("N=", summ$N, " ", tl(summ$tau), " ", summ$term)
# MM-QR is biased for the scale slope at fixed T (the within step shrinks
# the residual quantiles; the bias is the same at N = 200 and N = 1000), so
# coverage of the true value is reported only. The covariance is checked by
# the ratio of the mean SE to the Monte Carlo SD and by coverage of the
# Monte Carlo mean of the estimates.
mc_lab <- paste0("Monte Carlo AR(1) T=5 (", MC_REPS, " reps)")
add("panqr mmqr", mc_lab, "mean cluster SE / MC SD",
    stats::setNames(summ$se_cluster_ratio, lab), list(one = rep(1, nrow(summ))), 0.1)
add("panqr mmqr", mc_lab, "coverage of the MC mean, cluster",
    stats::setNames(summ$cover_cluster_mean, lab), list(nominal = rep(0.95, nrow(summ))), 0.03)


# ---- 2. one cluster per observation versus iid ------------------------------
m_iid <- suppressWarnings(panqr(f, hi, idx, tau = taus, model = "mmqr", vcov = "iid"))
m_obs <- suppressWarnings(panqr(f, hi, idx, tau = taus, model = "mmqr", cluster = "obs"))
m_cl <- suppressWarnings(panqr(f, hi, idx, tau = taus, model = "mmqr"))
ratio_obs <- se(m_obs) / se(m_iid)
ratio_id <- se(m_cl) / se(m_iid)
cat("HealthIns: SE with one cluster per observation / iid SE\n")
print(round(ratio_obs, 3))
cat("HealthIns: SE clustered by individual / iid SE\n")
print(round(ratio_id, 3))
# report only: heteroskedasticity makes these differ from one
add("panqr mmqr", "HealthIns", "SE obs-cluster / iid", as.vector(ratio_obs),
    list(one = rep(1, length(ratio_obs))), 0.2)

# ---- 3. cluster SE versus cluster bootstrap over individuals ----------------
m_bt <- suppressWarnings(panqr(f, hi, idx, tau = taus, model = "mmqr",
                               vcov = "bootstrap", boot = 500, seed = 20261026))
ratio_bt <- se(m_cl) / se(m_bt)
cat("HealthIns: cluster SE / bootstrap SE (boot = 500)\n")
print(round(ratio_bt, 3))
for (j in seq_along(taus)) {
  add("panqr mmqr", "HealthIns", paste("cluster SE / bootstrap SE", tl(taus[j])),
      ratio_bt[, j], list(one = rep(1, nrow(ratio_bt))), 0.1)
}

do.call(vsave, c(out, name = "qr-cluster"))
