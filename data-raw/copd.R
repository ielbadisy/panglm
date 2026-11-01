# Synthetic COPD follow-up panel used in examples, tests and the vignettes.
#
# Not real patient data: simulated to have a realistic panel structure
# (time-invariant baseline covariates, an individual-level frailty) so it
# exercises pooled, within and random effects across the families without
# pulling in a Suggests package.
#
# The file has two parts.
#
# 1. simulate_copd(): the original 10 clinical columns (id, visit,
#    treatment, age, smoker, bmi, crp, fev1, exacerbations, hospitalized),
#    generated exactly as in panglm <= 1.7.0 with set.seed(20240501), so the
#    shipped values (and every example and test using them) are unchanged.
#
# 2. add_economics(): health-economics columns added in a SEPARATE random
#    number stream (set.seed(20261101)) after the clinical columns:
#    clinic, adopt, programme, cost, adherence, symptoms. Their
#    data-generating values are listed in `truth` below and in the help
#    page, so students can check that the estimators recover them.
#
# Running the script rebuilds data/copd.rda, checks that the clinical
# columns are identical to the previous release, and (if panglm is
# installed) checks on a large simulated panel (6,000 patients) that the
# estimators recover the generating values.

simulate_copd <- function(n_patients = 60, n_visits = 5, seed = 20240501) {
  set.seed(seed)

  id <- rep(seq_len(n_patients), each = n_visits)
  visit <- rep(seq_len(n_visits), times = n_patients)

  treatment <- rep(rbinom(n_patients, 1, 0.5), each = n_visits)
  age <- rep(round(rnorm(n_patients, mean = 65, sd = 8)), each = n_visits)
  smoker <- rep(rbinom(n_patients, 1, 0.35), each = n_visits)
  bmi <- rep(round(rnorm(n_patients, mean = 27, sd = 4), 1), each = n_visits)

  # individual-level frailty (unobserved heterogeneity) shared across a
  # patient's visits, used to generate both the count and binary outcomes
  frailty <- rep(rgamma(n_patients, shape = 4, rate = 4), each = n_visits)
  re_int <- rep(rnorm(n_patients, sd = 0.4), each = n_visits)

  # crp: C-reactive protein (mg/L), an inflammation marker measured at each
  # visit, the one time-varying covariate, needed so fixed-effects
  # ("within") fits have something left to estimate after demeaning out the
  # time-invariant covariates above
  crp <- pmax(round(5 + 0.4 * (visit - 1) - 1.5 * treatment +
                      rnorm(n_patients * n_visits, sd = 2), 1), 0.1)

  # fev1: forced expiratory volume in 1s (litres), a continuous lung-function
  # outcome that declines with age/smoking/inflammation and improves under
  # treatment
  fev1 <- 3.2 - 0.02 * (age - 65) - 0.25 * smoker + 0.15 * treatment - 0.03 * crp +
    re_int + rnorm(n_patients * n_visits, sd = 0.2)
  fev1 <- round(pmax(fev1, 0.5), 2)

  # exacerbations: count of COPD exacerbations since the previous visit
  exacerbations <- rpois(
    n_patients * n_visits,
    lambda = frailty * exp(0.4 - 0.35 * treatment + 0.02 * (age - 65) + 0.3 * smoker + 0.05 * crp)
  )

  # hospitalized: binary indicator of a hospital admission in the interval
  lin_pred <- -1.5 - 0.6 * treatment + 0.03 * (age - 65) + 0.5 * smoker + 0.08 * crp + re_int
  p_hosp <- 1 / (1 + exp(-lin_pred))
  hospitalized <- rbinom(n_patients * n_visits, 1, p_hosp)

  data.frame(
    id = id,
    visit = visit,
    treatment = treatment,
    age = age,
    smoker = smoker,
    bmi = bmi,
    crp = crp,
    fev1 = fev1,
    exacerbations = exacerbations,
    hospitalized = hospitalized
  )
}

# Data-generating values of the economic columns (documented in ?copd).
truth <- list(
  # clinics and the staggered care-management programme
  n_clinics = 6L,
  adopt = c(`1` = 3L, `2` = 3L, `3` = 4L, `4` = 4L, `5` = 0L, `6` = 0L),
  clinic_effect = c(0.20, 0.10, 0.15, 0.00, -0.10, -0.15),
  # any cost: probit, P(cost > 0) = pnorm(index + a_i)
  any = c(`(Intercept)` = -0.60, exacerbations = 0.40, hospitalized = 1.50,
          smoker = 0.30, treatment = -0.20, crp = 0.05),
  sd_a = 0.50,
  # positive cost: lognormal, log(cost) = index + b_i + e, e ~ N(0, sigma_e^2)
  positive = c(`(Intercept)` = 5.80, exacerbations = 0.10, hospitalized = 0.70,
               smoker = 0.25, treatment = -0.20, crp = 0.04, visit = 0.05,
               programme = -0.50),
  sigma_e = 0.40,
  sd_b = 0.40,
  rho_ab = 0.50,
  # adherence: beta regression, logit mean = index + v_i, precision 25
  adherence = c(`(Intercept)` = 0.80, treatment = 0.60, smoker = -0.50, visit = -0.08),
  sd_v = 0.30,
  phi = 25,
  # symptoms: ordered logit, latent = index + w_i + logistic error
  symptoms = c(smoker = 0.50, treatment = -0.60, crp = 0.25),
  sd_w = 0.80,
  cuts = c(0.5, 2.0, 3.5)
)

add_economics <- function(d, seed = 20261101, tr = truth) {
  set.seed(seed)
  ids <- unique(d$id)
  n <- length(ids)
  nv <- as.integer(table(d$id)[1])
  N <- nrow(d)
  ind <- match(d$id, ids)

  # clinic: patients spread evenly over 6 clinics at random
  clinic_i <- sample(rep(seq_len(tr$n_clinics), length.out = n))
  clinic <- clinic_i[ind]
  adopt <- unname(tr$adopt[clinic])
  programme <- as.integer(adopt > 0 & d$visit >= adopt)

  # correlated patient random intercepts of the two cost parts
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  a <- (tr$sd_a * z1)[ind]
  b <- (tr$sd_b * (tr$rho_ab * z1 + sqrt(1 - tr$rho_ab^2) * z2))[ind]
  v <- rnorm(n, sd = tr$sd_v)[ind]
  w <- rnorm(n, sd = tr$sd_w)[ind]

  xa <- cbind(1, d$exacerbations, d$hospitalized, d$smoker, d$treatment, d$crp)
  any <- rbinom(N, 1, pnorm(drop(xa %*% tr$any) + a))

  xb <- cbind(1, d$exacerbations, d$hospitalized, d$smoker, d$treatment, d$crp,
              d$visit, 0)
  mu0 <- exp(drop(xb %*% tr$positive) + tr$clinic_effect[clinic] + b)
  draw <- exp(rnorm(N, sd = tr$sigma_e))
  cost0 <- any * mu0 * draw                                  # without programme
  cost1 <- any * mu0 * exp(tr$positive[["programme"]]) * draw  # with programme
  cost <- round(ifelse(programme == 1, cost1, cost0), 2)

  eta_v <- drop(cbind(1, d$treatment, d$smoker, d$visit) %*% tr$adherence) + v
  m <- 1 / (1 + exp(-eta_v))
  adherence <- rbeta(N, m * tr$phi, (1 - m) * tr$phi)
  adherence <- pmin(pmax(round(adherence, 3), 0.001), 0.999)

  lat <- drop(cbind(d$smoker, d$treatment, d$crp) %*% tr$symptoms) + w + rlogis(N)
  symptoms <- cut(lat, c(-Inf, tr$cuts, Inf),
                  labels = c("none", "mild", "moderate", "severe"),
                  ordered_result = TRUE)

  out <- d
  out$clinic <- as.integer(clinic)
  out$adopt <- as.integer(adopt)
  out$programme <- programme
  out$cost <- cost
  out$adherence <- adherence
  out$symptoms <- symptoms
  # sample average treatment effect on the treated (euros), known because
  # both potential costs are simulated from the same draws
  attr(out, "satt") <- mean((cost1 - cost0)[programme == 1])
  out
}

copd <- simulate_copd()

# the clinical columns must be identical to the previous release
if (file.exists("data/copd.rda")) {
  old <- new.env()
  load("data/copd.rda", envir = old)
  stopifnot(identical(old$copd[1:10], copd))
}

copd <- add_economics(copd)
satt <- attr(copd, "satt")
attr(copd, "satt") <- NULL

cat("share of zero costs:", round(mean(copd$cost == 0), 3), "\n")
cat("sample ATT of the programme on cost (euros):", round(satt, 1), "\n")
print(table(copd$symptoms))
print(summary(copd$cost[copd$cost > 0]))
stopifnot(mean(copd$cost == 0) > 0.25, mean(copd$cost == 0) < 0.35)

# Recovery of the generating values on a large panel drawn from the same
# design: the correctly specified estimator should land close to the truth.
if (requireNamespace("panglm", quietly = TRUE)) {
  library(panglm)
  idx <- c("id", "visit")
  big <- add_economics(simulate_copd(n_patients = 3000), seed = 1)
  big_satt <- attr(big, "satt")
  big$fclinic <- factor(big$clinic)
  check <- function(est, true, tol, what) {
    cat(sprintf("%-32s estimate %9.3f  truth %9.3f\n", what, est, true))
    stopifnot(abs(est - true) < tol)
  }

  # binary part: random-effects probit for any cost
  any <- panglm(I(cost > 0) ~ exacerbations + hospitalized + smoker + treatment + crp,
                data = big, index = idx, model = "random", family = "probit")
  for (v in c("exacerbations", "hospitalized", "smoker", "crp"))
    check(coef(any)[[v]], truth$any[[v]], 0.08, paste("any cost (RE probit):", v))

  # positive part alone, conditional on cost > 0 (selection on a_i only
  # shifts the intercept a little, slopes stay close to the truth)
  pos <- panglm(log(cost) ~ exacerbations + hospitalized + crp + visit + programme +
                  fclinic, data = big[big$cost > 0, ], index = idx, model = "random")
  for (v in c("exacerbations", "hospitalized", "crp", "programme"))
    check(coef(pos)[[v]], truth$positive[[v]], 0.05, paste("log positive cost:", v))

  ad <- panglm(adherence ~ treatment + smoker + visit, data = big, index = idx,
               model = "random", family = "beta")
  for (v in names(truth$adherence))
    check(coef(ad)[[v]], truth$adherence[[v]], 0.05, paste("adherence:", v))

  sy <- panglm(symptoms ~ crp, data = big, index = idx, model = "within",
               family = "ologit")
  check(coef(sy)[["crp"]], truth$symptoms[["crp"]], 0.03, "symptoms (FE ologit): crp")

  att <- summary(attgt(cost ~ 1, data = big, index = idx, treat = "adopt",
                       control = "notyet"))$overall
  check(att$estimate, big_satt, 3 * att$se, "attgt overall ATT (euros)")
}

usethis::use_data(copd, overwrite = TRUE, compress = "xz")
