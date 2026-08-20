# Synthetic COPD follow-up panel used in examples and the vignette.
#
# Not real patient data -- simulated to have realistic panel structure
# (unbalanced-friendly, time-invariant baseline covariates, an
# individual-level frailty) so it exercises pooled/within/random effects
# across all four families without pulling in a Suggests package.

set.seed(20240501)

n_patients <- 60
n_visits <- 5

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
# visit -- the one time-varying covariate, needed so fixed-effects
# ("within") fits have something left to estimate after demeaning out the
# time-invariant covariates above
crp <- pmax(round(5 + 0.4 * (visit - 1) - 1.5 * treatment + rnorm(n_patients * n_visits, sd = 2), 1), 0.1)

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

copd <- data.frame(
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

usethis::use_data(copd, overwrite = TRUE)
