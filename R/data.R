#' Synthetic COPD follow-up panel
#'
#' A simulated longitudinal (panel) dataset of chronic obstructive pulmonary
#' disease (COPD) patients followed over several clinic visits, used
#' throughout this package's examples and vignette. It is entirely
#' synthetic (not real patient data), generated to have a realistic panel
#' structure (time-invariant baseline covariates plus an individual-level
#' frailty or random intercept), so it exercises pooled, fixed-effects and
#' random-effects fits for continuous, count and binary outcomes without
#' depending on any suggested package. The generating code is in
#' `data-raw/copd.R` in the package's GitHub repository.
#'
#' @format A data frame with 300 rows (60 patients x 5 visits) and 10 columns:
#' \describe{
#'   \item{id}{integer patient identifier (1-60), the panel individual index}
#'   \item{visit}{integer visit number (1-5), the panel time index}
#'   \item{treatment}{integer treatment arm, time-invariant (0 = control,
#'     1 = active)}
#'   \item{age}{numeric patient age in whole years at baseline,
#'     time-invariant}
#'   \item{smoker}{integer current-smoker indicator, time-invariant (0/1)}
#'   \item{bmi}{numeric body mass index (kg/m^2) at baseline, one decimal,
#'     time-invariant; not used to generate the outcomes}
#'   \item{crp}{numeric C-reactive protein (mg/L) measured at each visit,
#'     one decimal and at least 0.1; the one time-varying covariate, needed
#'     for fixed-effects ("within") fits to have something left to estimate
#'     after demeaning}
#'   \item{fev1}{numeric forced expiratory volume in 1 second (litres), two
#'     decimals and at least 0.5; a continuous lung-function outcome, one
#'     value per visit}
#'   \item{exacerbations}{integer number of COPD exacerbations since the
#'     previous visit; a count outcome, one value per visit}
#'   \item{hospitalized}{integer indicator (0/1) of a hospital admission in
#'     the interval since the previous visit; a binary outcome, one value per
#'     visit}
#' }
#' @source Simulated; see `data-raw/copd.R` in the package's GitHub
#'   repository.
#' @examples
#' data(copd)
#' str(copd)
#' table(copd$visit)
"copd"
