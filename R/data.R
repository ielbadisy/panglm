#' Synthetic COPD follow-up panel
#'
#' A simulated longitudinal (panel) dataset of chronic obstructive pulmonary
#' disease (COPD) patients followed over several clinic visits, used
#' throughout this package's examples and vignette. It is entirely
#' synthetic (not real patient data), generated to have a realistic panel
#' structure -- time-invariant baseline covariates plus an individual-level
#' frailty/random intercept -- so it exercises pooled, fixed-effects, and
#' random-effects fits across all four families (`gaussian`, `poisson`,
#' `binomial`, `negbin`) without depending on any Suggested package. See
#' `data-raw/copd.R` for the generating code.
#'
#' @format A data frame with 300 rows (60 patients x 5 visits) and 9 columns:
#' \describe{
#'   \item{id}{patient identifier (1-60)}
#'   \item{visit}{visit number (1-5), the panel time index}
#'   \item{treatment}{treatment arm, time-invariant (0 = control, 1 = active)}
#'   \item{age}{patient age in years at baseline, time-invariant}
#'   \item{smoker}{current smoker, time-invariant (0/1)}
#'   \item{bmi}{body mass index at baseline, time-invariant}
#'   \item{crp}{C-reactive protein (mg/L) measured at each visit; the one
#'     time-varying covariate, needed for fixed-effects ("within") fits to
#'     have something left to estimate after demeaning}
#'   \item{fev1}{forced expiratory volume in 1 second (litres); a continuous
#'     lung-function outcome, one value per visit}
#'   \item{exacerbations}{number of COPD exacerbations since the previous
#'     visit; a count outcome, one value per visit}
#'   \item{hospitalized}{hospital admission in the interval since the
#'     previous visit; a binary outcome, one value per visit}
#' }
#' @source Simulated; see `data-raw/copd.R`.
"copd"
