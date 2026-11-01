#' Synthetic COPD follow-up panel with costs and a staggered programme
#'
#' A simulated longitudinal (panel) dataset of chronic obstructive pulmonary
#' disease (COPD) patients followed over five clinic visits, used
#' throughout this package's examples and vignettes. It is entirely
#' synthetic (not real patient data), generated to have a realistic panel
#' structure (time-invariant baseline covariates plus individual-level
#' frailties or random intercepts). Because every data-generating value is
#' known (see Details), it is the dataset on which students can check that
#' an estimator recovers the truth before using it on the real panels
#' shipped with the package ([docvisits], [meps], [medicaid],
#' [healthspend], [cigar]). The generating code is in `data-raw/copd.R`
#' in the package's GitHub repository.
#'
#' The first ten columns (`id` to `hospitalized`) are unchanged from
#' earlier releases. The last six (`clinic` to `symptoms`) add health
#' economics outcomes and were generated afterwards in a separate random
#' number stream, so the clinical columns did not change.
#'
#' @details
#' Data-generating values (patient random intercepts are normal with mean
#' zero and are independent of the covariates; `visit` enters linearly).
#'
#' * `exacerbations`: Poisson with mean
#'   `f_i exp(0.4 - 0.35 treatment + 0.02 (age - 65) + 0.3 smoker + 0.05 crp)`,
#'   with a gamma frailty `f_i` of mean 1 and variance 0.25.
#' * `hospitalized`: logit with index
#'   `-1.5 - 0.6 treatment + 0.03 (age - 65) + 0.5 smoker + 0.08 crp + r_i`,
#'   `sd(r_i) = 0.4`.
#' * `fev1`: `3.2 - 0.02 (age - 65) - 0.25 smoker + 0.15 treatment - 0.03 crp + r_i + e`,
#'   `sd(e) = 0.2` (same `r_i`).
#' * `cost`, any cost: probit with index
#'   `-0.6 + 0.4 exacerbations + 1.5 hospitalized + 0.3 smoker - 0.2 treatment + 0.05 crp + a_i`,
#'   `sd(a_i) = 0.5`.
#' * `cost`, positive part: lognormal, `log(cost) = 5.8 + 0.1 exacerbations
#'   + 0.7 hospitalized + 0.25 smoker - 0.2 treatment + 0.04 crp + 0.05 visit
#'   - 0.5 programme + c_k + b_i + e`, with `sd(b_i) = 0.4`,
#'   `corr(a_i, b_i) = 0.5`, `sd(e) = 0.4` and clinic effects
#'   `c_k = (0.20, 0.10, 0.15, 0.00, -0.10, -0.15)` for clinics 1 to 6.
#'   The programme therefore lowers expected cost among users by
#'   `1 - exp(-0.5)`, about 39 percent, without changing the probability of
#'   any cost. Clinics 1 and 2 adopt at visit 3, clinics 3 and 4 at visit 4,
#'   clinics 5 and 6 never; adopting clinics have higher baseline costs
#'   (selection on levels) but parallel trends hold. Because both potential
#'   costs were simulated from the same draws, the sample average effect on
#'   the treated rows is known: about -406.0 euros per visit interval (100 treated patient-visits).
#' * `adherence`: beta with precision 25 and logit mean
#'   `0.8 + 0.6 treatment - 0.5 smoker - 0.08 visit + v_i`, `sd(v_i) = 0.3`.
#' * `symptoms`: ordered logit with latent index
#'   `0.5 smoker - 0.6 treatment + 0.25 crp + w_i`, `sd(w_i) = 0.8`, and
#'   cut points 0.5, 2.0 and 3.5.
#'
#' Random-effects (conditional) estimators target these values; pooled
#' nonlinear fits estimate population-averaged coefficients, which are
#' attenuated for the binary, ordered and beta outcomes. With only 60
#' patients the estimates are noisy: `data-raw/copd.R` also refits the
#' correctly specified models on a 6,000-patient panel drawn from the same
#' design and checks that they land close to these values.
#'
#' @section Teaching use:
#' * Lung function and inflammation (`fev1` on `crp`): pooled, within,
#'   between and random effects, CRE and [panhausman()].
#' * Exacerbation counts: `family = "poisson"` or `"negbin"`, [panhurdle()],
#'   [panfmm()].
#' * Hospital admission: `family = "binomial"` or `"probit"`, with
#'   `dynamic = TRUE` for state dependence.
#' * What drives outpatient and inpatient cost? [twopart()] (probit plus
#'   lognormal or gamma, `correlated = TRUE`), `family = "gamma"` on
#'   positive costs, `family = "tobit"`, [panqr()] for cost quantiles.
#' * Did the care-management programme lower costs? [attgt()], [etwfe()],
#'   [eventstudy()] with `treat = "adopt"`, against the known effect above.
#' * Adherence (`family = "beta"`) and symptom severity (`family = "ologit"`,
#'   within = fixed-effects ordered logit).
#'
#' @format A data frame with 300 rows (60 patients x 5 visits) and 16 columns;
#'   the panel index is `c("id", "visit")`:
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
#'   \item{clinic}{integer clinic (1-6) treating the patient, time-invariant;
#'     10 patients per clinic}
#'   \item{adopt}{integer visit at which the patient's clinic adopted the
#'     care-management programme (3 or 4), 0 if the clinic never adopted;
#'     the `treat` column for [attgt()], [etwfe()] and [eventstudy()]}
#'   \item{programme}{integer indicator (0/1) that the programme is in place
#'     at this visit (`adopt > 0` and `visit >= adopt`)}
#'   \item{cost}{numeric outpatient plus inpatient cost (euros, two
#'     decimals) since the previous visit; 28.7 percent zeros and a
#'     right-skewed positive part}
#'   \item{adherence}{numeric proportion of prescribed inhaler doses taken
#'     since the previous visit, in (0, 1), three decimals}
#'   \item{symptoms}{ordered factor of symptom severity at the visit:
#'     none < mild < moderate < severe}
#' }
#' @source Simulated; see `data-raw/copd.R` in the package's GitHub
#'   repository.
#' @keywords datasets
#' @examples
#' data(copd)
#' str(copd)
#' mean(copd$cost == 0)
#' # two-part model of cost, programme effect on the positive part
#' # (clinic dummies absorb the clinic cost levels that differ between
#' # adopting and never-adopting clinics)
#' tp <- twopart(cost ~ exacerbations + hospitalized + smoker + treatment +
#'                 crp + visit + programme + factor(clinic),
#'               binary_formula = ~ exacerbations + hospitalized + smoker +
#'                 treatment + crp,
#'               data = copd, index = c("id", "visit"), positive = "lognormal")
#' summary(tp)
#' # staggered adoption: group-time effects of the programme on cost
#' att <- attgt(cost ~ 1, data = copd, index = c("id", "visit"),
#'              treat = "adopt", control = "notyet")
#' summary(att)
"copd"
