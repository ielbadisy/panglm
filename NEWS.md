# panglm 1.7.0

Plotting and teaching release: every fitted object has a plot method,
the diagnostics, marginal effects and difference-in-differences results
can be drawn directly with base graphics, and five real health economics
panels ship with the package. Estimates are unchanged.

## Diagnostics and coefficients

* `plot.panglm()` gains `which = "qq"` (normal Q-Q plot of the Pearson
  residuals), `"calibration"` (observed proportion against the mean
  predicted probability in quantile bins, for binary outcomes) and
  `"rootogram"` (hanging rootogram for Poisson and negative binomial
  fits, including random-intercept fits through the marginal count
  distribution). The residual view uses Pearson residuals by default
  (`type`), adds a smooth, and uses a log axis for log-link models.
  Coefficient plots take `parm`, omit the intercept by default and size
  the margin to the term names.
* `residuals.panglm()` gains `type = c("response", "pearson")`.
* New `coefplot()` compares the coefficients of several fitted models in
  one figure.
* New plot methods for `twopart()`, `panfmm()` (coefficients by class and
  the posterior class probabilities), `pangmm()`, `paniv()` and
  `panselect()` fits; `plot.panglm_hurdle()` gains `parm`.
* `plot.panglm()` gains `which = "density"`: the observed kernel density
  against the model-implied density, with the point masses of Tobit
  censoring and two-part zeros; also for `twopart()` fits.
* `plot.panglm()` gains `which = "effects"`: the distribution of the
  estimated individual effects against the assumed mixing distribution,
  with an intercept-slope scatter for random slopes.
* `which = "correlation"` draws the correlation matrix of the estimates
  as a heat map, for `panglm()`, `twopart()`, `panfmm()` and
  `panhurdle()` fits.

## Individual effects

* New `fixef()` returns the estimated fixed effects of within fits, and
  `ranef()` the posterior means or modes of the random effects of random
  fits, for random intercepts in every family and for random slopes.
  The methods also register on the `nlme` generics.

## Marginal effects

* New plot method for `partialeffects()` results.
* `effectcurve()` gains `scale = c("link", "response", "slope")`: the
  average predicted outcome, or the average marginal effect, as the
  regressor varies, with delta-method bands, also for factor regressors;
  and `by`, curves for each level of a factor. The plot adds a rug of the
  observed values.

## Difference-in-differences

* `attgt()`, `etwfe()` and `eventstudy()` plots share one event-study
  layout, with uniform confidence bands drawn as visible bars behind the
  pointwise intervals. All three take `min_e`, `max_e`, `ylim`, `main`
  and `legend`.
* Event-study plots and `eventplot()` gain `style = c("classic",
  "points", "line")`. The default classic style draws the estimates with
  capped intervals, the reference period as an open circle at zero and a
  dashed line labelled "Intervention" at event time 0; `"line"` joins
  the estimates over a shaded band. The arguments `ref`, `onset`,
  `onset_label` and `caps` control the markers.
* `plot(<attgt>)` and `plot(<etwfe>)` gain `type = "heatmap"`: the
  ATT(g, t) estimates by cohort and period on a diverging scale, with
  pre-treatment cells hatched and, with `cells = TRUE`, the estimates
  printed.
* `plot(<attgt>, type = "gt")` draws the group-time effects, one panel
  per cohort.
* New `treatmentplot()` shows the treatment timing of a staggered design.
* New `eventplot()` compares event studies from several estimators on
  one axis. Legends wrap into columns when the panel is narrow.
* New `trendplot()` draws the raw mean outcome by adoption cohort, by
  treated versus never treated, or on the event-time axis, with
  cluster-robust bands.
* New `panelplot()` draws a unit by period heat map of a panel: the
  observation pattern, a numeric or categorical variable, and units
  sorted by adoption cohort.

## Data

* `copd` gains clinics, staggered adoption of a care programme, annual
  costs with structural zeros, a medication adherence proportion and an
  ordered symptom score, all simulated with known parameters.
* New real panels for teaching: `docvisits` (German health care demand,
  Riphahn, Wambach and Million 2003), `meps` (two-year Medical
  Expenditure Panel Survey spending), `medicaid` (US states, Affordable
  Care Act Medicaid expansion and uninsurance, 2009-2019), `healthspend`
  (countries, health expenditure from the World Development Indicators)
  and `cigar` (US states, cigarette demand, Baltagi). Each help page
  states the source, the licence and the questions the data can teach.
* New vignette `vignette("panglm-datasets")`.
* Lazy-loaded data are compressed with xz.

## Bug fixes

* `panhausman(method = "classical")` no longer reports p = 1 when
  `V_fe - V_re` is not positive definite and the statistic is negative:
  the signed statistic is kept, the p-value is `NA` and a warning points
  to `method = "mundlak"`. The print method drops the stray space before
  commas.
* Print methods use one sample-size line, `n = <observations>  groups =
  <individuals>` (`pangmm()` printed the number of individuals as `N`),
  and keep every line within 80 characters (the `pangmm()` summary, the
  `attgt()` summary and the `etwfe()` footer were wider).
* `panglm()` detects columns that are collinear up to rounding, not only
  exactly. With `cre = TRUE`, an age that rises one per year but is
  stored with a small rounding error (e.g. `pglm::HealthIns`) gave
  offsetting coefficients of order 1e5 with no warning. Columns are now
  scaled by their own norm before the rank check, and a near dependency
  is dropped with the usual "rank-deficient model matrix" warning. Exact
  dependencies drop the same columns as `lm()`. This applies to every
  model and family.
* `pangmm(model = "onestep")` no longer warns that the second-step
  weighting matrix is singular: a one-step fit uses it only for the
  Hansen statistic, through a generalized inverse as in `plm::pgmm()`.
  Two-step fits keep the warning.
* `panhurdle(family = "negbin")` no longer stops with "NAs are not
  allowed in subscripted assignments" when the NB2 shape parameter runs
  off to infinity (counts no more dispersed than Poisson): `log P(Y = 0)`
  is computed as `-theta * log1p(mu / theta)`, which stays finite.
* `panhurdle()` drops time-invariant regressors with a warning, as
  `panglm(model = "within")` does, instead of failing with "system is
  computationally singular". Groups whose positive counts are all 1 (no
  finite intercept) are dropped from the count part with a message, and
  a binary outcome stops with a clear error.

## Performance

* `twopart(correlated = TRUE)` is much faster: the joint likelihood and
  its exact gradient (including the dependence of the adaptive quadrature
  mode and scale on the parameters) are computed in C++, replacing
  numerical differences of an R likelihood. A `meps` fit on 6,000
  person-years takes about 5 seconds instead of 5 minutes. Estimates are
  unchanged up to the optimizer tolerance, and the likelihood matches
  `GLMMadaptive::mixed_model(family = hurdle.lognormal())` at the same
  parameters.
* The `panhurdle()` count part is fit by Newton-Raphson with analytic
  derivatives and the group intercepts eliminated, instead of a dense
  dummy matrix (and, for NB2, BFGS with a numerical Hessian). A fit that
  took minutes on a few hundred `docvisits` individuals now takes well
  under a second, the full panel about one second. Estimates match
  `glmmTMB` truncated Poisson and NB2 fits with group dummies; the NB2
  standard errors, which came from a not fully converged fit, change.

* `panglm(model = "random", family = "tobit")` now converges on outcomes
  in raw units such as costs in dollars. The quasi-Newton search used to
  stop silently short of the maximum (on `meps$totexp` the log-likelihood
  was 29 units below it and `sigma` 22% off) because the slopes and the
  log standard deviations differed in scale by orders of magnitude. Both
  Tobit fits now maximize the likelihood of `y / c` (`c` a power of two
  near `sd(y)`) with standardized regressors and map the estimates back,
  so that coefficients, `sigma_e`, `sigma`, their covariance and the
  log-likelihood (which loses `log c` per uncensored observation) do not
  depend on the units of `y`. The random-effects likelihood and its exact
  scores are computed in C++: the `meps` fit takes 0.4 s instead of 210 s.

## Performance

* Random-intercept binomial, ordered and Tobit fits by adaptive quadrature are
  faster: parameter values with `sigma` outside `[1e-11, 7e10]`, which the
  first quasi-Newton steps can visit, are rejected without a mode search
  (each such evaluation used to cost as much as 200 ordinary ones),
  each individual's mode search starts from its mode at the previous
  evaluation and takes one pass over the rows per Newton step, and the
  binomial fit starts from a few BHHH steps. The probit model of `doctor`
  on `docvisits` (27,326 rows) takes 4 s instead of 37 s; the estimates
  agree to `1e-7`.

# panglm 1.6.2

Speed release for large panels. Estimates are unchanged: the
optimizations either reproduce the previous computations exactly or
change results only at the level of floating-point rounding, and the
tests and validation scripts pass unchanged.

## Performance

* Two-way demeaning computes the group means in one pass over the data
  instead of a parallel reduction that allocated one accumulator per
  thread. The two-way Gaussian fit demeans the outcome and the
  regressors once and reuses the result for the check of which columns
  are identified.
* The fixed-effects Poisson and random-effects Poisson and negative
  binomial Newton-Raphson fitters reuse the gradient and Hessian of the
  accepted step instead of evaluating the same point twice.
* Random-effects negative binomial: the gamma, digamma and trigamma
  differences of the likelihood are computed by their finite-sum
  recurrences for counts up to 50, skipped for zero counts, and spread
  over threads. A fit on 250,000 observations is about seven times
  faster.
* `pangmm()`: the individual sums of the GMM moments and the
  Windmeijer correction use block sums and two matrix products instead
  of one grouped sum per regressor.
* The check for duplicated (individual, time) pairs in `panglm()`,
  `pangmm()`, `attgt()` and `predict()` compares neighbouring sorted
  rows instead of pasting every row into a string.
* `attgt()` and `did()`: the checks that weights, clusters and the
  treatment timing are constant within individuals no longer call R once
  per individual.
* `attgt()` drops the row names of the outcome and design before
  subsetting them, and fits the propensity score with a compiled
  logistic IRLS that follows `glm.fit()` step by step. Unusual fits
  (nearly collinear design, no convergence, fitted probabilities of 0
  or 1) still go through `glm.fit()`. On 50,000 individuals the
  `"dr"` and `"ipw"` estimators are about twice as fast.
* Gaussian fits with the identity link are computed by one weighted
  least squares solve instead of running the IRLS loop to convergence.
* `panglm()` no longer copies the outcome and design when the data are
  already sorted by individual and time.

## Bug fixes

* `attgt()`: when the individual index is a factor, individuals whose
  rows had all been removed earlier (for example for missing values) were
  also counted among the individuals dropped to balance the panel
  (warning message and `dropped$unbalanced`). Only individuals still
  present are counted now. Estimates are unaffected.

# panglm 1.6.1

Correctness release. Every fix is covered by the tests and the
validation scripts.

## Changed behaviour

* Random-intercept binomial (logit, probit) and ordered (logit, probit)
  models are now integrated by adaptive Gauss-Hermite quadrature centred
  at each individual's mode, instead of non-adaptive quadrature.
  Estimates change, most for large random-effect variances or long
  panels. At the default `R = 21` the fits match
  `lme4::glmer(nAGQ = 21)` and `ordinal::clmm(nAGQ = 21)` to about 1e-5;
  `R = 1` is now the Laplace approximation (glmmTMB, `clmm(nAGQ = 1)`).
* `panqr(model = "mmqr")`: the analytic covariance (`vcov = "analytic"`,
  the default) is now cluster-robust by individual. The influence
  functions of the location, scale and quantile steps are summed within
  each cluster and `V = G/(G-1) sum_g psi_g psi_g'`, where `G` is the
  number of clusters. The formula of Machado and Santos Silva (2019),
  which assumes independent observations and can badly overstate the
  standard errors in short panels, is kept as `vcov = "iid"`. A new
  `cluster` argument names another clustering column. Standard errors
  of existing MM-QR fits change.
* Gaussian random effects on unbalanced panels now use the Swamy-Arora
  variance components of plm (`random.method = "swar"`), so estimates on
  unbalanced panels change and match `plm::plm(model = "random")`.
  Balanced fits are unchanged. A negative individual-effect variance is
  set to zero with a warning.
* The one-way and two-way fixed-effects logit are fitted by
  Newton-Raphson with the exact analytic Hessian of the conditional
  likelihood instead of BFGS with a finite-difference Hessian.
  Coefficients and standard errors match `survival::clogit()` to about
  1e-13, and the fit is faster than clogit at N = 10,000 and 50,000.

## Bug fixes

* `panglm()` accepts base R family objects (`gaussian()`, `poisson()`,
  `binomial("logit")`, `binomial("probit")`, `Gamma("log")`) and maps
  them onto the panglm families; other stats families stop with a
  message naming the panglm constructor.
* `attgt()`, `did()`, `eventstudy()` and `etwfe()` no longer fail when
  the individual index is a factor with unused levels (e.g. after
  `subset()`).
* `pangmm()` drops regressors that are not identified in the transformed
  equation (for example age in a balanced annual panel) with a warning,
  instead of stopping with a singular-system error.
* Two-way demeaning (`effect = "twoways"`) by alternating projections
  now warns if it reaches the iteration limit before converging.

## Documentation

* Two-way demeaning is described as iterative, converging to tolerance,
  rather than exact.
* The `cre` test is described as a Hausman-type (Mundlak) test that is
  robust only with `vcov = "cluster"` (or `"HC1"` for
  heteroskedasticity alone).
* The beta family validation is now in `validation/validate-beta.R`
  (betareg, GLMMadaptive, glmmTMB and numerical integration).
* The README, vignettes and help pages describe what each function fits
  rather than listing the packages it is checked against; the checks are
  summarised in the validation article on the package website.

# panglm 1.6.0

This release collects all changes since 1.1.4, the previous CRAN version.
Every new estimator is covered by the tests and the validation scripts;
the validation article on the package website summarises the checks.

## Deprecations and changed behaviour

* Exported names no longer use underscores. The families are now
  `gaussianfamily()`, `poissonfamily()`, `binomialfamily()`,
  `negbinfamily()`, `gammafamily()`, `betafamily()`, `orderedfamily()` and
  `tobitfamily()`; `panglm_hausman()`, `panglm_dispersiontest()` and
  `panglm_hurdle()` are now `panhausman()`, `pandispersiontest()` and
  `panhurdle()`. The old names still work and warn once per session; see
  `?"panglm-deprecated"`. Result classes are unchanged.
* `panhurdle()` takes `family =` instead of `count_family =`, and
  `panqr()` takes `boot =` instead of `R =` and gains `seed`. The old
  arguments still work with a deprecation warning.
* `predict.panfmm(type = "class")` now returns the modal class label; the
  class-specific means it returned before are `type = "components"`. A
  message flags the change once per session.

## New estimators

* `model = "between"`: the between (individual means) estimator for the
  Gaussian family (plm).
* `tobitfamily(left, right)`: censored outcomes, pooled or with a random
  intercept (AER, censReg).
* `orderedfamily()`: ordered logit and probit, pooled, random intercept, and
  the fixed-effects ordered logit by "blow-up and cluster" (Baetschmann,
  Staub and Winkelmann 2015).
* `betafamily()`: beta regression for outcomes in (0, 1), pooled or with a
  random intercept (betareg, glmmTMB).
* `gammafamily()`: gamma regression with log link for costs, pooled, fixed
  effects or random intercept.
* Two-way fixed-effects logit: `model = "within", family = "binomial",
  effect = "twoways"` (survival::clogit).
* Random slopes: `panglm(model = "random", random = ~ time)`, correlated
  random slopes with an unstructured covariance. Gaussian fits use exact ML
  with any number of slopes (lme4); Poisson, logit, probit, negative
  binomial and gamma fits use adaptive Gauss-Hermite quadrature with any
  number of slopes (glmmTMB, GLMMadaptive).
* `paniv()`: panel instrumental variables (fixed effects, G2SLS or EC2SLS
  random effects, pooled 2SLS) and the Hausman-Taylor estimator, with
  first-stage F statistics and Sargan/Hansen tests (plm).
* `pangmm()`: Arellano-Bond difference and Blundell-Bond system GMM, one- or
  two-step with Windmeijer standard errors, AR and Hansen tests, collapsed
  instruments and lag limits; it warns when instruments outnumber
  individuals (plm::pgmm).
* `panselect()`: the Wooldridge (1995) correction for attrition or death in
  cohort studies (glm, plm, sampleSelection).
* `panqr()`: panel quantile regression by Canay (2011), MM-QR (Machado and
  Santos Silva 2019) or pooling, with a `plot()` method (quantreg).
* `panfmm()` and `panfmmselect()`: latent class GLMs with class membership
  constant within individual, Poisson, NB2, Gaussian, logit and probit
  components, weights, offsets and a concomitant model (flexmix).
* `twopart()`: two-part models for semicontinuous costs, with correlated
  random intercepts (`correlated = TRUE`).
* `did()` and `eventstudy()`: two-way fixed-effects difference-in-differences
  and event studies for Gaussian, Poisson, negative binomial and binomial
  (conditional logit) outcomes (fixest, survival::clogit).
* `attgt()`: Callaway and Sant'Anna (2021) group-time ATTs with
  aggregations, pre-trend test, weights, extra clustering and a multiplier
  bootstrap with uniform bands (did).
* `etwfe()`: the Wooldridge (2023) extended two-way FE estimator for linear,
  Poisson and logit outcomes, with weights, `xvar` and `ivar` (etwfe).

## New features

* `lag()`, `lead()` and `diff()` in the formulas of `panglm()`, `paniv()`,
  `panqr()` and `panfmm()` are taken within each individual, with `NA` at
  the start of each panel (plm). Model functions that cannot support them
  stop with an informative error.
* `cre = TRUE`: correlated random effects (Mundlak) for pooled and
  random-effects fits of every family, with a joint test in `summary()`.
* `dynamic = TRUE`: lagged and initial outcomes with Mundlak means
  (Wooldridge 2005 initial conditions).
* `partialeffects()`: average partial effects with delta-method standard
  errors.
* Weights: pooled and random-effects fits of every family, `twopart()`,
  `panfmm()`, `attgt()` and `etwfe()` accept weights, given as a vector or a
  column name. `offset` also accepts a column name.
* Splines: `jointtest()`, `lineartest()`, `effectcurve()` with a `plot()`
  method, and `model.matrix()` that keeps the knots on new data.
* `panhausman(method = "mundlak")`: the regression-based Hausman test.
* `pangmm()` accepts `vcov = "cluster"`.
* `tidy()` and `glance()` for every model class, so modelsummary works;
  `glance()` reports theta for negative binomial fits and an adjusted
  R-squared that counts absorbed fixed effects (fixest).
* `predict()` for `pangmm`, `paniv` and `panselect`; `summary()`, `vcov()`,
  `nobs()`, `confint()` and `logLik()` for event studies; `confint()` for
  `twopart`; full methods for `panhurdle()` fits; one print style across
  classes.
* `options(panglm.threads)` caps the parallel threads, which are also capped
  by `OMP_THREAD_LIMIT` and at 2 under `R CMD check`. OpenMP is no longer
  enabled.

## Robustness

* Maximum likelihood fits standardize the design internally and warn when
  the Hessian is not negative definite, standard errors are missing or the
  gradient is not zero. Badly scaled regressors no longer break beta fits.
* The random-effects negative binomial warns when its dispersion diverges
  (no overdispersion beyond the individual effect, so the model reduces to
  the random-effects Poisson).
* Random-slope fits warn when the random-effects covariance is singular, for
  every family.
* The random-effects Tobit is stable in the far censored tail.
* The random-effects beta likelihood is about 3 times faster.

## Bug fixes

* Robust and cluster covariances: the pooled negative binomial used the
  Poisson score, and the two-way fixed-effects Poisson ignored the time
  effects; both now match sandwich and `fixest::fepois()`.
* A fit made with `vcov = "HC1"` or `"cluster"` overwrote its model-based
  covariance; `vcov(fit, type = "classical")` now returns it.
* `panhausman()` always uses the model-based covariances in the classical
  contrast.
* Negative binomial summaries print theta and its standard error instead of
  "Dispersion parameter: 1".
* `did()` and `eventstudy()` summaries show the user's call.
* `paniv(model = "random")` no longer warns about regressors it keeps.
* A factor individual index works.
* The `copd` help page documents all 10 columns.

## Dependencies and documentation

* `data.table` is no longer a dependency.
* Roxygen markdown is enabled and every help page was reviewed; functions
  are grouped into families with cross-links.
* New vignettes: a getting started guide, four task guides (outcomes,
  endogeneity, policy evaluation, heterogeneity) and a technical reference.
* A pkgdown site and a rewritten README.

# panglm 1.1.4

* Added `\value` documentation to `confint.panglm()`, `fitted.panglm()`,
  `summary.panglm()`, and `vcov.panglm()`, describing the structure and
  meaning of each return value (CRAN resubmission).

# panglm 1.1.3

* Added method references (Hausman 1978, Chamberlain 1980, Allison and
  Waterman 2002, Croissant and Millo 2008) to the `DESCRIPTION` as
  `<doi:...>` links.

# panglm 1.1.2

## Bug fixes

* `model = "random"`, `family = "gaussian"` always reported `logLik`/`AIC`/
  `BIC` as `NA`, even though `fit_random_gaussian()` already estimates the
  Swamy-Arora variance components (`sigma_v2`, `sigma_mu2`) needed to
  compute a real likelihood - they were simply discarded. Added the
  closed-form Gaussian random-intercept log-likelihood (Sherman-Morrison
  identity on `Omega_i = sigma_v2 * I + sigma_mu2 * J`), using the already-
  estimated coefficients and variance components. Cross-validated: on the
  bundled `copd` data, `logLik(random)` closely tracks `logLik(pooling)` on
  the same formula, as expected since random-effects GLS nests pooled OLS
  as `sigma_mu2 -> 0`.
* `model = "within"` (one-way and two-way) had the same `loglik = NA_real_`
  gap. By the Frisch-Waugh-Lovell theorem the demeaned-regression residuals
  are algebraically identical to the full dummy-variable-regression
  residuals, so the standard OLS ML Gaussian log-likelihood applies
  directly - added, using the `npar` (regressors + fixed-effect count + 1
  dispersion parameter) that `panglm_parameter_count()` already computed
  for these branches but previously had no likelihood to pair with.

## New features

* `glance.panglm()` now reports `r.squared`/`adj.r.squared` (computed from
  `fitted.values()`/`y`) for every model type; previously not reported at
  all.

## Bug fixes

* `panglm(model = "random", effect = "individual", family = "gaussian")`
  printed a raw `solve(): system is singular; attempting approx solution`
  from Armadillo whenever a covariate was time-invariant (e.g. baseline sex
  alongside a pre/post indicator) or the panel was balanced with a
  deterministic time indicator (e.g. any two-wave pre/post design) - both
  make the intermediate within/between regressions used to estimate the
  Swamy-Arora variance components exactly rank-deficient. The estimates
  were already correct (verified against `lme4::lmer`), but the warning
  bypassed R's condition system entirely. Switched those two intermediate
  regressions to `stats::lm.fit()`, which handles rank deficiency via
  pivoted QR the same way `lm()` does, with no near-singular solve and no
  raw stderr output.

# panglm 1.1.0

CRAN submission candidate. Adds a bundled synthetic dataset, `copd` (a
simulated COPD follow-up panel), and runnable `@examples` to every exported
function and documented method.

* Added `plot.panglm()` for coefficient, residual, and fitted-value views,
  and `plot.panglm_hurdle()` for the two hurdle components.
* Added zero-truncated NB2 positive-count models to `panglm_hurdle()` through
  `count_family = "negbin"`.
* Added HC1 and cluster-robust covariance estimates for one-way conditional
  binomial fixed-effects models.
* Reorganized the package vignette as a methodological account of model
  specification, estimation, inference, diagnostics, and scope.
* Completed fitted-value and prediction semantics across estimators.
  Random-effects response predictions integrate over the fitted mixing
  distribution, fixed-effects predictions include known panel effects, and
  exact conditional binomial fits expose conditional inclusion
  probabilities. Fitted values and residuals now preserve the input row
  order for unsorted panels.
* Corrected likelihood parameter accounting. `logLik()` now counts fitted
  dispersion, shape, variance, and absorbed fixed-effect parameters where
  they belong to the reported likelihood. Gaussian pooled likelihoods use
  the maximum-likelihood variance while retaining the unbiased residual
  variance for covariance estimation. `glance()` now reports AIC, BIC, and
  the likelihood degrees of freedom.
* Corrected sandwich inference for exact conditional binomial models. HC1
  now operates on independent stratum scores, and custom clusters are
  required to preserve complete conditioned strata. Custom cluster vectors
  are now correctly reordered when a panel was supplied out of index order.
* Revised the fixed-effects NB2 description to state that the
  dummy-variable estimator remains subject to finite-T incidental-parameter
  concerns in short panels.
* Hardened model-frame handling. Factor interactions are retained,
  unidentified columns are reported and removed before numerical fitting,
  `na.exclude` restores fitted values and residuals to the input length, and
  panel indexes are checked for missing or duplicate individual-time keys.
  Nontrivial weights and nonzero offsets now fail explicitly instead of
  risking an apparently successful fit that ignores them.

# panglm 0.5.0

First CRAN submission candidate. Consolidates the estimators, inference, and
zero-inflated-count support added since 0.1.0 (see below); no functional
changes since the last version bump beyond package-metadata and
documentation cleanup for submission (`URL`/`BugReports` fields, vignette
housekeeping).

# panglm 0.1.0

Initial release.

## Estimators

* `model = "pooling"`: gaussian/poisson/binomial GLM via IRLS, and genuine
  NB2 (negative binomial) via theta-profiled IRLS -- matches
  `MASS::glm.nb()` exactly. (`family = "negbin"` previously silently
  aliased to Poisson at every model level; this is fixed everywhere below.)
* `model = "within"`: gaussian (exact demeaning, one-way or `effect =
  "twoways"`), poisson (conditional MLE, one-way or `effect = "twoways"`
  via outer-IRLS/inner-weighted-FWL -- the algorithm behind
  `fixest::feglm()`), negative binomial (Allison & Waterman 2002
  unconditional dummy-variable estimator -- matches `fixest::fenegbin()`
  exactly; tractable for small-to-moderate N, not intended for very large
  panels), and binomial (exact conditional logistic regression, Chamberlain
  1980 -- matches `survival::clogit(method = "exact")`, no `pglm`
  equivalent). `effect = "twoways"` is available for gaussian, poisson, and
  negbin (the latter two via outer-IRLS/inner-weighted-FWL, matching
  `fixest::fenegbin()` exactly); not yet for binomial (no general
  closed-form two-way conditional logit exists).
* Groups with an all-zero outcome are now dropped (with an informative
  `message()`) before fitting the Allison-Waterman FE-NB2 estimator -- their
  dummy-variable intercept is unbounded under the log link, the same
  boundary case `model = "within", family = "binomial"` already screened
  for. Previously left in, this could contaminate the shared covariate
  coefficients via the joint solve, not just the offending group's own
  intercept.
* `model = "random"`: gaussian (Swamy-Arora), poisson (Poisson-Gamma, a
  single dispersion parameter), negative binomial (a genuinely different,
  2-parameter beta-negative-binomial mixture matching `pglm`'s
  `lnl.negbin.R` exactly -- previously aliased to the Poisson-Gamma model),
  and binomial (Gauss-Hermite quadrature) random effects. The
  random-effects Poisson fit now seeds from the pooled-Poisson MLE plus a
  method-of-moments dispersion estimate, and uses Levenberg-Marquardt/ridge
  damping (with a warning when it engages) instead of erroring on a
  singular Hessian -- verified on heavily zero-inflated, extremely skewed
  panels that previously crashed.

## Inference

* `vcov(fit, type = "HC1" | "cluster")` for pooled/within models, including
  `effect = "twoways"` (gaussian) and within negbin (Allison-Waterman): the
  sandwich correction uses the full covariate+dummy score and bread (all
  (k+G) parameters), not just the covariate block, since the dummy
  intercepts are jointly estimated and contribute sampling variability of
  their own. Random-effects models intentionally keep model-based
  (information-matrix) SEs only, matching `lme4`/`glmmTMB` convention.
* `confint.panglm()` -- Wald intervals using whichever vcov the fit
  currently carries.
* `panglm_dispersiontest()` -- Pearson chi-squared/df overdispersion test
  for poisson/negbin fits with available `fitted.values`. Observations with
  a zero or unavailable fitted value (e.g. rows in an all-zero group
  screened out of a fixed-effects negbin fit) are excluded from the
  statistic rather than producing a `0/0` = `NaN`; the excluded count is
  reported as `n_excluded` on the returned object.
* `panglm_hausman()` for FE-vs-RE specification testing.
* `tidy()`/`glance()` methods registered against `generics::tidy`/`glance`
  for `broom`/`modelsummary` compatibility.

## Zero-inflated counts

* `panglm_hurdle()`: a fixed-effects hurdle model for panels with
  structural zeros (a pattern that shows up as extreme, family-invariant
  overdispersion under either `family = "poisson"` or `family = "negbin"`
  in `panglm()`). Fits `1(y>0)` via the existing conditional-logit FE
  estimator, and a zero-truncated Poisson or NB2 FE model (an
  Allison-Waterman-style dummy-variable estimator, Newton-Raphson with
  exact Fisher information) on the `y>0` subsample. The truncated-Poisson
  math (no fixed effects) matches `pscl::hurdle()`'s count part exactly.
  Both zero-truncated likelihoods are validated against `pscl::hurdle()`.

## Known limitations

* `effect = "twoways"` not yet available for binomial (no general
  closed-form two-way conditional logit exists).
* No between effects.
* No ordinal/tobit families.
* Not yet on CRAN.
