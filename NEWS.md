# panglm (development version)

* Between estimator: `model = "between"` (OLS on individual means) for the
  Gaussian family, validated against `plm(model = "between")`.
* Tobit: `tobit_family(left, right)` for outcomes censored at a floor and/or
  a ceiling, pooled or with a random intercept (adaptive Gauss-Hermite
  quadrature). `predict()` and `partialeffects()` work on the observed,
  censored-mean scale. Validated against `AER::tobit()` and `censReg`.
* Panel instrumental variables: `paniv()` for fixed-effects, random-effects
  (G2SLS or EC2SLS), pooled 2SLS and the Hausman-Taylor estimator, with
  first-stage F statistics, Sargan/Hansen tests and cluster covariances.
  Validated against `plm(inst)` and `plm::pht()`.
* Fixed-effects Poisson: the HC1 covariance (and clusters that cut across
  individuals) now use the correctly centred scores and match
  `fixest::fepois()`.
* Two-way fixed-effects logit: `panglm(model = "within", family = "binomial",
  effect = "twoways")` conditions out the individual effects and estimates
  time effects as dummies in the conditional likelihood, as in
  `survival::clogit(y ~ x + factor(time) + strata(id), method = "exact")`.
  `coef()` holds the slopes; `fit$time_effects` and `fit$time_effects_se`
  hold the time effects.
* Dynamic panel GMM: `pangmm()` for Arellano-Bond difference GMM and
  Blundell-Bond system GMM (one-step and two-step, robust and
  Windmeijer-corrected SEs), with AR(1)/AR(2) tests, Sargan and Hansen J
  tests, collapsed instruments (`collapse`) and lag-depth limits (`maxlag`).
  Uses the `plm::pgmm()` formula syntax and is validated against it.
* Staggered difference-in-differences: `attgt()` gives the Callaway and
  Sant'Anna (2021) group-time ATTs (never- or not-yet-treated controls;
  outcome regression, IPW or doubly robust; analytic influence-function
  SEs), aggregated by `summary(type = "simple", "dynamic", "group" or
  "calendar")` with a pre-trend Wald test and an event-study `plot()`.
  Validated against `did::att_gt()` and `did::aggte()`.
* Nonlinear difference-in-differences: `etwfe()` fits the Wooldridge (2023)
  extended two-way FE model for linear, Poisson and logit outcomes and
  reports ATTs on the response scale with delta-method SEs. Validated
  against `etwfe::etwfe()` and `etwfe::emfx()`.
* Panel sample selection: `panselect()` implements the Wooldridge (1995)
  correction for attrition or death in cohort studies: per-period probits
  with Mundlak means, inverse Mills ratio terms (by period or common) in a
  fixed-effects or pooled Mundlak outcome equation, panel bootstrap SEs by
  default and a Wald test of selection bias. Validated against `glm`, `plm`
  and, for one period, `sampleSelection::heckit()`.
* `pangmm()` warns when there are at least as many instruments as individuals
  (too many instruments bias GMM toward OLS and weaken the Hansen test) and
  says so in `summary()`; a singular weighting matrix now uses a generalized
  inverse, as `plm::pgmm()` does, instead of failing; the AR test returns `NA`
  instead of `NaN` when undefined. The help page documents the finite-sample
  bias of system GMM (Roodman 2009, Bun and Windmeijer 2010) and how to limit
  it (`collapse = TRUE`, smaller `maxlag`, FE/OLS bounds).
* New `panqr()`: panel quantile regression with the Canay (2011) two-step
  fixed-effects estimator, the Machado and Santos Silva (2019) MM-QR
  location-scale estimator and pooled quantile regression, for one or several
  `tau`, with bootstrap or analytic standard errors. Validated against
  `quantreg::rq()` and `plm`.
* `panglm(model = "random")` gains `random`, a one-sided formula for
  correlated random slopes (for example `random = ~ time`, unstructured
  covariance). Gaussian models use exact ML with any number of slopes,
  matching `lme4::lmer(REML = FALSE)`; Poisson, logit and probit models take
  one slope and use adaptive Gauss-Hermite quadrature (`R` nodes per
  dimension, default 7; `R = 1` is the Laplace approximation). Fits store
  `varcomp` and `ranef`; cluster-robust `vcov()`, `summary()`, `predict()`,
  `partialeffects()` and `cre = TRUE` are supported. Validated against lme4,
  nlme with clubSandwich, glmmTMB and GLMMadaptive.
* `validation/`: per-model scripts comparing panglm with established
  packages, summarized in `validation/REPORT.md`.

# panglm 1.4.0

* Ordered outcomes: `ordered_family("logit")` and `ordered_family("probit")`
  for pooled and random-intercept fits, and the fixed-effects ordered logit
  by "blow-up and cluster" (Baetschmann, Staub and Winkelmann 2015) with
  `model = "within"`. `predict()` returns category probabilities;
  `partialeffects()` gives one effect per category.
* Beta regression: `beta_family()` (logit or probit link) for outcomes in
  (0, 1), pooled or with a random intercept; precision in `fit$precision`.
* Dynamic models: `dynamic = TRUE` adds the lagged and initial outcomes and
  the Mundlak means (Wooldridge 2005 initial conditions), for pooled and
  random-effects fits of every family.
* Random-effects fits of every family accept weights constant within
  individuals and `vcov = "cluster"` (sandwich from the individual scores of
  the integrated likelihood). The weighted Gaussian model uses maximum
  likelihood.
* Regression splines: new `jointtest()` (joint Wald test of coefficients),
  `lineartest()` (Wald test that a spline effect is linear),
  `effectcurve()` with a `plot()` method (fitted curve with pointwise
  confidence band) and `model.matrix()` for `panglm` fits, which keeps the
  spline knots on new data.
* New `did()` (two-way fixed-effects difference-in-differences) and
  `eventstudy()` (event-time coefficients, joint pre-trend test, `plot()`
  method), for Gaussian, Poisson and negative binomial outcomes.
* `twopart(..., correlated = TRUE)`: random-intercept lognormal two-part
  model with correlated intercepts, fitted jointly (Van Ourti 2004), with
  population-averaged predictions and partial effects.
* Bug fix: robust (`HC1`, `cluster`) covariances of the two-way
  fixed-effects Poisson model ignored the time effects in the score; they
  now match `fixest::fepois()`. The same covariances are now available for
  the two-way negative binomial model.

# panglm 1.3.0

* Observation weights: pooled models accept nonnegative prior `weights`
  (for example inverse probability weights for attrition), for every
  family. Zero weights drop observations. Robust and cluster covariances use
  the weighted scores. Fixed- and random-effects models still reject
  nontrivial weights.
* New `gamma_family()` (`family = "gamma"`): gamma regression with log link
  for positive, skewed outcomes such as costs. Pooled (IRLS), fixed effects
  (one intercept per individual, consistent for fixed T with a
  multiplicative effect) and random intercept (adaptive Gauss-Hermite
  quadrature, `R = 1` is the Laplace approximation) estimators, with robust
  and cluster covariances for pooled and fixed-effects fits.
* New `cre` argument: correlated random effects (Mundlak) for pooled and
  random-effects fits of every family. Individual means of the
  time-varying regressors are added as `mean_<variable>` columns;
  `summary()` reports their joint Wald test. `predict()` recomputes the
  means on new data.
* New `partialeffects()`: average partial effects (continuous derivatives,
  binary and factor contrasts) with delta-method standard errors, for
  `panglm` and `twopart` fits. Random-effects means integrate over the
  random intercept.
* New `twopart()`: two-part model for semicontinuous costs, probit or logit
  first part times a gamma or lognormal positive part, pooled or random
  intercept, optional correlated random effects, with smearing or normal
  retransformation. Methods: `predict()`, `fitted()`, `coef()`, `vcov()`,
  `logLik()`, `nobs()`, `summary()`.
* Bug fix: the robust (`HC1`, `cluster`) covariance of the pooled negative
  binomial model used the Poisson score `(y - mu) x` instead of the NB2
  score `(y - mu) theta / (theta + mu) x`. It now matches
  `sandwich::vcovHC()` on `MASS::glm.nb()`.
* The pooled IRLS covariance is now evaluated at the final coefficients.

# panglm 1.2.0

* Bug fix: a fit made with `vcov = "HC1"` or `vcov = "cluster"` replaced
  its model-based covariance, so `vcov(fit, type = "classical")` returned
  the robust matrix. The model-based covariance is now kept on the fit.
  `vcov()` gains `type = NULL` (the new default), which returns the
  covariance the fit was made with, so `vcov(fit)` is unchanged.
* Bug fix: `panglm_hausman()` used the covariance stored on the FE fit, so
  an FE fit with a robust `vcov` entered the classical contrast
  `V_fe - V_re`, which is valid only for model-based covariances. The
  classical test now always uses the model-based covariances.
* `panglm_hausman()` gains `method = "mundlak"`: the regression-based
  (correlated random effects) version of the test, a Wald test of the
  individual-mean coefficients in the pooled model `y ~ x + xbar`, with a
  cluster-robust (default), HC1 or classical covariance. With the
  cluster-robust covariance the reference is `F(k, G - 1)`. Gaussian
  family, individual effects only.
* Test of `panglm_hausman()` against `plm::phtest()` now passes fitted
  `plm` models, since `phtest()`'s formula method calls `plm()` unqualified.

# panglm 1.1.5

* Dropped the `data.table` dependency entirely. `build_panel_index()`'s
  panel sort-by-(id, time) now uses plain base R data.frame indexing;
  its run-length grouping uses `base::rle()` instead of
  `data.table::rleidv()` (no `as.character()` coercion needed, since the
  id column is already atomic and pre-sorted). `demean_twoway()`'s dense
  time-rank encoding uses a new internal `dense_rank0()` helper
  (`match(x, sort(unique(x))) - 1L`) instead of `data.table::frankv()`.
  No basetable equivalent exists for either primitive, so this ended up
  as a full drop rather than a basetable swap. Also removed several
  dead `@importFrom data.table`/`utils::globalVariables()` declarations
  for NSE symbols (`.SD`, `..cols`, `..xnames`, `:=`) that were already
  unused before this change.

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
