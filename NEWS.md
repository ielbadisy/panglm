# panglm 0.6.0

CRAN submission candidate. Adds a bundled synthetic dataset, `copd` (a
simulated COPD follow-up panel), and runnable `@examples` to every exported
function and documented method -- previously undocumented, which is a CRAN
submission blocker. No functional changes to the estimators.

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
  estimator, and a zero-truncated Poisson FE model (an
  Allison-Waterman-style dummy-variable estimator, Newton-Raphson with
  exact Fisher information) on the `y>0` subsample. The truncated-Poisson
  math (no fixed effects) matches `pscl::hurdle()`'s count part exactly.
  Truncated NB2 for the count part is not yet implemented (documented as a
  known limitation below, not silently unsupported).

## Known limitations

* `effect = "twoways"` not yet available for binomial (no general
  closed-form two-way conditional logit exists).
* `panglm_hurdle()`'s count part is zero-truncated Poisson only; a
  zero-truncated NB2 count part is not yet implemented.
* No between effects.
* No ordinal/tobit families.
* Not yet on CRAN.
