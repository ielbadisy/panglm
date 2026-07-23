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
  equivalent). `effect = "twoways"` for poisson/gaussian is parallelized in
  C++ (RcppParallel); not yet available for negbin/binomial.
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
  `effect = "twoways"` (gaussian); not yet for within negbin (the
  Allison-Waterman dummy-variable score isn't wired up for this yet --
  raises an informative error rather than a silently wrong answer).
  Random-effects models intentionally keep model-based (information-matrix)
  SEs only, matching `lme4`/`glmmTMB` convention.
* `confint.panglm()` -- Wald intervals using whichever vcov the fit
  currently carries.
* `panglm_dispersiontest()` -- Pearson chi-squared/df overdispersion test
  for poisson/negbin fits with available `fitted.values`.
* `panglm_hausman()` for FE-vs-RE specification testing.
* `tidy()`/`glance()` methods registered against `generics::tidy`/`glance`
  for `broom`/`modelsummary` compatibility.

## Known limitations

* `effect = "twoways"` not yet available for negbin/binomial.
* No between effects.
* No ordinal/tobit families.
* No cluster/HC1 vcov for within negbin.
* Not yet on CRAN.
