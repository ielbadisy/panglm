# Fit a generalized linear model for panel data

Fast pooled, fixed-effects ("within") and random-effects panel GLMs,
backed by an Rcpp/RcppArmadillo/RcppParallel numerical core. See the
package vignette for the estimators currently implemented per family.

## Usage

``` r
panglm(
  formula,
  data,
  index,
  model = c("pooling", "within", "random", "between"),
  family = "gaussian",
  weights = NULL,
  offset = NULL,
  na.action = stats::na.omit,
  effect = c("individual", "twoways"),
  vcov = c("classical", "HC1", "cluster"),
  cre = dynamic,
  dynamic = FALSE,
  random = ~1,
  R = 21,
  maxit = 100,
  tol = 1e-10
)
```

## Arguments

- formula:

  a model formula, e.g. `y ~ x1 + x2`. `lag(x, k)`, `lead(x, k)` and
  `diff(x)` are taken within each individual on the grid of sorted
  distinct time values (see "Lags and differences" in `panglm()`; needs
  `index = c(id, time)`); the first waves of each individual get `NA`
  and leave the sample, and `lag(x, 0:2)` gives several lags

- data:

  a data.frame containing the variables in `formula` and `index`

- index:

  length-1 or length-2 character vector identifying the individual (and
  optionally time) columns, e.g. `c("id", "time")`

- model:

  one of `"pooling"`, `"within"`, `"random"`, `"between"`. The between
  estimator (`family = "gaussian"` only) is least squares on the
  individual means of the outcome and the regressors, one row per
  individual, as in `plm(model = "between")`; it uses only the variation
  between individuals and is the cross-sectional counterpart of the
  within estimator. Its robust covariances (`"HC1"`, `"cluster"`) are
  the same, since each individual is one observation of the means
  regression.

- family:

  one of `"gaussian"`, `"poisson"`, `"binomial"`, `"probit"`,
  `"negbin"`, `"gamma"`, `"ologit"`, `"oprobit"`, `"beta"`, `"tobit"`,
  or a family spec from
  [`gaussianfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  /
  [`poissonfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  /
  [`binomialfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  /
  [`negbinfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  /
  [`gammafamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  /
  [`orderedfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  /
  [`betafamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  /
  [`tobitfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)

- weights:

  optional nonnegative observation (prior) weights: the name of a column
  of `data` or a numeric vector with one value per row of `data`, e.g.
  inverse probability weights for attrition. Nontrivial weights are
  supported for `model = "pooling"` (all families); the fit then
  maximizes the weighted (pseudo) log-likelihood and robust covariances
  use the weighted scores. For `model = "random"` the weights must be
  constant within individuals (one weight per individual, e.g. a
  baseline sampling weight) and the fit maximizes the weighted sum of
  the individual integrated log-likelihoods; the Gaussian model is then
  estimated by maximum likelihood instead of Swamy-Arora FGLS. Use
  `vcov = "cluster"` for design-based standard errors. Weighted
  fixed-effects fits are rejected.

- offset:

  optional model offset (the name of a column of `data` or a numeric
  vector with one value per row of `data`). Nonzero offsets, including
  offsets specified in `formula`, are rejected explicitly rather than
  ignored.

- na.action:

  missing-data action passed to
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)

- effect:

  one of `"individual"` or `"twoways"` (individual + time fixed
  effects). `"twoways"` is implemented for `model = "within"` only:
  `family = "gaussian"` (exact alternating-projections demeaning),
  `"poisson"` and `"negbin"` (both via outer-IRLS / inner-weighted-FWL,
  the algorithm behind
  [`fixest::feglm()`](https://lrberge.github.io/fixest/reference/feglm.html)),
  and `"binomial"`, the two-way conditional (Chamberlain) logit: the
  individual effects are conditioned out exactly and the time effects
  are estimated as time dummies inside the conditional likelihood, as in
  `survival::clogit(y ~ x + factor(time) + strata(id), method = "exact")`.
  For the two-way logit [`coef()`](https://rdrr.io/r/stats/coef.html)
  holds the slopes only; the time effects (first time level as
  reference, `NA` when not identified) and their standard errors are in
  `fit$time_effects` and `fit$time_effects_se`, and the classical and
  robust covariances of the slopes account for their estimation.

- vcov:

  one of `"classical"`, `"HC1"` (heteroskedasticity-robust) or
  `"cluster"` (cluster-robust, clustered by the panel individual by
  default). For `model = "random"` only `"classical"` and `"cluster"`
  are available; the cluster covariance uses the individual scores of
  the integrated likelihood, so clusters must contain whole individuals.
  See
  [`vcov.panglm()`](https://elbadisyimad.com/panglm/reference/vcov.panglm.md)
  to recompute a different type after fitting.

- cre:

  logical. If `TRUE`, fit a correlated random effects (Mundlak) model:
  the individual means of every time-varying regressor are added to the
  design matrix, named `mean_<column>`. Available for
  `model = "pooling"` and `model = "random"` and every family. The
  slopes of the time-varying regressors then have a within
  interpretation, and the joint test of the mean terms is a robust
  Hausman-type test (reported by
  [`summary.panglm()`](https://elbadisyimad.com/panglm/reference/summary.panglm.md)).
  Defaults to `dynamic`.

- dynamic:

  logical. If `TRUE`, fit a dynamic model with the Wooldridge (2005)
  solution to the initial conditions problem: the outcome of the
  previous wave (`lag_<y>`) and the outcome of the first observed wave
  (`initial_<y>`) are added to the regressors, the first wave is
  dropped, and the individual means of the time-varying regressors are
  added (`cre = TRUE` by default; the lag is not averaged). The second
  element of `index` must be a numeric time variable; a lag is only
  defined when the previous wave `t - 1` is observed, other rows are
  dropped. Factor and logical outcomes (ordered and binary models) give
  factor or logical lags, so their partial effects are discrete changes.
  With `model = "random"` the coefficient of the lag measures true state
  dependence net of unobserved heterogeneity;
  [`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md)
  gives it on the probability or mean scale. Available for
  `model = "pooling"` and `model = "random"`.

- random:

  one-sided formula for the random effects of `model = "random"`. The
  default `~ 1` is the random intercept. A formula with variables, e.g.
  `~ time`, adds correlated random slopes:
  `g(E[y_it | b_i]) = x_it' beta + z_it' b_i` with `z_it = (1, time)`
  and `b_i ~ N(0, Sigma)`, `Sigma` unstructured (Cholesky
  parameterized), the model of `lme4::lmer(y ~ x + (1 + time | id))`.
  The slope variables usually also appear in `formula` (their average
  slope). Available for `family = "gaussian"` (exact maximum likelihood,
  any number of slopes), and for `"poisson"`, `"negbin"` (NB2, the model
  of
  [`glmmTMB::nbinom2()`](https://rdrr.io/pkg/glmmTMB/man/nbinom2.html),
  with `theta` stored as `negbin_theta`), `"gamma"` (log link, the model
  of `glmmTMB` with `Gamma(link = "log")`, shape stored as `shape`) and
  `"binomial"` (logit or probit), with any number of slopes, by adaptive
  Gauss-Hermite product quadrature with `R` nodes per dimension, `R^q`
  nodes for `q` random effects; `R = 1` is the Laplace approximation of
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) and
  `glmmTMB`. With random slopes the random effects are always normal on
  the link scale, so for `"poisson"` and `"negbin"` the multiplicative
  effect `exp(z_it' b_i)` is lognormal. These are different models from
  the random-intercept fits of `random = ~ 1`: the random-intercept
  `"poisson"` model uses gamma mixing (the conjugate model of Hausman,
  Hall and Griliches, whose likelihood is in closed form, as Stata's
  `xtpoisson, re`) and the random-intercept `"negbin"` model is the
  beta-mixed negative binomial of Hausman, Hall and Griliches (Stata's
  `xtnbreg, re`), whereas with random slopes `"negbin"` is NB2 with a
  lognormal effect, as
  [`glmmTMB::nbinom2()`](https://rdrr.io/pkg/glmmTMB/man/nbinom2.html).
  Adding a slope whose variance is estimated at zero therefore does not
  reproduce the `random = ~ 1` fit for these two families: the slope
  coefficients keep their conditional (rate ratio) meaning, but the
  intercept (`E[exp(b_i)] = exp(sigma^2 / 2)` under lognormal mixing
  against 1 under gamma mixing), the variance parameters and the
  log-likelihood differ. For `"binomial"` and `"gamma"` the
  random-intercept models are already normal on the link scale and are
  the zero-slope special case. The fit stores `varcomp` (the covariance
  matrix `Sigma`), `ranef` (posterior modes, the BLUPs in the Gaussian
  model, one row per individual) and, for the Gaussian model, the
  residual variance `sigma_v2`. Weights, `dynamic = TRUE` and the other
  families are not supported with random slopes; `cre = TRUE` is.
  Response predictions and
  [`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md)
  are population averages over the random effects.

- R:

  number of Gauss-Hermite quadrature nodes (random-effects binomial and
  gamma models only). With random slopes, the number of nodes per
  dimension of the adaptive product rule. When `R` is not given it is 7
  for one slope, 5 for two, 3 for three to five and 1 (Laplace) beyond,
  so that `R^q` stays moderate; a warning is given when `R^q`
  exceeds 2000. `R = 1` gives the Laplace approximation.

- maxit:

  maximum IRLS/Newton iterations

- tol:

  convergence tolerance

## Value

an object of class `"panglm"` (with the additional class `"panordered"`
for the ordered families): a list with `coefficients`, `vcov`,
`vcov_type`, `fitted.values`, `linear.predictors`, `loglik`, the
variance components of random-effects fits (`sigma_mu2`, `sigma_v2`,
`theta`), the panel structure, `nobs`, `n_groups`, `call` and `formula`.
Methods: `print`, `summary`, `coef`, `vcov`, `confint`, `predict`,
`fitted`, `residuals`, `logLik`, `nobs`, `model.matrix`, `plot`,
[`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md),
[`effectcurve()`](https://elbadisyimad.com/panglm/reference/effectcurve.md),
`tidy` and `glance`.

## Lags and differences

Inside the formula, `lag(x, k)`, `lead(x, k)` and
`diff(x, lag, differences)` are taken within each individual, never
across individuals. Shifts step through the sorted distinct values of
the time index, so `lag(x, 1)` is the value at the previous observed
wave. With gaps in the time index this differs from plm: for waves coded
2010, 2012, 2014, `lag(x, 1)` at 2012 is the 2010 value in panglm, while
[`plm::lag()`](https://rdrr.io/r/stats/lag.html) (which shifts by the
numeric time value) gives `NA`. To obtain plm's behaviour, keep a row
for every period, with `NA` outcomes for the missing waves.
[`pangmm()`](https://elbadisyimad.com/panglm/reference/pangmm.md) uses
the same rule. The first `k` waves of each individual get `NA` and leave
the estimation sample. [`diff()`](https://rdrr.io/r/base/diff.html)
follows base R's `diff(x, lag, differences)`; plm's third argument of
[`diff()`](https://rdrr.io/r/base/diff.html) is `shift`, not
`differences`. The same rule applies in
[`paniv()`](https://elbadisyimad.com/panglm/reference/paniv.md),
[`panqr()`](https://elbadisyimad.com/panglm/reference/panqr.md) and
[`panfmm()`](https://elbadisyimad.com/panglm/reference/panfmm.md).

## Validation

Each estimator is checked against an established implementation of the
same model, on the same data, in the package's validation suite on
GitHub and in the package tests:

- pooled models: [`stats::lm()`](https://rdrr.io/r/stats/lm.html),
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html),
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html),
  [`plm::plm()`](https://rdrr.io/pkg/plm/man/plm.html) and fixest, with
  `sandwich` robust and cluster covariances;

- fixed effects: [`plm::plm()`](https://rdrr.io/pkg/plm/man/plm.html),
  [`fixest::feols()`](https://lrberge.github.io/fixest/reference/feols.html),
  [`fixest::fepois()`](https://lrberge.github.io/fixest/reference/feglm.html),
  [`fixest::fenegbin()`](https://lrberge.github.io/fixest/reference/femlm.html),
  [`pglm::pglm()`](https://rdrr.io/pkg/pglm/man/pglm.html),
  [`glm()`](https://rdrr.io/r/stats/glm.html) with individual dummies,
  and
  [`survival::clogit()`](https://rdrr.io/pkg/survival/man/clogit.html)
  for the conditional and two-way logit;

- between estimator: `plm::plm(model = "between")`;

- random effects: [`plm::plm()`](https://rdrr.io/pkg/plm/man/plm.html)
  and [`pglm::pglm()`](https://rdrr.io/pkg/pglm/man/pglm.html) (random
  intercept), [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html),
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) and glmmTMB
  (random intercept and slopes), and GLMMadaptive (adaptive quadrature
  with random slopes);

- ordered outcomes:
  [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html),
  [`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html),
  [`ordinal::clmm()`](https://rdrr.io/pkg/ordinal/man/clmm.html), and a
  blow-up and cluster fit with
  [`survival::clogit()`](https://rdrr.io/pkg/survival/man/clogit.html)
  for the fixed-effects ordered logit;

- Tobit: [`AER::tobit()`](https://rdrr.io/pkg/AER/man/tobit.html) and
  censReg; beta:
  [`betareg::betareg()`](https://rdrr.io/pkg/betareg/man/betareg.html)
  and glmmTMB; gamma: [`stats::glm()`](https://rdrr.io/r/stats/glm.html)
  and glmmTMB.

## See also

[panglm-families](https://elbadisyimad.com/panglm/reference/panglm-families.md)
for the outcome families,
[`summary.panglm()`](https://elbadisyimad.com/panglm/reference/summary.panglm.md),
[`vcov.panglm()`](https://elbadisyimad.com/panglm/reference/vcov.panglm.md),
[`predict.panglm()`](https://elbadisyimad.com/panglm/reference/predict.panglm.md)
and
[`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md)
for the fitted model,
[tidy-methods](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
for `tidy()` and `glance()`.

Other panel models:
[`panfmm()`](https://elbadisyimad.com/panglm/reference/panfmm.md),
[`panfmmselect()`](https://elbadisyimad.com/panglm/reference/panfmmselect.md),
[`pangmm()`](https://elbadisyimad.com/panglm/reference/pangmm.md),
[`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md),
[`paniv()`](https://elbadisyimad.com/panglm/reference/paniv.md),
[`panqr()`](https://elbadisyimad.com/panglm/reference/panqr.md),
[`panselect()`](https://elbadisyimad.com/panglm/reference/panselect.md),
[`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md)

## Examples

``` r
data(copd)

# pooled OLS on lung function (ignores the panel structure)
fit_pool <- panglm(fev1 ~ treatment + age + smoker + crp, data = copd,
                    index = c("id", "visit"), model = "pooling", family = "gaussian")
summary(fit_pool)
#> 
#> Call:
#> panglm(formula = fev1 ~ treatment + age + smoker + crp, data = copd, 
#>     index = c("id", "visit"), model = "pooling", family = "gaussian")
#> 
#> Family: gaussian  Link: identity  Model: pooling  Effect: individual 
#> N = 300  groups = 60  vcov: classical 
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  4.091479   0.225597  18.136  < 2e-16 ***
#> treatment    0.208556   0.056578   3.686 0.000228 ***
#> age         -0.014239   0.003133  -4.545  5.5e-06 ***
#> smoker      -0.355207   0.055905  -6.354  2.1e-10 ***
#> crp         -0.032895   0.013135  -2.504 0.012269 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -189.6 
#> Dispersion parameter: 0.2107 
#> Residual df: 295  Iterations: 2
#> 

# fixed-effects ("within") estimator: only the time-varying regressor
# (crp) survives demeaning; time-invariant covariates like treatment
# or age are not identified by a fixed-effects estimator
fit_fe <- panglm(fev1 ~ crp, data = copd,
                  index = c("id", "visit"), model = "within", family = "gaussian")
coef(fit_fe)
#>         crp 
#> -0.03530482 

# random-effects (Swamy-Arora) estimator
fit_re <- panglm(fev1 ~ treatment + age + smoker + crp, data = copd,
                  index = c("id", "visit"), model = "random", family = "gaussian")
coef(fit_re)
#> (Intercept)   treatment         age      smoker         crp 
#>  4.10943890  0.20549733 -0.01430869 -0.35583462 -0.03519775 

# fixed-effects Poisson on exacerbation counts
fit_pois <- panglm(exacerbations ~ crp, data = copd,
                    index = c("id", "visit"), model = "within", family = "poisson")
coef(fit_pois)
#>        crp 
#> 0.05112674 

# correlated random effects (Mundlak) probit for hospitalization
fit_cre <- panglm(hospitalized ~ treatment + age + crp, data = copd,
                  index = c("id", "visit"), model = "random",
                  family = "probit", cre = TRUE)
summary(fit_cre)
#> 
#> Call:
#> panglm(formula = hospitalized ~ treatment + age + crp, data = copd, 
#>     index = c("id", "visit"), model = "random", family = "probit", 
#>     cre = TRUE)
#> 
#> Family: binomial  Link: probit  Model: random  Effect: individual 
#> N = 300  groups = 60  vcov: classical 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -2.72956    1.07818  -2.532 0.011353 *  
#> treatment   -0.47330    0.22838  -2.072 0.038225 *  
#> age          0.03989    0.01126   3.541 0.000398 ***
#> crp          0.05887    0.04543   1.296 0.194953    
#> mean_crp    -0.12856    0.11532  -1.115 0.264921    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -159.8 
#> Dispersion parameter: 1 
#> Random-intercept SD: 0.2689 
#> Mundlak test (individual means jointly zero): chi2(1) = 1.243, p = 0.2649
#> Residual df: 294  Iterations: 43
#> 

# patient trajectories: random intercept and random slope of FEV1
# decline per visit, correlated
fit_rs <- panglm(fev1 ~ visit + treatment + age + smoker, data = copd,
                 index = c("id", "visit"), model = "random",
                 random = ~ visit)
#> Warning: the random-effects covariance matrix is singular (on the boundary of the parameter space): the variance of 'visit' given the preceding random effects is estimated at zero; consider a simpler 'random' formula
summary(fit_rs)
#> 
#> Call:
#> panglm(formula = fev1 ~ visit + treatment + age + smoker, data = copd, 
#>     index = c("id", "visit"), model = "random", random = ~visit)
#> 
#> Family: gaussian  Link: identity  Model: random  Effect: individual 
#> N = 300  groups = 60  vcov: classical 
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  3.872439   0.414554   9.341  < 2e-16 ***
#> visit       -0.004717   0.008262  -0.571  0.56808    
#> treatment    0.251411   0.110762   2.270  0.02322 *  
#> age         -0.013627   0.006396  -2.131  0.03312 *  
#> smoker      -0.342166   0.114817  -2.980  0.00288 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -38.65 
#> Dispersion parameter: 0.04087 
#> Random effects (SD): (Intercept) 0.42204, visit 0.00297 
#> Random-effects correlation: visit ~ (Intercept) -1 
#> Residual df: 295  Iterations: 19
#> 
fit_rs$varcomp
#>              (Intercept)         visit
#> (Intercept)  0.178115843 -1.253264e-03
#> visit       -0.001253264  8.818254e-06

# two-way conditional logit: smoking across survey waves, with year
# shocks (e.g. a tobacco tax) absorbed by the time effects
set.seed(1)
smk <- expand.grid(wave = 1:5, id = 1:200)
smk$price <- rnorm(nrow(smk)) + 0.2 * smk$wave
smk$smoking <- rbinom(nrow(smk), 1, plogis(rnorm(200)[smk$id] - 0.5 * smk$price +
                                               c(0, 0.2, -0.3, -0.6, -0.8)[smk$wave]))
fit_2w <- panglm(smoking ~ price, data = smk, index = c("id", "wave"),
                 model = "within", family = "binomial", effect = "twoways",
                 vcov = "cluster")
#> 37 group(s) with no within-group outcome variation dropped (all-0 or all-1) -- they carry no information for the conditional likelihood.
summary(fit_2w)
#> 
#> Call:
#> panglm(formula = smoking ~ price, data = smk, index = c("id", 
#>     "wave"), model = "within", family = "binomial", effect = "twoways", 
#>     vcov = "cluster")
#> 
#> Family: binomial  Link: logit  Model: within  Effect: twoways 
#> N = 1000  groups = 200  vcov: cluster 
#> 
#> Coefficients:
#>       Estimate Std. Error z value Pr(>|z|)    
#> price -0.57501    0.08186  -7.024 2.15e-12 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -282.3 
#> Dispersion parameter: 1 
#> Iterations: 34
#> 
fit_2w$time_effects
#>          1          2          3          4          5 
#>  0.0000000  0.3259111 -0.2093176 -0.4856443 -0.5288287 

# nonlinear effect through a natural spline: linearity test and curve
fit_ns <- panglm(fev1 ~ splines::ns(crp, 3), data = copd,
                 index = c("id", "visit"), model = "within", vcov = "cluster")
lineartest(fit_ns, "crp")
#> 
#> Joint Wald test (vcov: cluster)
#> H0: the effect of crp is linear 
#> Chi-squared = 1.338, df = 2, p-value = 0.5122
#> 
plot(effectcurve(fit_ns, "crp"))
```
