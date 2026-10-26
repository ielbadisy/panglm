# Average partial effects

Average partial effects (APE, also called average marginal effects) of
each regressor on the expected outcome, with delta-method standard
errors. For a continuous regressor the APE is the sample average of the
derivative of the mean with respect to that regressor; for a binary (0/1
or logical) regressor and for each non-reference level of a factor it is
the sample average of the discrete change in the mean, all other
regressors held at their observed values. Derivatives and changes are
taken through the model formula, so interactions, polynomial terms and
factor codings are handled automatically.

## Usage

``` r
# S3 method for class 'panordered'
partialeffects(object, variables = NULL, level = 0.95, ...)

partialeffects(object, ...)

# S3 method for class 'panglm'
partialeffects(object, variables = NULL, level = 0.95, ...)

# S3 method for class 'twopart'
partialeffects(
  object,
  variables = NULL,
  level = 0.95,
  part = c("total", "binary", "positive"),
  ...
)
```

## Arguments

- object:

  a `"panglm"` or `"twopart"` fit

- variables:

  character vector of variable names (as they appear in the data);
  defaults to every regressor in the model formula(s)

- level:

  confidence level for the Wald intervals

- ...:

  unused

- part:

  for two-part models: `"total"` (the unconditional mean), `"binary"`
  (the probability of a positive outcome) or `"positive"` (the mean
  conditional on a positive outcome)

## Value

a data frame of class `"panglm_ape"` with columns `term`, `contrast`,
`estimate`, `std.error`, `statistic`, `p.value`, `conf.low` and
`conf.high`, one row per continuous or binary regressor and one row per
non-reference factor level.

## Details

The mean used is the one returned by `predict(type = "response")`: the
GLM mean for pooled models; the mean including the estimated individual
(and time) effects for fixed-effects models; and the population-averaged
mean (integrating over the random intercept, and the random slopes when
`random` has slope variables) for random-effects models, which is the
average structural function of Wooldridge (2005). For correlated random
effects fits (`cre = TRUE`) the individual means are held fixed, so the
APE of a time-varying regressor is the within-individual effect. For
two-part models the mean is `P(y > 0) * E[y | y > 0]` (or one of its
parts, see `part`).

Standard errors use the covariance stored on the fit (so a fit made with
`vcov = "cluster"` gives cluster-robust APE standard errors) and treat
variance components, dispersion and fixed-effect intercepts as fixed.
Exact conditional logit fits are not supported, because their individual
intercepts are conditioned out and the mean is not identified.

## Validation

Checked against
[`marginaleffects::avg_slopes()`](https://rdrr.io/pkg/marginaleffects/man/slopes.html)
(pooled probit) and analytic values (linear and fixed-effects Poisson)
in the package tests, and against hand computations for the Tobit and
the population-averaged random-slope models in the package's validation
suite on GitHub.

## References

Wooldridge, J. M. (2005). Unobserved heterogeneity and estimation of
average partial effects. In *Identification and Inference for
Econometric Models*, Cambridge University Press.

## See also

[`effectcurve()`](https://ielbadisy.github.io/panglm/reference/effectcurve.md)
for the effect of a spline term across its range.

Other panglm methods:
[`confint.panglm()`](https://ielbadisy.github.io/panglm/reference/confint.panglm.md),
[`effectcurve()`](https://ielbadisy.github.io/panglm/reference/effectcurve.md),
[`fitted.panglm()`](https://ielbadisy.github.io/panglm/reference/fitted.panglm.md),
[`glance.panglm()`](https://ielbadisy.github.io/panglm/reference/glance.panglm.md),
[`logLik.panglm()`](https://ielbadisy.github.io/panglm/reference/logLik.panglm.md),
[`model.matrix.panglm()`](https://ielbadisy.github.io/panglm/reference/model.matrix.panglm.md),
[`plot.panglm()`](https://ielbadisy.github.io/panglm/reference/plot.panglm.md),
[`plot.panglm_effectcurve()`](https://ielbadisy.github.io/panglm/reference/plot.panglm_effectcurve.md),
[`predict.panglm()`](https://ielbadisy.github.io/panglm/reference/predict.panglm.md),
[`predict.panordered()`](https://ielbadisy.github.io/panglm/reference/predict.panordered.md),
[`summary.panglm()`](https://ielbadisy.github.io/panglm/reference/summary.panglm.md),
[`tidy.panglm()`](https://ielbadisy.github.io/panglm/reference/tidy.panglm.md),
[`vcov.panglm()`](https://ielbadisy.github.io/panglm/reference/vcov.panglm.md)

## Examples

``` r
data(copd)
fit <- panglm(hospitalized ~ treatment + age + crp, data = copd,
              index = c("id", "visit"), model = "random",
              family = "probit", cre = TRUE)
partialeffects(fit)
#> 
#> Average partial effects
#> Delta-method standard errors, vcov: classical 
#> 
#>       term contrast estimate std.error statistic   p.value  conf.low conf.high
#>  treatment    1 - 0 -0.13560  0.063190    -2.146 0.0318600 -0.259400  -0.01176
#>        age    dY/dX  0.01161  0.003044     3.813 0.0001372  0.005640   0.01757
#>        crp    dY/dX  0.01713  0.013100     1.308 0.1909000 -0.008539   0.04280
#> 
```
