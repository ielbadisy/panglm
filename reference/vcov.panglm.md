# Recompute the coefficient covariance matrix for a panglm fit

With `type = NULL` (the default) the covariance the object was fit with
is returned (set by the `vcov` argument of
[`panglm()`](https://ielbadisy.github.io/panglm/reference/panglm.md)).
`type = "classical"` returns the model-based covariance, also for a fit
made with a robust `vcov`. `type = "HC1"` (heteroskedasticity-robust /
White) and `type = "cluster"` (cluster-robust sandwich, clustered by the
panel individual unless `cluster` is given) are available for pooled and
fixed-effects models. For random-effects models only `"cluster"` is
available: the units of the integrated likelihood are individuals, so
the sandwich is built from individual scores (clusters must contain
whole individuals) with the factor G / (G - 1). It guards against a
misspecified random-effects distribution and is the right choice with
sampling weights. For the unweighted Gaussian (Swamy-Arora) fit it is
the cluster sandwich of the quasi-demeaned regression, as in Stata's
`xtreg, re vce(cluster id)`. For exact conditional binomial models, the
independent likelihood units are panel strata rather than observations.
Their HC1 covariance therefore uses stratum-level scores. A custom
cluster may combine complete strata, but it may not split observations
from one conditioned stratum. For the two-way conditional logit
(`effect = "twoways"`) the time dummies are parameters of the
conditional likelihood: the sandwich is built from the full (slopes +
time dummies) stratum scores and bread, its small-sample factors count
all of them, and the slope block is returned.

## Usage

``` r
# S3 method for class 'panglm'
vcov(object, type = NULL, cluster = NULL, ...)
```

## Arguments

- object:

  a `"panglm"` fit

- type:

  `NULL` (the covariance stored on the fit) or one of `"classical"`,
  `"HC1"`, `"cluster"`

- cluster:

  optional vector giving the clustering variable (one value per original
  observation, i.e. `nrow(data)` long); defaults to the panel individual
  identifier used to fit the model

- ...:

  unused

## Value

a numeric covariance matrix for the estimated coefficients, with row and
column names matching `names(coef(object))`. For `type = "classical"`
this is the model-based covariance; for `"HC1"` and `"cluster"` it is a
sandwich (robust or cluster-robust) covariance estimate of the same
dimension.

## See also

Other panglm methods:
[`confint.panglm()`](https://ielbadisy.github.io/panglm/reference/confint.panglm.md),
[`effectcurve()`](https://ielbadisy.github.io/panglm/reference/effectcurve.md),
[`fitted.panglm()`](https://ielbadisy.github.io/panglm/reference/fitted.panglm.md),
[`glance.panglm()`](https://ielbadisy.github.io/panglm/reference/glance.panglm.md),
[`logLik.panglm()`](https://ielbadisy.github.io/panglm/reference/logLik.panglm.md),
[`model.matrix.panglm()`](https://ielbadisy.github.io/panglm/reference/model.matrix.panglm.md),
[`partialeffects.panordered()`](https://ielbadisy.github.io/panglm/reference/partialeffects.md),
[`plot.panglm()`](https://ielbadisy.github.io/panglm/reference/plot.panglm.md),
[`plot.panglm_effectcurve()`](https://ielbadisy.github.io/panglm/reference/plot.panglm_effectcurve.md),
[`predict.panglm()`](https://ielbadisy.github.io/panglm/reference/predict.panglm.md),
[`predict.panordered()`](https://ielbadisy.github.io/panglm/reference/predict.panordered.md),
[`summary.panglm()`](https://ielbadisy.github.io/panglm/reference/summary.panglm.md),
[`tidy.panglm()`](https://ielbadisy.github.io/panglm/reference/tidy.panglm.md)

## Examples

``` r
data(copd)
fit <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
              model = "within", family = "gaussian")
vcov(fit)
#>             crp
#> crp 3.56989e-05
vcov(fit, type = "cluster")
#>              crp
#> crp 3.686711e-05
```
