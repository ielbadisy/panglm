# Confidence intervals for panglm coefficients

Wald-type intervals `estimate +/- z * se`, using whichever
[`vcov()`](https://rdrr.io/r/stats/vcov.html) the object currently
carries (classical by default; HC1/cluster if the fit or a prior
[`vcov.panglm()`](https://ielbadisy.github.io/panglm/reference/vcov.panglm.md)
call set it).

## Usage

``` r
# S3 method for class 'panglm'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  a `"panglm"` fit

- parm:

  which parameters (names or indices); defaults to all

- level:

  confidence level

- ...:

  unused

## Value

a numeric matrix with one row per requested parameter and two columns
giving the lower and upper confidence limits, named after the quantiles
implied by `level` (e.g. `"2.5 %"` and `"97.5 %"`). Row names are the
parameter names.

## See also

Other panglm methods:
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
[`tidy.panglm()`](https://ielbadisy.github.io/panglm/reference/tidy.panglm.md),
[`vcov.panglm()`](https://ielbadisy.github.io/panglm/reference/vcov.panglm.md)

## Examples

``` r
data(copd)
fit <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
              model = "within", family = "gaussian")
confint(fit)
#>           2.5 %      97.5 %
#> crp -0.04701532 -0.02359432
```
