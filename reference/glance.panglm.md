# One-row model summary of a panglm fit

Implements the
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
generic.

## Usage

``` r
# S3 method for class 'panglm'
glance(x, ...)
```

## Arguments

- x:

  a `"panglm"` fit

- ...:

  unused

## Value

a one-row `data.frame` including log-likelihood, AIC, BIC, the
likelihood parameter count when a likelihood is available, R-squared and
adjusted R-squared (see Details), and `theta`, the NB2 dispersion
parameter of pooled, fixed-effects and random-slope negative binomial
fits (`NA` otherwise)

## Details

`r.squared` is `1 - RSS / TSS` on the scale of
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and the
outcome. For Gaussian fixed-effects fits the fitted values include the
absorbed individual (and, with `effect = "twoways"`, time) effects, so
`r.squared` is the R-squared of the dummy-variable regression,
`fixest::r2(type = "r2")`, not the within R-squared reported by `plm`.
`adj.r.squared` is `1 - (1 - R^2) (n - 1) / df`, with `df` the residual
degrees of freedom of the fit, which count the absorbed effects, as
`fixest::r2(type = "ar2")`. For non-Gaussian families both are
Gaussian-style pseudo-R-squared measures, not likelihood based.

## See also

[tidy-methods](https://ielbadisy.github.io/panglm/reference/tidy-methods.md)
for the other model classes.

Other panglm methods:
[`confint.panglm()`](https://ielbadisy.github.io/panglm/reference/confint.panglm.md),
[`effectcurve()`](https://ielbadisy.github.io/panglm/reference/effectcurve.md),
[`fitted.panglm()`](https://ielbadisy.github.io/panglm/reference/fitted.panglm.md),
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
fit <- panglm(fev1 ~ treatment + age + smoker + crp, data = copd,
              index = c("id", "visit"), model = "pooling", family = "gaussian")
if (requireNamespace("generics", quietly = TRUE)) {
  generics::glance(fit)
}
#>     model     effect   family     link vcov.type    logLik      AIC      BIC
#> 1 pooling individual gaussian identity classical -189.5639 391.1278 413.3505
#>   df.logLik df.residual r.squared adj.r.squared nobs n.groups theta
#> 1         6         295 0.2417961     0.2315154  300       60    NA
```
