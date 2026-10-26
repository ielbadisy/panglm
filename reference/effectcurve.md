# Estimated effect curve of a continuous regressor

For a regressor entered through a regression spline
([`splines::ns()`](https://rdrr.io/r/splines/ns.html),
[`splines::bs()`](https://rdrr.io/r/splines/bs.html)), a polynomial or
linearly, returns the change in the linear predictor when the regressor
moves from the reference value `ref` to each of `values`, with
delta-method pointwise confidence limits. Centering at a reference is
what makes the curve comparable across models: in fixed-effects fits the
level is absorbed by the individual effects. The scale is the link
scale: a difference in mean for the identity link, a log rate ratio for
Poisson and negative binomial models, a log odds ratio for the logit
link. Use [`plot()`](https://rdrr.io/r/graphics/plot.default.html) to
draw it.

## Usage

``` r
effectcurve(
  object,
  variable,
  values = NULL,
  ref = NULL,
  vcov = NULL,
  level = 0.95
)
```

## Arguments

- object:

  a `"panglm"` fit

- variable:

  name of the regressor

- values:

  values at which to evaluate the curve; default: 100 points between the
  1st and 99th percentile of the regressor in the fit

- ref:

  reference value; default: the median of the regressor

- vcov:

  `NULL` (the covariance stored on the fit) or a type passed to
  [`vcov.panglm()`](https://ielbadisy.github.io/panglm/reference/vcov.panglm.md)

- level:

  confidence level of the pointwise band

## Value

a data frame of class `"panglm_effectcurve"` with columns `value`,
`estimate`, `se`, `lower` and `upper`

## See also

Other panglm methods:
[`confint.panglm()`](https://ielbadisy.github.io/panglm/reference/confint.panglm.md),
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
fit <- panglm(fev1 ~ splines::ns(crp, 3), data = copd,
              index = c("id", "visit"), model = "within", vcov = "cluster")
curve_crp <- effectcurve(fit, "crp")
head(curve_crp)
#>       value  estimate         se      lower     upper
#> 1 0.1000000 0.1765138 0.05858363 0.06169198 0.2913356
#> 2 0.1989899 0.1724357 0.05652094 0.06165674 0.2832148
#> 3 0.2979798 0.1683591 0.05447000 0.06159989 0.2751184
#> 4 0.3969697 0.1642853 0.05243396 0.06151665 0.2670540
#> 5 0.4959596 0.1602158 0.05041601 0.06140218 0.2590293
#> 6 0.5949495 0.1561518 0.04841943 0.06125152 0.2510522
plot(curve_crp)
```
