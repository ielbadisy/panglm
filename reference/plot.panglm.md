# Plot a fitted panglm model

Provides base-graphics coefficient and diagnostic plots. The coefficient
view shows Wald confidence intervals. Diagnostic views are available
when the fitted object stores fitted values.

## Usage

``` r
# S3 method for class 'panglm'
plot(
  x,
  which = c("coefficients", "residuals", "fitted"),
  level = 0.95,
  main = NULL,
  xlab = NULL,
  ...
)
```

## Arguments

- x:

  a fitted `"panglm"` object

- which:

  one of `"coefficients"`, `"residuals"`, or `"fitted"`

- level:

  confidence level for coefficient intervals

- main:

  optional plot title

- xlab:

  optional horizontal-axis label

- ...:

  additional arguments passed to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html)

## Value

`x`, invisibly

## See also

Other panglm methods:
[`confint.panglm()`](https://ielbadisy.github.io/panglm/reference/confint.panglm.md),
[`effectcurve()`](https://ielbadisy.github.io/panglm/reference/effectcurve.md),
[`fitted.panglm()`](https://ielbadisy.github.io/panglm/reference/fitted.panglm.md),
[`glance.panglm()`](https://ielbadisy.github.io/panglm/reference/glance.panglm.md),
[`logLik.panglm()`](https://ielbadisy.github.io/panglm/reference/logLik.panglm.md),
[`model.matrix.panglm()`](https://ielbadisy.github.io/panglm/reference/model.matrix.panglm.md),
[`partialeffects.panordered()`](https://ielbadisy.github.io/panglm/reference/partialeffects.md),
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
plot(fit)
```
