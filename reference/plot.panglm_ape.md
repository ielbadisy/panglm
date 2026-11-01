# Plot average partial effects

Dot-and-whisker plot of the average partial effects returned by
[`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md):
one row per continuous or binary regressor and per non-reference factor
level, the estimate as a dot, its Wald confidence interval as a
horizontal segment and a dotted reference line at zero. The left margin
is widened to fit the row labels and restored on exit. For two-part
models the part of the mean (`"total"`, `"binary"` or `"positive"`) is
shown in the axis label; plot one object per part to compare them.

## Usage

``` r
# S3 method for class 'panglm_ape'
plot(x, labels = NULL, xlab = NULL, main = NULL, col = "black", pch = 19, ...)
```

## Arguments

- x:

  a `"panglm_ape"` object from
  [`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md)

- labels:

  row labels; by default the variable name for continuous regressors and
  `"variable: contrast"` for discrete ones

- xlab:

  horizontal-axis label

- main:

  plot title

- col:

  colour of the dots and whiskers

- pch:

  plotting symbol of the estimates

- ...:

  passed to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html)

## Value

`x`, invisibly

## See also

Other panglm methods:
[`coefplot()`](https://elbadisyimad.com/panglm/reference/coefplot.md),
[`confint.panglm()`](https://elbadisyimad.com/panglm/reference/confint.panglm.md),
[`effectcurve()`](https://elbadisyimad.com/panglm/reference/effectcurve.md),
[`fitted.panglm()`](https://elbadisyimad.com/panglm/reference/fitted.panglm.md),
[`fixef()`](https://elbadisyimad.com/panglm/reference/fixef.md),
[`glance.panglm()`](https://elbadisyimad.com/panglm/reference/glance.panglm.md),
[`logLik.panglm()`](https://elbadisyimad.com/panglm/reference/logLik.panglm.md),
[`model.matrix.panglm()`](https://elbadisyimad.com/panglm/reference/model.matrix.panglm.md),
[`partialeffects.panordered()`](https://elbadisyimad.com/panglm/reference/partialeffects.md),
[`plot.panglm()`](https://elbadisyimad.com/panglm/reference/plot.panglm.md),
[`plot.panglm_effectcurve()`](https://elbadisyimad.com/panglm/reference/plot.panglm_effectcurve.md),
[`predict.panglm()`](https://elbadisyimad.com/panglm/reference/predict.panglm.md),
[`predict.panordered()`](https://elbadisyimad.com/panglm/reference/predict.panordered.md),
[`residuals.panglm()`](https://elbadisyimad.com/panglm/reference/residuals.panglm.md),
[`summary.panglm()`](https://elbadisyimad.com/panglm/reference/summary.panglm.md),
[`tidy.panglm()`](https://elbadisyimad.com/panglm/reference/tidy.panglm.md),
[`vcov.panglm()`](https://elbadisyimad.com/panglm/reference/vcov.panglm.md)

## Examples

``` r
data(copd)
fit <- panglm(hospitalized ~ treatment + age + crp + factor(smoker),
              data = copd, index = c("id", "visit"), family = "probit")
plot(partialeffects(fit))
```
