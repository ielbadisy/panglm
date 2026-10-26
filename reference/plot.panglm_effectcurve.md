# Plot an effect curve

Plot an effect curve

## Usage

``` r
# S3 method for class 'panglm_effectcurve'
plot(
  x,
  xlab = attr(x, "variable"),
  ylab = "Effect on the linear predictor",
  band_col = "grey85",
  ...
)
```

## Arguments

- x:

  a `"panglm_effectcurve"` object from
  [`effectcurve()`](https://elbadisyimad.com/panglm/reference/effectcurve.md)

- xlab, ylab:

  axis labels

- band_col:

  colour of the confidence band

- ...:

  passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html)

## Value

`x`, invisibly

## See also

Other panglm methods:
[`confint.panglm()`](https://elbadisyimad.com/panglm/reference/confint.panglm.md),
[`effectcurve()`](https://elbadisyimad.com/panglm/reference/effectcurve.md),
[`fitted.panglm()`](https://elbadisyimad.com/panglm/reference/fitted.panglm.md),
[`glance.panglm()`](https://elbadisyimad.com/panglm/reference/glance.panglm.md),
[`logLik.panglm()`](https://elbadisyimad.com/panglm/reference/logLik.panglm.md),
[`model.matrix.panglm()`](https://elbadisyimad.com/panglm/reference/model.matrix.panglm.md),
[`partialeffects.panordered()`](https://elbadisyimad.com/panglm/reference/partialeffects.md),
[`plot.panglm()`](https://elbadisyimad.com/panglm/reference/plot.panglm.md),
[`predict.panglm()`](https://elbadisyimad.com/panglm/reference/predict.panglm.md),
[`predict.panordered()`](https://elbadisyimad.com/panglm/reference/predict.panordered.md),
[`summary.panglm()`](https://elbadisyimad.com/panglm/reference/summary.panglm.md),
[`tidy.panglm()`](https://elbadisyimad.com/panglm/reference/tidy.panglm.md),
[`vcov.panglm()`](https://elbadisyimad.com/panglm/reference/vcov.panglm.md)

## Examples

``` r
data(copd)
fit <- panglm(fev1 ~ splines::ns(crp, 3), data = copd,
              index = c("id", "visit"), model = "within")
plot(effectcurve(fit, "crp"), band_col = "lightblue")
```
