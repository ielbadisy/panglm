# Plot an effect curve

Draws the estimate as a line over a shaded pointwise confidence band, a
rug of the percentiles of the observed regressor, and dotted reference
lines: at zero and at the reference value for link-scale curves, at zero
for slope curves. Curves by group share one panel, one colour each, with
a legend. Curves of a discrete regressor are drawn as points with
confidence whiskers.

## Usage

``` r
# S3 method for class 'panglm_effectcurve'
plot(
  x,
  xlab = attr(x, "variable"),
  ylab = NULL,
  band_col = "grey85",
  rug = TRUE,
  col = NULL,
  legend = "auto",
  ...
)
```

## Arguments

- x:

  a `"panglm_effectcurve"` object from
  [`effectcurve()`](https://elbadisyimad.com/panglm/reference/effectcurve.md)

- xlab, ylab:

  axis labels; `ylab` defaults to a label for the scale

- band_col:

  colour of the confidence band (single curve)

- rug:

  if `TRUE`, mark the percentiles of the observed regressor on the
  horizontal axis

- col:

  colour(s) of the curve(s); black for a single curve and the Okabe-Ito
  palette for curves by group by default

- legend:

  position of the legend for curves by group: `"auto"` (the corner of
  the plot least covered by the curves and bands), a keyword accepted by
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html), or
  `FALSE` for none

- ...:

  passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html)

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
[`plot.panglm_ape()`](https://elbadisyimad.com/panglm/reference/plot.panglm_ape.md),
[`predict.panglm()`](https://elbadisyimad.com/panglm/reference/predict.panglm.md),
[`predict.panordered()`](https://elbadisyimad.com/panglm/reference/predict.panordered.md),
[`residuals.panglm()`](https://elbadisyimad.com/panglm/reference/residuals.panglm.md),
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
