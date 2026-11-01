# Extract residuals from a panglm model

Response residuals `y - mu` or Pearson residuals
`(y - mu) / sqrt(V(mu))`, with `mu` the stored fitted values and `V` the
variance function of the family, unscaled by the dispersion as in
[`stats::residuals.glm()`](https://rdrr.io/r/stats/glm.summaries.html):
`1` (Gaussian), `mu` (Poisson), `mu + mu^2 / theta` (NB2 negative
binomial), `mu (1 - mu)` (binomial), `mu^2` (gamma) and
`mu (1 - mu) / (1 + phi)` (beta). For random-effects Poisson-gamma and
beta-negative-binomial random-intercept fits, whose fitted values are
marginal means, `V` is the marginal variance implied by the mixing
distribution; where it is not finite (beta parameter `a <= 2`) the
Poisson variance is used. For random-slope fits `V` is the conditional
variance function evaluated at the marginal mean. Pearson residuals are
not defined for the ordered and Tobit families.

## Usage

``` r
# S3 method for class 'panglm'
residuals(object, type = c("response", "pearson"), ...)
```

## Arguments

- object:

  a fitted `"panglm"` object

- type:

  `"response"` (default) or `"pearson"`

- ...:

  unused

## Value

a numeric vector, one value per row of the original data (in that row
order), `NA` at rows dropped for missingness.

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
[`plot.panglm_effectcurve()`](https://elbadisyimad.com/panglm/reference/plot.panglm_effectcurve.md),
[`predict.panglm()`](https://elbadisyimad.com/panglm/reference/predict.panglm.md),
[`predict.panordered()`](https://elbadisyimad.com/panglm/reference/predict.panordered.md),
[`summary.panglm()`](https://elbadisyimad.com/panglm/reference/summary.panglm.md),
[`tidy.panglm()`](https://elbadisyimad.com/panglm/reference/tidy.panglm.md),
[`vcov.panglm()`](https://elbadisyimad.com/panglm/reference/vcov.panglm.md)

## Examples

``` r
data(copd)
fit <- panglm(exacerbations ~ crp, data = copd,
              index = c("id", "visit"), model = "pooling",
              family = "poisson")
head(residuals(fit))
#>          1          2          3          4          5          6 
#> -0.1661605  1.0956670  2.9187342 -1.1280102  0.9553894 -2.0858936 
head(residuals(fit, type = "pearson"))
#>          1          2          3          4          5          6 
#> -0.1128970  0.7939755  2.0231625 -0.7732611  0.6681517 -1.4442623 
```
