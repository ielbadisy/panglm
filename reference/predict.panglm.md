# Predict from a fitted panglm model

Pooled predictions use the ordinary GLM mean. Random-effects response
predictions are population averages over the fitted latent-effect
distribution. Fixed-effects predictions include the estimated individual
and, where applicable, time effects. New fixed-effect levels produce an
error unless `allow.new.levels = TRUE`, in which case their contribution
is set to zero.

## Usage

``` r
# S3 method for class 'panglm'
predict(
  object,
  newdata = NULL,
  type = c("link", "response"),
  allow.new.levels = FALSE,
  ...
)
```

## Arguments

- object:

  a fitted `"panglm"` object

- newdata:

  optional data frame; when omitted, return in-sample values

- type:

  either `"link"` or `"response"`

- allow.new.levels:

  if `TRUE`, assign zero to unseen fixed-effect levels; otherwise unseen
  levels are an error

- ...:

  unused

## Value

a numeric vector of predictions

## Details

Exact conditional binomial fits provide in-sample conditional inclusion
probabilities. They cannot provide predictions for arbitrary `newdata`
because the individual intercept is conditioned out and the probability
depends on the complete response set and its success total.

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
[`plot.panglm_effectcurve()`](https://elbadisyimad.com/panglm/reference/plot.panglm_effectcurve.md),
[`predict.panordered()`](https://elbadisyimad.com/panglm/reference/predict.panordered.md),
[`summary.panglm()`](https://elbadisyimad.com/panglm/reference/summary.panglm.md),
[`tidy.panglm()`](https://elbadisyimad.com/panglm/reference/tidy.panglm.md),
[`vcov.panglm()`](https://elbadisyimad.com/panglm/reference/vcov.panglm.md)

## Examples

``` r
data(copd)
fit <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
              model = "random", family = "gaussian")
predict(fit, newdata = copd[1:4, ], type = "response")
#> [1] 2.921594 3.125174 2.984774 2.949674
```
