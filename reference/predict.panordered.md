# Predict from an ordered panel model

Predict from an ordered panel model

## Usage

``` r
# S3 method for class 'panordered'
predict(object, newdata = NULL, type = c("prob", "class", "link"), ...)
```

## Arguments

- object:

  a `"panordered"` fit (from
  [`panglm()`](https://ielbadisy.github.io/panglm/reference/panglm.md)
  with an ordered family)

- newdata:

  optional data frame; defaults to the fitted data

- type:

  `"prob"` for the matrix of category probabilities (population-averaged
  over the random intercept for `model = "random"`), `"class"` for the
  most likely category, or `"link"` for `x'beta`

- ...:

  unused

## Value

a matrix with one column per category (`"prob"`), a factor (`"class"`)
or a numeric vector (`"link"`)

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
[`summary.panglm()`](https://ielbadisy.github.io/panglm/reference/summary.panglm.md),
[`tidy.panglm()`](https://ielbadisy.github.io/panglm/reference/tidy.panglm.md),
[`vcov.panglm()`](https://ielbadisy.github.io/panglm/reference/vcov.panglm.md)

## Examples

``` r
data(copd)
copd$lf <- cut(copd$fev1, 3, labels = c("low", "mid", "high"))
fit <- panglm(lf ~ crp + age, data = copd, index = c("id", "visit"),
              model = "pooling", family = "ologit")
head(predict(fit, type = "prob"))
#>            low       mid       high
#> [1,] 0.3622632 0.5380532 0.09968362
#> [2,] 0.2083513 0.5987692 0.19287953
#> [3,] 0.3091034 0.5676442 0.12325240
#> [4,] 0.3381230 0.5522571 0.10961992
#> [5,] 0.2869118 0.5779029 0.13518524
#> [6,] 0.2197194 0.5977054 0.18257514
predict(fit, newdata = copd[1:3, ], type = "class")
#> [1] mid mid mid
#> Levels: low mid high
```
