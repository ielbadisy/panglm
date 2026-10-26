# Extract fitted values from a panglm model

Values are returned in the row order of the data supplied to
[`panglm()`](https://ielbadisy.github.io/panglm/reference/panglm.md),
even though estimation internally sorts observations by panel index.
Fixed-effects binomial values are conditional inclusion probabilities,
conditional on each individual's observed number of successes.

## Usage

``` r
# S3 method for class 'panglm'
fitted(object, ...)
```

## Arguments

- object:

  a fitted `"panglm"` object

- ...:

  unused

## Value

a numeric vector of fitted values, one per row of the original data
supplied to
[`panglm()`](https://ielbadisy.github.io/panglm/reference/panglm.md) (in
that row order), or `NA` at rows dropped for missingness, following the
model's `na.action`.

## See also

Other panglm methods:
[`confint.panglm()`](https://ielbadisy.github.io/panglm/reference/confint.panglm.md),
[`effectcurve()`](https://ielbadisy.github.io/panglm/reference/effectcurve.md),
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
fit <- panglm(exacerbations ~ crp + treatment, data = copd,
              index = c("id", "visit"), model = "random", family = "poisson")
head(fitted(fit))
#> [1] 1.832969 1.472463 1.712529 1.778428 1.661572 1.719008
```
