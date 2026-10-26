# Extract the model log-likelihood

The `df` attribute counts every parameter associated with the reported
likelihood, including fitted shape, dispersion, variance, and explicit
fixed-effect parameters. For conditional Poisson and conditional
logistic likelihoods, only parameters remaining in the conditional
likelihood are counted. Models estimated without a likelihood return
`NA`.

## Usage

``` r
# S3 method for class 'panglm'
logLik(object, ...)
```

## Arguments

- object:

  a fitted `"panglm"` object

- ...:

  unused

## Value

an object of class `"logLik"`

## See also

Other panglm methods:
[`confint.panglm()`](https://elbadisyimad.com/panglm/reference/confint.panglm.md),
[`effectcurve()`](https://elbadisyimad.com/panglm/reference/effectcurve.md),
[`fitted.panglm()`](https://elbadisyimad.com/panglm/reference/fitted.panglm.md),
[`glance.panglm()`](https://elbadisyimad.com/panglm/reference/glance.panglm.md),
[`model.matrix.panglm()`](https://elbadisyimad.com/panglm/reference/model.matrix.panglm.md),
[`partialeffects.panordered()`](https://elbadisyimad.com/panglm/reference/partialeffects.md),
[`plot.panglm()`](https://elbadisyimad.com/panglm/reference/plot.panglm.md),
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
logLik(fit)
#> 'log Lik.' -604.8321 (df=2)
AIC(fit)
#> [1] 1213.664
BIC(fit)
#> [1] 1221.072
```
