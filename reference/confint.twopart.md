# Confidence intervals for a two-part model

Wald intervals for the coefficients of both parts of a
[`twopart()`](https://ielbadisy.github.io/panglm/reference/twopart.md)
fit, named as in [`coef()`](https://rdrr.io/r/stats/coef.html)
(`"binary:..."`, `"positive:..."`), from
[`vcov()`](https://rdrr.io/r/stats/vcov.html) (block diagonal for
independent parts, the joint covariance for `correlated = TRUE`).

## Usage

``` r
# S3 method for class 'twopart'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  a `"twopart"` fit

- parm:

  names or positions of the coefficients; all by default

- level:

  confidence level

- ...:

  unused

## Value

a two-column matrix

## See also

Other twopart methods:
[`predict.twopart()`](https://ielbadisy.github.io/panglm/reference/predict.twopart.md)

## Examples

``` r
data(copd)
copd$cost <- ifelse(copd$exacerbations > 0, exp(6 + 0.1 * copd$crp), 0)
fit <- twopart(cost ~ treatment + crp, data = copd, index = c("id", "visit"))
#> Warning: NaNs produced
confint(fit)
#>                              2.5 %       97.5 %
#> binary:(Intercept)    3.470918e-01 1.350479e+00
#> binary:treatment     -4.799854e-01 1.997329e-01
#> binary:crp           -7.260436e-02 8.511146e-02
#> positive:(Intercept)  6.000000e+00 6.000000e+00
#> positive:treatment    8.387221e-14 1.298222e-13
#> positive:crp          1.000000e-01 1.000000e-01
```
