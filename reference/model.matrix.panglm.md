# Design matrix of a panglm fit

Returns the model matrix used for the fit (`data = NULL`) or rebuilt for
new data with the stored terms, factor levels, contrasts and spline
knots ([`splines::ns()`](https://rdrr.io/r/splines/ns.html),
[`splines::bs()`](https://rdrr.io/r/splines/bs.html) and
[`poly()`](https://rdrr.io/r/stats/poly.html) keep the basis of the
fit). Its columns match `coef(object)`. With a regression spline,
`model.matrix(fit, grid) %*% coef(fit)` traces the fitted curve, and the
covariance `vcov(fit)` gives pointwise confidence bands; for
fixed-effects fits the curve is identified up to a constant, so center
it at a reference value.

## Usage

``` r
# S3 method for class 'panglm'
model.matrix(object, data = NULL, cre_means = NULL, ...)
```

## Arguments

- object:

  a `"panglm"` fit

- data:

  optional new data containing the regressors

- cre_means:

  for correlated random effects fits with new data: a matrix of the
  `mean_*` columns, one row per row of `data`; when `NULL` they are
  recomputed as individual means within `data`, which then needs the
  individual index

- ...:

  unused

## Value

a numeric matrix with one column per coefficient

## See also

Other panglm methods:
[`confint.panglm()`](https://elbadisyimad.com/panglm/reference/confint.panglm.md),
[`effectcurve()`](https://elbadisyimad.com/panglm/reference/effectcurve.md),
[`fitted.panglm()`](https://elbadisyimad.com/panglm/reference/fitted.panglm.md),
[`glance.panglm()`](https://elbadisyimad.com/panglm/reference/glance.panglm.md),
[`logLik.panglm()`](https://elbadisyimad.com/panglm/reference/logLik.panglm.md),
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
fit <- panglm(fev1 ~ splines::ns(crp, 3), data = copd,
              index = c("id", "visit"), model = "within")
grid <- data.frame(crp = seq(min(copd$crp), max(copd$crp), length.out = 5))
X <- model.matrix(fit, grid)
X0 <- model.matrix(fit, data.frame(crp = median(copd$crp)))
Xc <- sweep(X, 2, X0)            # centered at the median CRP
cbind(grid, effect = drop(Xc %*% coef(fit)),
      se = sqrt(rowSums((Xc %*% vcov(fit)) * Xc)))
#>    crp      effect         se
#> 1  0.1  0.17651379 0.06824028
#> 2  3.4  0.04928458 0.01512769
#> 3  6.7 -0.04884615 0.01859332
#> 4 10.0 -0.21194570 0.05430171
#> 5 13.3 -0.42355337 0.14824802
```
