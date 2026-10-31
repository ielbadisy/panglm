# Predictions from a panel instrumental-variables fit

Linear predictions from a
[`paniv()`](https://elbadisyimad.com/panglm/reference/paniv.md) fit. For
`model = "pooling"`, `"random"` and `"hausmantaylor"` the prediction is
the linear index `x'b` (with the intercept), the conditional mean given
the regressors that does not involve the individual effect; this is what
[`predict()`](https://rdrr.io/r/stats/predict.html) returns for
[`AER::ivreg()`](https://rdrr.io/pkg/AER/man/ivreg.html) and
[`plm::plm()`](https://rdrr.io/pkg/plm/man/plm.html) fits of these
models. For `model = "within"` the estimated individual effects
`a_i = mean_t(y_it - x_it'b)` are added (and, for `effect = "twoways"`,
the time effects, from the two-way decomposition of `y - x'b`), as in
[`predict.panglm()`](https://elbadisyimad.com/panglm/reference/predict.panglm.md).

## Usage

``` r
# S3 method for class 'paniv'
predict(object, newdata = NULL, allow.new.levels = FALSE, ...)
```

## Arguments

- object:

  a `"paniv"` fit

- newdata:

  optional data frame

- allow.new.levels:

  for within fits, if `TRUE` an individual (or time period) of `newdata`
  that is not in the fit gets an effect of zero, so that its prediction
  is the linear index `x'b`; if `FALSE` (default) it is an error

- ...:

  unused

## Value

a numeric vector of predictions, one per row of `newdata` (`NA` where a
regressor is missing)

## Details

Without `newdata` the in-sample fitted values are returned in the row
order of the fitting data (with `NA` for the rows dropped by
`na.action`). With `newdata`, the regressors of the first part of the
formula are evaluated on `newdata`; factor regressors must carry the
levels used in the fit. For within fits `newdata` must contain the panel
index columns.

## Validation

Numerical accuracy is checked in the package tests and in the validation
suite; see the validation article on the package website.

## See also

Other paniv methods:
[`summary.paniv()`](https://elbadisyimad.com/panglm/reference/summary.paniv.md),
[`vcov.paniv()`](https://elbadisyimad.com/panglm/reference/vcov.paniv.md)

## Examples

``` r
set.seed(1)
d <- data.frame(id = rep(1:50, each = 4), t = rep(1:4, 50))
d$z <- rnorm(200); d$u <- rnorm(200)
d$x <- d$z + d$u + rnorm(200)
d$y <- 1 + d$x + d$u + rep(rnorm(50), each = 4) + rnorm(200)
fit <- paniv(y ~ x | z, data = d, index = c("id", "t"), model = "within")
head(predict(fit))
#> [1]  2.1897524  5.4910216  1.3851551  2.2080306 -0.3994349  3.7668555
predict(fit, newdata = data.frame(id = 1:2, t = 5, x = 0))
#> [1] 1.217343 2.290716
```
