# Predictions from a dynamic panel GMM fit

Linear predictions from a
[`pangmm()`](https://elbadisyimad.com/panglm/reference/pangmm.md) fit on
`newdata`, a long data frame with the panel index columns and the
variables of the first part of the formula. Lags are taken within
individuals over the sorted periods of `newdata`, as in the fit.

## Usage

``` r
# S3 method for class 'pangmm'
predict(object, newdata, type = c("level", "difference"), ...)
```

## Arguments

- object:

  a `"pangmm"` fit

- newdata:

  a data frame in long format; the fit does not store its data, so
  `newdata` is required (pass the fitting data for in-sample
  predictions)

- type:

  `"level"` or `"difference"`

- ...:

  unused

## Value

a numeric vector with one prediction per row of `newdata`

## Details

- `type = "level"`: `x_it'b + d_t (+ c)`, the regressors (with their
  lags, including the lagged dependent variable) times their
  coefficients, plus the time effect and, for system GMM, the constant.
  The individual effect is not identified by GMM and is excluded, so for
  difference GMM (no constant) this is the level up to an individual
  intercept. On the fitting data these are the fitted values of the
  level equations of system GMM.

- `type = "difference"`: the prediction of the first difference
  `y_it - y_i,t-1`, the difference of the level predictions, in which
  the individual effect cancels. On the fitting data these are the
  fitted values of the differenced equations.

Periods that precede the first time effect of the fit have a time effect
of zero (the reference); a later period that is not in the fit gives a
missing prediction, with a warning. A prediction is missing whenever a
required lag is unavailable.

## Validation

Numerical accuracy is checked in the package tests and in the validation
suite; see the validation article on the package website.

## See also

Other pangmm methods:
[`summary.pangmm()`](https://elbadisyimad.com/panglm/reference/summary.pangmm.md),
[`vcov.pangmm()`](https://elbadisyimad.com/panglm/reference/vcov.pangmm.md)

## Examples

``` r
set.seed(3)
N <- 100; TT <- 6
d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
d$x <- rnorm(N * TT)
d$y <- ave(rep(rnorm(N), each = TT) + 0.3 * d$x + rnorm(N * TT), d$id,
           FUN = function(e) stats::filter(e, 0.4, method = "recursive"))
fit <- pangmm(y ~ lag(y, 1) + x | lag(y, 2:99), data = d, index = c("id", "t"),
              transformation = "ld")
head(predict(fit, newdata = d))
#> [1]         NA -0.6878707 -0.8761411 -1.0115302  0.1196619 -0.5116195
head(predict(fit, newdata = d, type = "difference"))
#> [1]         NA         NA -0.1882704 -0.1353892  1.1311921 -0.6312814
```
