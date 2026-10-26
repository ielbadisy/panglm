# Methods for event-study fits

[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html),
[`confint()`](https://rdrr.io/r/stats/confint.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) for the fits of
[`eventstudy()`](https://ielbadisy.github.io/panglm/reference/eventstudy.md).
The coefficients are labelled by event time, as in
[`coef()`](https://rdrr.io/r/stats/coef.html). The reference period
enters as a coefficient fixed at zero, with zero variance and covariance
and a degenerate interval at zero; an event time that the data do not
identify has missing values throughout. Intervals and tests are Wald
(normal) based, from the variance of the underlying two-way
fixed-effects fit (`object$fit`).

## Usage

``` r
# S3 method for class 'panglm_eventstudy'
vcov(object, ...)

# S3 method for class 'panglm_eventstudy'
nobs(object, ...)

# S3 method for class 'panglm_eventstudy'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'panglm_eventstudy'
summary(object, level = object$level, ...)

# S3 method for class 'summary.panglm_eventstudy'
print(x, digits = max(3, getOption("digits") - 3), signif.stars = TRUE, ...)

# S3 method for class 'panglm_eventstudy'
logLik(object, ...)
```

## Arguments

- object:

  a `"panglm_eventstudy"` fit from
  [`eventstudy()`](https://ielbadisy.github.io/panglm/reference/eventstudy.md)

- ...:

  unused

- parm:

  event times (as numbers or names) or positions of the coefficients;
  all by default

- level:

  confidence level

- x:

  a `"summary.panglm_eventstudy"` object

- digits:

  significant digits to print

- signif.stars:

  passed to
  [`stats::printCoefmat()`](https://rdrr.io/r/stats/printCoefmat.html)

## Value

[`vcov()`](https://rdrr.io/r/stats/vcov.html): a square matrix with one
row per event time; [`nobs()`](https://rdrr.io/r/stats/nobs.html): the
number of observations of the fit;
[`confint()`](https://rdrr.io/r/stats/confint.html): a two-column
matrix; [`summary()`](https://rdrr.io/r/base/summary.html): an object of
class `"summary.panglm_eventstudy"` with the coefficient table
`coefficients` (estimate, standard error, z value, p-value), the
intervals `conf.int`, the pre-trend test `pretrend` and the fit
description.

## See also

Other eventstudy methods:
[`plot.panglm_eventstudy()`](https://ielbadisy.github.io/panglm/reference/plot.panglm_eventstudy.md)

## Examples

``` r
set.seed(1)
d <- expand.grid(id = 1:200, year = 2015:2022)
d$first <- ifelse(d$id <= 100, 2019, NA)
d$visits <- 3 + rnorm(200)[d$id] + 0.1 * (d$year - 2015) +
  ifelse(!is.na(d$first) & d$year >= d$first, 0.5, 0) + rnorm(nrow(d))
es <- eventstudy(visits ~ 1, data = d, index = c("id", "year"),
                 treat = "first", leads = 3, lags = 2)
summary(es)
#> 
#> Event study (two-way fixed effects, gaussian)
#> N = 1600  groups = 200  vcov: cluster 
#> Reference event time: -1 (end points binned) 
#> 
#> Event-time coefficients:
#>    Estimate Std. Error z value Pr(>|z|)    
#> -3 -0.01428    0.17845  -0.080  0.93620    
#> -2 -0.20585    0.21162  -0.973  0.33068    
#> -1  0.00000    0.00000                     
#> 0   0.63560    0.21502   2.956  0.00312 ** 
#> 1   0.90970    0.20631   4.409 1.04e-05 ***
#> 2   0.53888    0.17496   3.080  0.00207 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Pre-trend test: Chi-squared = 1.24, df = 2, p-value = 0.5379
confint(es)
#>         2.5 %    97.5 %
#> -3 -0.3640365 0.3354688
#> -2 -0.6206051 0.2089111
#> -1  0.0000000 0.0000000
#> 0   0.2141671 1.0570314
#> 1   0.5053292 1.3140620
#> 2   0.1959669 0.8817952
vcov(es)
#>            -3         -2 -1          0          1          2
#> -3 0.03184387 0.02070525  0 0.02334587 0.02237814 0.02124967
#> -2 0.02070525 0.04478098  0 0.02352476 0.02424383 0.02216843
#> -1 0.00000000 0.00000000  0 0.00000000 0.00000000 0.00000000
#> 0  0.02334587 0.02352476  0 0.04623375 0.02622856 0.02300849
#> 1  0.02237814 0.02424383  0 0.02622856 0.04256513 0.02063024
#> 2  0.02124967 0.02216843  0 0.02300849 0.02063024 0.03061079
```
