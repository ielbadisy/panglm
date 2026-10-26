# Covariance of a paniv fit

Covariance of a paniv fit

## Usage

``` r
# S3 method for class 'paniv'
vcov(object, type = NULL, cluster = NULL, ...)
```

## Arguments

- object:

  a `"paniv"` fit

- type:

  `NULL` (the covariance stored on the fit), `"classical"` or
  `"cluster"`

- cluster:

  optional clustering vector with one value per observation used in the
  fit (in the row order of the data); defaults to the panel individual

- ...:

  unused

## Value

the covariance matrix of the coefficients

## See also

Other paniv methods:
[`predict.paniv()`](https://ielbadisy.github.io/panglm/reference/predict.paniv.md),
[`summary.paniv()`](https://ielbadisy.github.io/panglm/reference/summary.paniv.md)

## Examples

``` r
set.seed(2)
d <- data.frame(id = rep(1:50, each = 4), t = rep(1:4, 50))
d$z <- rnorm(200); d$u <- rnorm(200)
d$x <- d$z + d$u + rnorm(200); d$y <- 1 + d$x + d$u + rnorm(200)
fit <- paniv(y ~ x | z, data = d, index = c("id", "t"), model = "pooling")
vcov(fit, type = "cluster")
#>              (Intercept)            x
#> (Intercept)  0.006622603 -0.001131347
#> x           -0.001131347  0.004784531
```
