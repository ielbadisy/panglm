# Covariance of a panqr fit

Covariance of a panqr fit

## Usage

``` r
# S3 method for class 'panqr'
vcov(object, tau = NULL, ...)
```

## Arguments

- object:

  a `"panqr"` fit

- tau:

  `NULL` or one of the fitted quantiles

- ...:

  unused

## Value

the covariance matrix of the coefficients at quantile `tau`; when
several quantiles were fitted and `tau` is `NULL`, a list of matrices,
one per quantile

## See also

Other panqr methods:
[`plot.panqr()`](https://ielbadisy.github.io/panglm/reference/plot.panqr.md),
[`summary.panqr()`](https://ielbadisy.github.io/panglm/reference/summary.panqr.md)

## Examples

``` r
set.seed(2)
d <- data.frame(id = rep(1:100, each = 5), t = rep(1:5, 100))
d$x <- rnorm(500) + rep(rnorm(100), each = 5)
d$y <- d$x + rep(rnorm(100), each = 5) + (1 + 0.2 * abs(d$x)) * rnorm(500)
fit <- panqr(y ~ x, data = d, index = c("id", "t"), tau = c(0.25, 0.75),
             model = "mmqr")
vcov(fit, tau = 0.75)
#>             x
#> x 0.003559133
```
