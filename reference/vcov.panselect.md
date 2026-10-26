# Covariance of a panselect fit

Covariance of a panselect fit

## Usage

``` r
# S3 method for class 'panselect'
vcov(object, type = NULL, ...)
```

## Arguments

- object:

  a `"panselect"` fit

- type:

  `NULL` (the default covariance of the fit), `"boot"` (panel
  bootstrap), `"cluster"` (cluster-robust by individual, ignoring the
  estimation of the first step) or `"classical"` (homoskedastic, also
  ignoring the first step)

- ...:

  unused

## Value

the covariance matrix of the outcome-equation coefficients

## See also

Other panselect methods:
[`predict.panselect()`](https://elbadisyimad.com/panglm/reference/predict.panselect.md),
[`summary.panselect()`](https://elbadisyimad.com/panglm/reference/summary.panselect.md)

## Examples

``` r
set.seed(3)
d <- data.frame(id = rep(1:100, each = 3), t = rep(1:3, 100))
d$x <- rnorm(300); d$z <- rnorm(300); v <- rnorm(300)
d$s <- as.numeric(0.5 + d$x + d$z + v > 0)
d$y <- ifelse(d$s == 1, 1 + d$x + 0.5 * v + rnorm(300), NA)
fit <- panselect(s ~ x + z, y ~ x, data = d, index = c("id", "t"), boot = 0)
vcov(fit, type = "cluster")
#>                   x   lambda:1   lambda:2   lambda:3
#> x        0.01969311 0.02560123 0.02623059 0.03407757
#> lambda:1 0.02560123 0.22417871 0.09592534 0.09934233
#> lambda:2 0.02623059 0.09592534 0.23396494 0.06628292
#> lambda:3 0.03407757 0.09934233 0.06628292 0.24408854
```
