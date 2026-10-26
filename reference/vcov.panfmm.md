# Covariance of a panfmm fit

Covariance of a panfmm fit

## Usage

``` r
# S3 method for class 'panfmm'
vcov(object, type = NULL, cluster = NULL, ...)
```

## Arguments

- object:

  a `"panfmm"` fit

- type:

  `NULL` (the covariance stored on the fit), `"classical"` (inverse of
  the negative Hessian of the observed-data log-likelihood) or
  `"cluster"` (sandwich with individual scores summed within clusters,
  factor `G / (G - 1)`; with weights the individual scores are
  multiplied by the individual weights)

- cluster:

  optional clustering vector with one value per observation used in the
  fit (in the row order of the data), constant within individuals;
  defaults to the individual

- ...:

  unused

## Value

the covariance matrix of all parameters

## See also

Other panfmm methods:
[`predict.panfmm()`](https://ielbadisy.github.io/panglm/reference/predict.panfmm.md),
[`summary.panfmm()`](https://ielbadisy.github.io/panglm/reference/summary.panfmm.md)

## Examples

``` r
set.seed(3)
d <- data.frame(id = rep(1:150, each = 3), t = rep(1:3, 150), x = rnorm(450))
d$y <- rpois(450, exp(ifelse(rep(rbinom(150, 1, 0.4), each = 3) == 1, 1.5, 0) + 0.3 * d$x))
fit <- panfmm(y ~ x, data = d, index = c("id", "t"), nstart = 2, seed = 1)
sqrt(diag(vcov(fit, type = "cluster")))
#> class1:(Intercept)           class1:x class2:(Intercept)           class2:x 
#>         0.07390899         0.05899464         0.04201631         0.03487727 
#> prior2:(Intercept) 
#>         0.17663282 
```
