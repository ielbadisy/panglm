# Covariance of a pangmm fit

Covariance of a pangmm fit

## Usage

``` r
# S3 method for class 'pangmm'
vcov(object, type = NULL, ...)
```

## Arguments

- object:

  a `"pangmm"` fit

- type:

  `NULL` (the covariance stored on the fit), `"robust"` (the one-step
  sandwich, or the Windmeijer-corrected two-step covariance) or
  `"classical"` (the homoskedastic one-step covariance, or the
  uncorrected two-step covariance); `"cluster"` is a synonym of
  `"robust"`

- ...:

  unused

## Value

the covariance matrix of the coefficients

## See also

Other pangmm methods:
[`predict.pangmm()`](https://elbadisyimad.com/panglm/reference/predict.pangmm.md),
[`summary.pangmm()`](https://elbadisyimad.com/panglm/reference/summary.pangmm.md)

## Examples

``` r
set.seed(3)
N <- 200; TT <- 6
d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
mu <- rep(rnorm(N), each = TT)
d$y <- ave(mu + rnorm(N * TT), d$id, FUN = function(e) stats::filter(e, 0.4,
  method = "recursive"))
fit <- pangmm(y ~ lag(y, 1) | lag(y, 2:99), data = d, index = c("id", "t"),
              model = "twostep")
sqrt(diag(vcov(fit, type = "classical")))
#>  lag(y, 1)          3          4          5          6 
#> 0.12466147 0.09332815 0.10048407 0.10112917 0.09340707 
sqrt(diag(vcov(fit, type = "robust")))
#>  lag(y, 1)          3          4          5          6 
#> 0.12523344 0.09553760 0.10970517 0.11368133 0.09659418 
```
