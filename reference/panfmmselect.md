# Compare latent class panel models with different numbers of classes

Fits [`panfmm()`](https://elbadisyimad.com/panglm/reference/panfmm.md)
for each number of classes in `k` and tabulates the log-likelihood,
number of parameters, AIC and BIC (with `n` the number of observations),
the usual criteria for choosing the number of classes (Deb and Trivedi
2002 prefer BIC). Mixtures with too many classes may fail to converge or
return near-empty classes; such fits are reported with `NA`.

## Usage

``` r
panfmmselect(
  formula,
  data,
  index,
  k = 1:3,
  family = "poisson",
  concomitant = NULL,
  ...
)
```

## Arguments

- formula, data, index, family, concomitant, ...:

  passed to
  [`panfmm()`](https://elbadisyimad.com/panglm/reference/panfmm.md)

- k:

  integer vector of numbers of classes

## Value

a data.frame with columns `k`, `logLik`, `df`, `AIC`, `BIC` and
`smallest` (the smallest mean class probability), with the fits in the
attribute `"fits"`

## Validation

Each fit is a
[`panfmm()`](https://elbadisyimad.com/panglm/reference/panfmm.md) fit,
whose log-likelihood and BIC are checked against flexmix in the
package's validation suite on GitHub.

## See also

Other panel models:
[`panfmm()`](https://elbadisyimad.com/panglm/reference/panfmm.md),
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md),
[`pangmm()`](https://elbadisyimad.com/panglm/reference/pangmm.md),
[`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md),
[`paniv()`](https://elbadisyimad.com/panglm/reference/paniv.md),
[`panqr()`](https://elbadisyimad.com/panglm/reference/panqr.md),
[`panselect()`](https://elbadisyimad.com/panglm/reference/panselect.md),
[`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md)

## Examples

``` r
set.seed(3)
d <- data.frame(id = rep(1:150, each = 3), t = rep(1:3, 150), x = rnorm(450))
d$y <- rpois(450, exp(ifelse(rep(rbinom(150, 1, 0.4), each = 3) == 1, 1.5, 0) + 0.3 * d$x))
panfmmselect(y ~ x, data = d, index = c("id", "t"), k = 1:2, nstart = 2, seed = 1)
#>   k     logLik df      AIC      BIC  smallest
#> 1 1 -1004.9506  2 2013.901 2022.120 1.0000000
#> 2 2  -799.9038  5 1609.808 1630.354 0.3791781
# \donttest{
# up to three classes (a few seconds)
panfmmselect(y ~ x, data = d, index = c("id", "t"), k = 1:3, nstart = 2, seed = 1)
#>   k     logLik df      AIC      BIC  smallest
#> 1 1 -1004.9506  2 2013.901 2022.120 1.0000000
#> 2 2  -799.9038  5 1609.808 1630.354 0.3791781
#> 3 3  -798.7190  8 1613.438 1646.312 0.1170190
# }
```
