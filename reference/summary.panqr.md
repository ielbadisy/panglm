# Summarize a panqr fit

Summarize a panqr fit

## Usage

``` r
# S3 method for class 'panqr'
summary(object, ...)
```

## Arguments

- object:

  a `"panqr"` fit

- ...:

  unused

## Value

an object of class `"summary.panqr"`: a list with one coefficient table
per quantile (`coefficients`: estimate, standard error, z value and
p-value), the location and scale tables and the quantiles `q` of the
standardized errors for MM-QR, and the fit information (`call`, `model`,
`tau`, `nobs`, `n_groups`, `vcov_type`, `R`, `cluster`, `n_clusters`,
`dropped`).

## See also

Other panqr methods:
[`plot.panqr()`](https://elbadisyimad.com/panglm/reference/plot.panqr.md),
[`vcov.panqr()`](https://elbadisyimad.com/panglm/reference/vcov.panqr.md)

## Examples

``` r
set.seed(2)
d <- data.frame(id = rep(1:100, each = 5), t = rep(1:5, 100))
d$x <- rnorm(500) + rep(rnorm(100), each = 5)
d$y <- d$x + rep(rnorm(100), each = 5) + (1 + 0.2 * abs(d$x)) * rnorm(500)
summary(panqr(y ~ x, data = d, index = c("id", "t"), tau = 0.9, model = "mmqr"))
#> 
#> Call:
#> panqr(formula = y ~ x, data = d, index = c("id", "t"), tau = 0.9, 
#>     model = "mmqr")
#> 
#> Estimator: MM-QR fixed-effects quantile regression (Machado and Santos Silva 2019) 
#> N = 500  groups = 100  vcov: analytic (cluster by individual, G = 100) 
#> 
#> Location (within regression, cluster-robust):
#>   Estimate Std. Error z value Pr(>|z|)    
#> x  0.90372    0.05381    16.8   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Scale (within regression, cluster-robust):
#>   Estimate Std. Error z value Pr(>|z|)
#> x 0.001334   0.031056   0.043    0.966
#> 
#> Quantile tau = 0.9 (q(tau) = 1.685):
#>   Estimate Std. Error z value Pr(>|z|)    
#> x  0.90597    0.07482   12.11   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
```
