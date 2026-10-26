# Summarize a panselect fit

Summarize a panselect fit

## Usage

``` r
# S3 method for class 'panselect'
summary(object, ...)
```

## Arguments

- object:

  a `"panselect"` fit

- ...:

  unused

## Value

an object of class `"summary.panselect"`: a list with the coefficient
table (`coefficients`: estimate, standard error, z value and p-value,
using the default covariance of the fit), the selection test
(`selection_test`), the probit coefficients (`probit`) and the fit
information.

## See also

Other panselect methods:
[`predict.panselect()`](https://ielbadisy.github.io/panglm/reference/predict.panselect.md),
[`vcov.panselect()`](https://ielbadisy.github.io/panglm/reference/vcov.panselect.md)

## Examples

``` r
set.seed(3)
d <- data.frame(id = rep(1:100, each = 3), t = rep(1:3, 100))
d$x <- rnorm(300); d$z <- rnorm(300); v <- rnorm(300)
d$s <- as.numeric(0.5 + d$x + d$z + v > 0)
d$y <- ifelse(d$s == 1, 1 + d$x + 0.5 * v + rnorm(300), NA)
summary(panselect(s ~ x + z, y ~ x, data = d, index = c("id", "t"),
                  boot = 29, seed = 1))
#> 
#> Call:
#> panselect(selection = s ~ x + z, outcome = y ~ x, data = d, index = c("id", 
#>     "t"), boot = 29, seed = 1)
#> 
#> Wooldridge (1995) panel selection correction: fixed effects with Mundlak means, inverse Mills ratio by period 
#> Selection sample: N = 300  groups = 100 
#> Outcome sample:   N = 193  groups = 94 
#> Probits fitted for periods: 1, 2, 3 
#> 
#> Outcome equation (vcov: panel bootstrap, 29 replicates):
#>          Estimate Std. Error z value Pr(>|z|)    
#> x          0.9704     0.1243   7.809 5.78e-15 ***
#> lambda:1   0.3577     0.5316   0.673    0.501    
#> lambda:2   0.1731     0.5164   0.335    0.737    
#> lambda:3   0.5822     0.4667   1.247    0.212    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Test of selection bias (H0: all inverse Mills ratio coefficients are zero):
#>   cluster-robust Wald chi2(3) = 1.473, p = 0.6885
#>   bootstrap Wald chi2(3) = 1.621, p = 0.6545
#> Residual df: 95 
#> 
```
