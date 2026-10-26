# Summarize a panfmm fit

Summarize a panfmm fit

## Usage

``` r
# S3 method for class 'panfmm'
summary(object, ...)
```

## Arguments

- object:

  a `"panfmm"` fit

- ...:

  unused

## Value

an object of class `"summary.panfmm"`: a list with one coefficient table
per class (`classes`: estimate, standard error, z value and p-value,
ancillary parameter on the log scale), the table of the
class-probability (multinomial logit) parameters (`prior`), the mean
class probabilities, the class sizes by modal posterior class, and the
fit information (`call`, `loglik`, `AIC`, `BIC`, `nobs`, `n_groups`,
`vcov_type`).

## See also

Other panfmm methods:
[`predict.panfmm()`](https://ielbadisy.github.io/panglm/reference/predict.panfmm.md),
[`vcov.panfmm()`](https://ielbadisy.github.io/panglm/reference/vcov.panfmm.md)

## Examples

``` r
set.seed(3)
d <- data.frame(id = rep(1:150, each = 3), t = rep(1:3, 150), x = rnorm(450))
d$y <- rpois(450, exp(ifelse(rep(rbinom(150, 1, 0.4), each = 3) == 1, 1.5, 0) + 0.3 * d$x))
summary(panfmm(y ~ x, data = d, index = c("id", "t"), nstart = 2, seed = 1))
#> 
#> Call:
#> panfmm(formula = y ~ x, data = d, index = c("id", "t"), nstart = 2, 
#>     seed = 1)
#> 
#> Model: 2-class latent class Poisson model 
#> N = 450  groups = 150  vcov: classical 
#> 
#> Class 1 (prior 0.621, 94 individuals by modal posterior, mean outcome 0.9923):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -0.04987    0.06919  -0.721    0.471    
#> x            0.24535    0.05846   4.197 2.71e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Class 2 (prior 0.379, 56 individuals by modal posterior, mean outcome 4.71):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.48145    0.04063  36.460   <2e-16 ***
#> x            0.32210    0.03908   8.241   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Class probabilities (multinomial logit, class 1 = reference):
#>                    Estimate Std. Error z value Pr(>|z|)   
#> prior2:(Intercept)  -0.4930     0.1758  -2.804  0.00504 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> logLik = -799.9038  df = 5  AIC = 1609.808  BIC = 1630.354 
#> 
```
