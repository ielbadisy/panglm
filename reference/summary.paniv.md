# Summarize a paniv fit

Summarize a paniv fit

## Usage

``` r
# S3 method for class 'paniv'
summary(object, ...)
```

## Arguments

- object:

  a `"paniv"` fit

- ...:

  unused

## Value

an object of class `"summary.paniv"`: a list with the coefficient table
(`coefficients`: estimate, standard error, z value and p-value),
`first_stage`, `overid`, `ercomp`, `endogenous`, `instruments`,
`varlist` (Hausman-Taylor) and the fit information (`call`, `model`,
`nobs`, `n_groups`, `df.residual`, `vcov_type`).

## See also

Other paniv methods:
[`predict.paniv()`](https://elbadisyimad.com/panglm/reference/predict.paniv.md),
[`vcov.paniv()`](https://elbadisyimad.com/panglm/reference/vcov.paniv.md)

## Examples

``` r
set.seed(2)
d <- data.frame(id = rep(1:50, each = 4), t = rep(1:4, 50))
d$z <- rnorm(200); d$u <- rnorm(200)
d$x <- d$z + d$u + rnorm(200); d$y <- 1 + d$x + d$u + rnorm(200)
summary(paniv(y ~ x | z, data = d, index = c("id", "t"), model = "within"))
#> 
#> Call:
#> paniv(formula = y ~ x | z, data = d, index = c("id", "t"), model = "within")
#> 
#> Estimator: fixed-effects 2SLS 
#> N = 200  groups = 50  vcov: classical 
#> Endogenous: x 
#> Excluded instruments: z 
#> 
#> Coefficients:
#>   Estimate Std. Error z value Pr(>|z|)    
#> x  1.14805    0.09381   12.24   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> First stage (excluded instruments):
#>   partial.R2     F df1 df2   p.value
#> x     0.3623 84.66   1 149 2.971e-16
#> Residual df: 149 
#> 
```
