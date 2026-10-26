# Summarize a pangmm fit

Summarize a pangmm fit

## Usage

``` r
# S3 method for class 'pangmm'
summary(object, time.dummies = FALSE, ...)
```

## Arguments

- object:

  a `"pangmm"` fit

- time.dummies:

  logical: show the time effects (and the system-GMM constant) in the
  coefficient table

- ...:

  unused

## Value

an object of class `"summary.pangmm"`: a list with the coefficient table
(`coefficients`: estimate, standard error, z value and p-value), the
Arellano-Bond tests (`artest`), `sargan`, `hansen`, `wald` and the fit
information (`call`, `label`, `vcov_type`, `nobs`, `n_groups`,
`n_periods`, `n_instruments`).

## See also

Other pangmm methods:
[`predict.pangmm()`](https://ielbadisy.github.io/panglm/reference/predict.pangmm.md),
[`vcov.pangmm()`](https://ielbadisy.github.io/panglm/reference/vcov.pangmm.md)

## Examples

``` r
set.seed(3)
N <- 200; TT <- 6
d <- data.frame(id = rep(seq_len(N), each = TT), t = rep(seq_len(TT), N))
mu <- rep(rnorm(N), each = TT)
d$y <- ave(mu + rnorm(N * TT), d$id, FUN = function(e) stats::filter(e, 0.4,
  method = "recursive"))
summary(pangmm(y ~ lag(y, 1) | lag(y, 2:99), data = d, index = c("id", "t")))
#> 
#> Call:
#> pangmm(formula = y ~ lag(y, 1) | lag(y, 2:99), data = d, index = c("id", 
#>     "t"))
#> 
#> Estimator: two-way one-step difference GMM (Arellano-Bond) 
#> N = 200  T = 6  equations used = 800  instruments = 14  vcov: robust 
#> 
#> Coefficients:
#>           Estimate Std. Error z value Pr(>|z|)   
#> lag(y, 1)   0.3645     0.1268   2.875  0.00404 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Arellano-Bond AR(1) test: z = -4.564, p = 5.018e-06
#> Arellano-Bond AR(2) test: z = -0.4769, p = 0.6334
#> Sargan test (not robust): chi2(9) = 7.899, p = 0.5443
#> Hansen J test (robust): chi2(9) = 9.133, p = 0.4251
#> Wald test for coefficients: chi2(1) = 8.264, p = 0.004043
#> 
```
