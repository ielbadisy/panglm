# Coefficient methods for the hurdle model

[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) for hurdle fits
([`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md)).
The coefficients of the zero part (conditional logit of `1(y > 0)`) and
of the count part (zero-truncated count model on the positive
observations) are named `"zero:..."` and `"count:..."`. The two parts
are estimated from separate factors of the hurdle likelihood, so
[`vcov()`](https://rdrr.io/r/stats/vcov.html) is block diagonal. No
[`logLik()`](https://rdrr.io/r/stats/logLik.html) method is provided:
the zero part maximises a conditional likelihood and the count part an
unconditional likelihood with individual intercepts, so their sum is not
the log-likelihood of one model.

## Usage

``` r
# S3 method for class 'panglm_hurdle'
coef(object, ...)

# S3 method for class 'panglm_hurdle'
vcov(object, ...)

# S3 method for class 'panglm_hurdle'
nobs(object, ...)

# S3 method for class 'panglm_hurdle'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'panglm_hurdle'
summary(object, ...)

# S3 method for class 'summary.panglm_hurdle'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- object:

  a `"panglm_hurdle"` fit

- ...:

  unused

- parm:

  names or positions of the coefficients; all by default

- level:

  confidence level

- x:

  a `"summary.panglm_hurdle"` object

- digits:

  significant digits to print

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html): a named vector;
[`vcov()`](https://rdrr.io/r/stats/vcov.html): a matrix;
[`confint()`](https://rdrr.io/r/stats/confint.html): a two-column
matrix; [`nobs()`](https://rdrr.io/r/stats/nobs.html): the number of
observations; [`summary()`](https://rdrr.io/r/base/summary.html): an
object of class `"summary.panglm_hurdle"` with one coefficient table per
part (`zero`, `count`: estimate, standard error, z value and p-value,
from the inverse Hessian of each part) and the fit information.

## See also

Other panhurdle methods:
[`plot.panglm_hurdle()`](https://elbadisyimad.com/panglm/reference/plot.panglm_hurdle.md)

## Examples

``` r
data(copd)
fit <- panhurdle(exacerbations ~ crp, data = copd, index = c("id", "visit"))
#> 26 group(s) with no within-group outcome variation dropped (all-0 or all-1) -- they carry no information for the conditional likelihood.
coef(fit)
#>   zero:crp  count:crp 
#> 0.10889286 0.04246066 
confint(fit)
#>                 2.5 %     97.5 %
#> zero:crp  -0.05129363 0.26907935
#> count:crp -0.00624235 0.09116367
summary(fit)
#> 
#> Formula: exacerbations ~ crp
#> 
#> Estimator: fixed-effects hurdle model (conditional logit and zero-truncated Poisson)
#> N = 300  positive = 238  vcov: model-based (inverse Hessian)
#> 
#> Zero part, 1(y > 0), conditional logit (groups used = 34):
#>     Estimate Std. Error z value Pr(>|z|)
#> crp  0.10889    0.08173   1.332    0.183
#> 
#> Count part, y | y > 0, zero-truncated Poisson (groups used = 59):
#>     Estimate Std. Error z value Pr(>|z|)  
#> crp  0.04246    0.02485   1.709   0.0875 .
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
```
