# Summarize a panglm fit

Produces a classic coefficient table (Estimate, Std. Error, z value,
Pr(\>\|z\|)) in the style of
[`summary.glm()`](https://rdrr.io/r/stats/summary.glm.html).

## Usage

``` r
# S3 method for class 'panglm'
summary(object, ...)
```

## Arguments

- object:

  a `"panglm"` object

- ...:

  unused

## Value

an object of class `"summary.panglm"`, a list with components `call`
(the matched call), `model`, `effect`, and `family` (as stored on the
fit), `coefficients` (a numeric matrix with one row per parameter and
columns `Estimate`, `Std. Error`, `z value`, `Pr(>|z|)`), `nobs`,
`n_groups`, `df.residual`, `loglik`, `dispersion` (`NULL` when not
applicable, as for negative binomial fits), `iterations`, and
`vcov_type` (which covariance estimator produced the reported standard
errors). It has an associated `print.summary.panglm()` method and is
normally printed rather than used programmatically. Also included are
`sigma` (random-intercept standard deviation) and `shape` (gamma shape)
when estimated, `varcomp` (covariance matrix of the random intercept and
slopes) for random-slope fits, for `family = "negbin"` the NB2
`negbin_theta` (variance `mu + mu^2 / theta`) of pooled, fixed-effects
and random-slope fits, with its standard error `negbin_theta_se` (pooled
and fixed-effects fits, from the observed information at the fitted
means, as in
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html)), or
`negbin_ab`, the beta mixing parameters `a` and `b` of the
random-intercept model, and, for `cre = TRUE` fits, `cre_test`: the Wald
(Mundlak) test that the individual-mean coefficients are jointly zero,
computed with the covariance stored on the fit.

## See also

Other panglm methods:
[`confint.panglm()`](https://elbadisyimad.com/panglm/reference/confint.panglm.md),
[`effectcurve()`](https://elbadisyimad.com/panglm/reference/effectcurve.md),
[`fitted.panglm()`](https://elbadisyimad.com/panglm/reference/fitted.panglm.md),
[`glance.panglm()`](https://elbadisyimad.com/panglm/reference/glance.panglm.md),
[`logLik.panglm()`](https://elbadisyimad.com/panglm/reference/logLik.panglm.md),
[`model.matrix.panglm()`](https://elbadisyimad.com/panglm/reference/model.matrix.panglm.md),
[`partialeffects.panordered()`](https://elbadisyimad.com/panglm/reference/partialeffects.md),
[`plot.panglm()`](https://elbadisyimad.com/panglm/reference/plot.panglm.md),
[`plot.panglm_effectcurve()`](https://elbadisyimad.com/panglm/reference/plot.panglm_effectcurve.md),
[`predict.panglm()`](https://elbadisyimad.com/panglm/reference/predict.panglm.md),
[`predict.panordered()`](https://elbadisyimad.com/panglm/reference/predict.panordered.md),
[`tidy.panglm()`](https://elbadisyimad.com/panglm/reference/tidy.panglm.md),
[`vcov.panglm()`](https://elbadisyimad.com/panglm/reference/vcov.panglm.md)

## Examples

``` r
data(copd)
fit <- panglm(fev1 ~ treatment + age + smoker + crp, data = copd,
              index = c("id", "visit"), model = "pooling", family = "gaussian")
summary(fit)
#> 
#> Call:
#> panglm(formula = fev1 ~ treatment + age + smoker + crp, data = copd, 
#>     index = c("id", "visit"), model = "pooling", family = "gaussian")
#> 
#> Family: gaussian  Link: identity  Model: pooling  Effect: individual 
#> N = 300  groups = 60  vcov: classical 
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  4.091479   0.225597  18.136  < 2e-16 ***
#> treatment    0.208556   0.056578   3.686 0.000228 ***
#> age         -0.014239   0.003133  -4.545  5.5e-06 ***
#> smoker      -0.355207   0.055905  -6.354  2.1e-10 ***
#> crp         -0.032895   0.013135  -2.504 0.012269 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -189.6 
#> Dispersion parameter: 0.2107 
#> Residual df: 295  Iterations: 2
#> 
```
