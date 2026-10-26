# Two-part model for panel cost data

Fits the two-part model of Duan et al. (1983), the standard model for
health care expenditures with a mass of zeros: a binary model for
`P(y > 0)` on all observations, and a model for `E[y | y > 0]` on the
positive observations. The unconditional mean is
`E[y] = P(y > 0) * E[y | y > 0]`.

## Usage

``` r
twopart(
  formula,
  data,
  index,
  binary_formula = NULL,
  model = c("pooling", "random"),
  link = c("probit", "logit"),
  positive = c("gamma", "lognormal"),
  retransform = c("smearing", "normal"),
  correlated = FALSE,
  cre = FALSE,
  weights = NULL,
  vcov = c("classical", "HC1", "cluster"),
  R = 21,
  maxit = 100,
  tol = 1e-10
)
```

## Arguments

- formula:

  model formula for the positive part, e.g. `cost ~ age + treatment`;
  also used for the binary part unless `binary_formula` is given

- data:

  a data.frame

- index:

  length-1 or length-2 character vector of panel index columns

- binary_formula:

  optional right-hand-side formula for the binary part, e.g.
  `~ age + treatment + smoker`

- model:

  `"pooling"` or `"random"`, used for both parts

- link:

  link of the binary part, `"probit"` or `"logit"`

- positive:

  `"gamma"` (gamma GLM, log link) or `"lognormal"` (linear model for
  `log(y)`)

- retransform:

  for `positive = "lognormal"`: `"smearing"` or `"normal"`

- correlated:

  logical; with `model = "random"` and `positive = "lognormal"`, let the
  random intercepts of the two parts be correlated and fit them jointly

- cre:

  logical; add individual means of the time-varying regressors to both
  parts

- weights:

  optional observation weights (`model = "pooling"` only): the name of a
  column of `data` or a numeric vector with one value per row of `data`

- vcov:

  covariance type passed to
  [`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md) for
  both parts

- R, maxit, tol:

  passed to
  [`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md)

## Value

an object of class `"twopart"` with components `binary` and `positive`
(the two
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md) fits),
`retransformation` (the factor `phi`, 1 for the gamma part) and the
settings used. Methods: `print`, `summary`, `coef`, `vcov`, `confint`,
`predict`, `fitted`, `logLik`, `nobs`,
[`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md),
`tidy` and `glance`.

## Details

The positive part is either a gamma GLM with a log link
(`positive = "gamma"`, no retransformation needed) or a linear model for
`log(y)` (`positive = "lognormal"`). For the log model,
`E[y | y > 0] = exp(x'beta) * phi` needs a retransformation factor
`phi`: Duan's smearing estimate, the (weighted) mean of the
exponentiated log-scale residuals (`retransform = "smearing"`), or the
normal-theory factor `exp(sigma^2 / 2)` (`retransform = "normal"`; with
`model = "random"`, `sigma^2` is the sum of the individual and
idiosyncratic variance components). Both assume homoskedastic log-scale
errors.

With `model = "random"` each part has its own random intercept and the
parts are estimated separately; the combined mean multiplies the two
population-averaged means, which assumes the two random intercepts are
independent. With `correlated = TRUE` (random effects, lognormal
positive part) the two random intercepts are instead bivariate normal
with correlation `rho` and both parts are estimated jointly by maximum
likelihood (Van Ourti 2004): individuals with a high propensity to use
care may also have high spending when they do. The positive-part random
intercept is integrated out in closed form and the binary-part one by
adaptive Gauss-Hermite quadrature. Predictions then integrate over the
joint distribution (normal-theory retransformation). A likelihood-ratio
test of `rho = 0` compares
[`logLik()`](https://rdrr.io/r/stats/logLik.html) with the separate fit.
`cre = TRUE` adds the individual means of the time-varying regressors
(computed over all observations of the individual, including zeros) to
both parts, the Mundlak correlated random effects device.

## Validation

Both parts are
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md) fits
and share its validation. The package tests also check that the pooled
two-part prediction is the product of the two
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) fits, the smearing
retransformation of the lognormal part, and the correlated
random-intercept likelihood against two-dimensional numerical
integration.

## References

Duan, N., Manning, W. G., Morris, C. N. and Newhouse, J. P. (1983). A
comparison of alternative models for the demand for medical care.
*Journal of Business and Economic Statistics*, 1, 115-126.
[doi:10.1080/07350015.1983.10509330](https://doi.org/10.1080/07350015.1983.10509330)

Duan, N. (1983). Smearing estimate: a nonparametric retransformation
method. *Journal of the American Statistical Association*, 78, 605-610.
[doi:10.1080/01621459.1983.10478017](https://doi.org/10.1080/01621459.1983.10478017)

Van Ourti, T. (2004). Measuring horizontal inequity in Belgian health
care using a Gaussian random effects two part count data model. *Health
Economics*, 13, 705-724.
[doi:10.1002/hec.869](https://doi.org/10.1002/hec.869)

## See also

[`predict.twopart()`](https://elbadisyimad.com/panglm/reference/predict.twopart.md),
[`confint.twopart()`](https://elbadisyimad.com/panglm/reference/confint.twopart.md)
and
[`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md)
for the fitted model,
[tidy-methods](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
for `tidy()` and `glance()`.

Other panel models:
[`panfmm()`](https://elbadisyimad.com/panglm/reference/panfmm.md),
[`panfmmselect()`](https://elbadisyimad.com/panglm/reference/panfmmselect.md),
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md),
[`pangmm()`](https://elbadisyimad.com/panglm/reference/pangmm.md),
[`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md),
[`paniv()`](https://elbadisyimad.com/panglm/reference/paniv.md),
[`panqr()`](https://elbadisyimad.com/panglm/reference/panqr.md),
[`panselect()`](https://elbadisyimad.com/panglm/reference/panselect.md)

## Examples

``` r
data(copd)
set.seed(1)
copd$cost <- ifelse(copd$exacerbations > 0,
                    exp(6 + 0.1 * copd$crp + rnorm(nrow(copd), 0, 0.5)), 0)
tp <- twopart(cost ~ treatment + age + crp, data = copd,
              index = c("id", "visit"), positive = "gamma")
summary(tp)
#> 
#> Call:
#> twopart(formula = cost ~ treatment + age + crp, data = copd, 
#>     index = c("id", "visit"), positive = "gamma")
#> 
#> Estimator: two-part model (pooling)
#> N = 300  groups = 60  positive = 238  vcov: classical 
#> 
#> Part 1: P(y > 0), binomial probit 
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept) -0.99307    0.72303  -1.373  0.16960   
#> treatment   -0.11262    0.17573  -0.641  0.52159   
#> age          0.02725    0.01002   2.721  0.00651 **
#> crp          0.02321    0.04127   0.563  0.57374   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Part 2: E[y | y > 0], gamma GLM with log link 
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) 6.1231335  0.2905320  21.076  < 2e-16 ***
#> treatment   0.0480158  0.0755347   0.636    0.525    
#> age         0.0003397  0.0039080   0.087    0.931    
#> crp         0.0937723  0.0171136   5.479 4.27e-08 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Log-likelihood: -1882 
#> 
head(predict(tp))
#> [1] 704.2566 390.6360 586.9989 649.5538 541.2372 557.9110
partialeffects(tp)
#> 
#> Average partial effects (total mean)
#> Delta-method standard errors, vcov: classical 
#> 
#>       term contrast estimate std.error statistic   p.value  conf.low conf.high
#>  treatment    1 - 0    5.320    61.550   0.08642 9.311e-01 -115.3000    126.00
#>        age    dY/dX    6.135     3.206   1.91400 5.567e-02   -0.1486     12.42
#>        crp    dY/dX   63.770    14.670   4.34600 1.386e-05   35.0100     92.54
#> 
```
