
# panglm <img src="man/figures/logo.svg" align="right" height="139" alt="panglm logo" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/panglm)](https://CRAN.R-project.org/package=panglm)
[![R-CMD-check](https://github.com/ielbadisy/panglm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ielbadisy/panglm/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

panglm fits regression models for panel data (pooled, fixed effects,
random effects and the estimators built on them) for the outcomes and
designs of applied health econometrics, from one formula interface and
one set of post-estimation tools (robust and cluster covariances,
partial effects, specification tests, prediction, `tidy()`/`glance()`
output and base-graphics plots for every fitted object).

## Installation

``` r
# release from CRAN
install.packages("panglm")

# development version from GitHub
# install.packages("remotes")
remotes::install_github("ielbadisy/panglm")
```

## What it fits

| Task                                                                                                                  | Function                                                                                       |
|-----------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------|
| Continuous, binary, count, ordered, positive, bounded and censored outcomes (pooled, within, between, random effects) | `panglm()`                                                                                     |
| Random slopes                                                                                                         | `panglm(model = "random", random = ~ x)`                                                       |
| Costs with zeros (two-part model)                                                                                     | `twopart()`                                                                                    |
| Counts with structural zeros                                                                                          | `panhurdle()`                                                                                  |
| Difference-in-differences and event studies                                                                           | `did()`, `eventstudy()`                                                                        |
| Staggered adoption                                                                                                    | `attgt()`, `etwfe()`                                                                           |
| Instrumental variables and Hausman-Taylor                                                                             | `paniv()`                                                                                      |
| Dynamic panels (difference and system GMM)                                                                            | `pangmm()`                                                                                     |
| Sample selection                                                                                                      | `panselect()`                                                                                  |
| Quantile regression                                                                                                   | `panqr()`                                                                                      |
| Latent classes                                                                                                        | `panfmm()`, `panfmmselect()`                                                                   |
| Specification tests                                                                                                   | `panhausman()`, `pandispersiontest()`, `jointtest()`, `lineartest()`                           |
| Marginal effects and individual effects                                                                               | `partialeffects()`, `effectcurve()`, `fixef()`, `ranef()`                                      |
| Plots                                                                                                                 | `plot()` methods, `coefplot()`, `eventplot()`, `trendplot()`, `treatmentplot()`, `panelplot()` |

Families: `gaussianfamily()`, `binomialfamily()`, `poissonfamily()`,
`negbinfamily()`, `orderedfamily()`, `gammafamily()`, `betafamily()` and
`tobitfamily()`. Names with underscores (`poisson_family()`,
`panglm_hausman()`, ...) are deprecated; see `?"panglm-deprecated"`.

Numerical accuracy is documented in the [validation
article](https://elbadisyimad.com/panglm/articles/validation.html).

## Example

Exacerbation counts of 60 COPD patients over five visits, with a
random-effects Poisson model, cluster-robust standard errors, a Hausman
test against fixed effects and average partial effects:

``` r
library(panglm)
data(copd)

re <- panglm(exacerbations ~ crp + treatment + age, data = copd,
             index = c("id", "visit"), model = "random",
             family = "poisson", vcov = "cluster")
summary(re)
#> 
#> Call:
#> panglm(formula = exacerbations ~ crp + treatment + age, data = copd, 
#>     index = c("id", "visit"), model = "random", family = "poisson", 
#>     vcov = "cluster")
#> 
#> Family: poisson  Link: log  Model: random  Effect: individual 
#> n = 300  groups = 60  vcov: cluster 
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -1.840104   0.525745  -3.500 0.000465 ***
#> crp          0.040731   0.020385   1.998 0.045711 *  
#> treatment   -0.309009   0.168576  -1.833 0.066795 .  
#> age          0.037769   0.007959   4.745 2.08e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Log-likelihood: -517.8 
#> Dispersion parameter: 4.215 
#> Residual df: 295  Iterations: 4

fe <- panglm(exacerbations ~ crp, data = copd,
             index = c("id", "visit"), model = "within", family = "poisson")
panhausman(fe, re)
#> 
#>  Hausman specification test (fixed vs. random effects)
#> 
#> chisq = 3.211, df = 1, p-value = 0.07312
#> alternative hypothesis: one model is inconsistent

partialeffects(re)
#> 
#> Average partial effects
#> Delta-method standard errors, vcov: cluster 
#> 
#>       term contrast estimate std.error statistic   p.value  conf.low conf.high
#>        crp    dY/dX  0.08660   0.04280     2.023 0.0430500  0.002708   0.17050
#>  treatment    1 - 0 -0.63650   0.34350    -1.853 0.0638500 -1.310000   0.03667
#>        age    dY/dX  0.08031   0.01837     4.372 0.0000123  0.044310   0.11630
```

Every fitted object has a `plot()` method. Here, a hanging rootogram
checks the count distribution, and `effectcurve()` traces the average
predicted count as CRP varies:

``` r
op <- par(mfrow = c(1, 2))
plot(re, which = "rootogram")
plot(effectcurve(re, "crp", scale = "response"))
```

<img src="man/figures/README-plots-1.png" alt="" width="100%" />

``` r
par(op)
```

## Datasets

`copd` is a simulated panel with known parameters (counts, costs with
zeros, an ordered score, a proportion and a staggered care programme).
Five real health economics panels ship for teaching: `docvisits` (German
health care demand), `meps` (Medical Expenditure Panel Survey spending),
`medicaid` (US states, Medicaid expansion), `healthspend` (countries,
health expenditure) and `cigar` (US states, cigarette demand). Each help
page gives the source, the licence and the questions the data can teach.

## Learn more

- [Getting started with
  panglm](https://elbadisyimad.com/panglm/articles/panglm-intro.html)
- [Counts, costs, bounded and censored
  outcomes](https://elbadisyimad.com/panglm/articles/panglm-outcomes.html)
- [Endogeneity: IV, Hausman-Taylor, GMM and
  selection](https://elbadisyimad.com/panglm/articles/panglm-endogeneity.html)
- [Policy evaluation: DiD, event studies, attgt and
  etwfe](https://elbadisyimad.com/panglm/articles/panglm-policy.html)
- [Heterogeneity: random slopes, latent classes and quantile
  regression](https://elbadisyimad.com/panglm/articles/panglm-heterogeneity.html)
- [Teaching
  datasets](https://elbadisyimad.com/panglm/articles/panglm-datasets.html)
- [Technical
  reference](https://elbadisyimad.com/panglm/articles/panglm-methods.html)
- [Function reference](https://elbadisyimad.com/panglm/reference/)

panglm is the software of a forthcoming book on panel data methods for
health econometrics, whose examples it runs throughout.
