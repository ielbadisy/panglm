# Two-way fixed-effects difference-in-differences

Adds the indicator `treated = 1(t >= first treated period)` to the model
and fits it with individual and time fixed effects (two-way fixed
effects, TWFE). With a single treatment date (the canonical two-group,
two-period design and its multi-period extension) the coefficient on
`treated` is the average treatment effect on the treated under parallel
trends. With staggered adoption and effects that vary over time or
across cohorts the TWFE coefficient is a weighted average that can put
negative weight on some effects (Goodman-Bacon 2021; de Chaisemartin and
D'Haultfoeuille 2020): use
[`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md)
to look at the dynamics and a heterogeneity-robust estimator
([`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md) for
Callaway and Sant'Anna 2021, or
[`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md) for the
extended TWFE of Wooldridge 2021, 2023) when adoption is staggered.

## Usage

``` r
did(formula, data, index, treat, family = "gaussian", vcov = "cluster", ...)
```

## Arguments

- formula:

  outcome and control variables, e.g. `y ~ age` or `y ~ 1`; the
  treatment indicator is added

- data:

  a data.frame in long format

- index:

  `c(individual, time)`; the time index must be numeric

- treat:

  name of the column holding each individual's first treated period (in
  the units of the time index); `NA`, `0` or `Inf` for units that are
  never treated

- family:

  `"gaussian"`, `"poisson"`, `"negbin"` or `"binomial"` (conditional
  logit, logit link only; see Details)

- vcov:

  covariance type, cluster-robust by individual by default

- ...:

  further arguments passed to
  [`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md)

## Value

a [`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md) fit
whose coefficient `treated` is the DiD estimate (a log odds ratio for
`family = "binomial"`), with a component `did` (number of treated and
never-treated individuals, treatment cohorts). Its `call` is the `did()`
call. All
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md)
methods apply.

## Details

The family can be Gaussian (effect on the mean), Poisson or negative
binomial (effect on the log rate: `exp(coef)` is a ratio, and Poisson FE
DiD needs only the conditional mean to be correct, Wooldridge 2023), or
binomial for a binary outcome. The binomial model is the two-way
conditional (Chamberlain) logit of
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md): the
individual effects are conditioned out exactly and the time effects
enter as dummies inside the conditional likelihood, as in
`survival::clogit(y ~ treated + factor(time) + strata(id))`. Its
coefficient `treated` is then a log odds ratio, conditional on the
individual effect (`exp(coef)` is the odds ratio), parallel trends is
assumed on the log-odds scale, and individuals whose outcome never
changes carry no information and are dropped. It is not an effect on the
probability scale: for average effects on the treated in probability
points use
[`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md) with
`family = "binomial"`. Controls that do not vary within individuals are
absorbed by the individual effects and dropped.

## Validation

Numerical accuracy is checked in the package tests and in the validation
suite; see the validation article on the package website.

## References

Goodman-Bacon, A. (2021). Difference-in-differences with variation in
treatment timing. *Journal of Econometrics*, 225, 254-277.
[doi:10.1016/j.jeconom.2021.03.014](https://doi.org/10.1016/j.jeconom.2021.03.014)

de Chaisemartin, C. and D'Haultfoeuille, X. (2020). Two-way fixed
effects estimators with heterogeneous treatment effects. *American
Economic Review*, 110, 2964-2996.
[doi:10.1257/aer.20181169](https://doi.org/10.1257/aer.20181169)

Wooldridge, J. M. (2023). Simple approaches to nonlinear
difference-in-differences with panel data. *The Econometrics Journal*,
26, C31-C66.
[doi:10.1093/ectj/utad016](https://doi.org/10.1093/ectj/utad016)

## See also

[`summary.panglm()`](https://elbadisyimad.com/panglm/reference/summary.panglm.md)
and the other
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md)
methods for the fitted model.

Other policy evaluation:
[`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md),
[`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md),
[`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md)

## Examples

``` r
set.seed(1)
d <- expand.grid(id = 1:200, year = 2015:2022)
d$first <- ifelse(d$id <= 100, 2019, NA)        # policy in 2019 for half
a <- rnorm(200)
d$visits <- 3 + a[d$id] + 0.1 * (d$year - 2015) +
  1.5 * (!is.na(d$first) & d$year >= d$first) + rnorm(nrow(d))
fit <- did(visits ~ 1, data = d, index = c("id", "year"), treat = "first")
summary(fit)
#> 
#> Call:
#> did(formula = visits ~ 1, data = d, index = c("id", "year"), 
#>     treat = "first")
#> 
#> Family: gaussian  Link: identity  Model: within  Effect: twoways 
#> N = 1600  groups = 200  vcov: cluster 
#> 
#> Coefficients:
#>         Estimate Std. Error z value Pr(>|z|)    
#> treated  1.71437    0.09866   17.38   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -2221 
#> Dispersion parameter: 1.081 
#> Residual df: 1392  Iterations: 2
#> 
# binary outcome: conditional logit, the coefficient is a log odds ratio
d$admitted <- rbinom(nrow(d), 1, plogis(-1 + a[d$id] + 0.05 * (d$year - 2015) +
  0.8 * (!is.na(d$first) & d$year >= d$first)))
fitb <- did(admitted ~ 1, data = d, index = c("id", "year"), treat = "first",
            family = "binomial")
#> 16 group(s) with no within-group outcome variation dropped (all-0 or all-1) -- they carry no information for the conditional likelihood.
exp(coef(fitb))                                  # odds ratio
#>  treated 
#> 2.020067 
```
