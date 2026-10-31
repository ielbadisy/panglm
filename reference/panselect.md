# Sample-selection correction for linear panel models

Two-step selection correction of Wooldridge (1995) for linear panel
models in which the outcome is observed only for a selected subsample
(attrition in a cohort, survival, response to a wave), when selection
may depend on the idiosyncratic error of the outcome equation.

## Usage

``` r
panselect(
  selection,
  outcome,
  data,
  index,
  model = c("within", "pooling"),
  mundlak = TRUE,
  lambda = c("period", "common"),
  boot = 199L,
  seed = NULL,
  na.action = stats::na.omit
)
```

## Arguments

- selection:

  a two-sided formula `s ~ z1 + x1 + ...` for the selection equation;
  `s` is a 0/1 (or logical) indicator equal to one when the outcome is
  observed

- outcome:

  a two-sided formula `y ~ x1 + ...` for the outcome equation, fitted on
  the observations with `s == 1`

- data:

  a data.frame in long format, with one row per individual and period
  whether or not the outcome is observed

- index:

  length-2 character vector giving the individual and period columns

- model:

  `"within"` (fixed effects in the outcome equation) or `"pooling"`
  (pooled OLS, with Mundlak means when `mundlak = TRUE`)

- mundlak:

  logical; add the individual means of the time-varying selection
  regressors to the probits (and to the pooled outcome equation)

- lambda:

  `"period"` (inverse Mills ratio interacted with period dummies) or
  `"common"` (a single inverse Mills ratio term)

- boot:

  number of panel bootstrap replicates; `0` skips the bootstrap and
  makes the cluster-robust covariance the default

- seed:

  optional random seed for the bootstrap; the global random number state
  is restored afterwards

- na.action:

  missing-data action passed to
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html) for
  each formula

## Value

an object of class `"panselect"`, a list with components `coefficients`,
`vcov` (the default covariance), `vcov_type` (`"boot"` or `"cluster"`),
`vcov_boot`, `vcov_cluster`, `vcov_classical`, `boot` (the replicate
estimates, the number of requested and failed replicates), `probit` (a
list with the probit coefficients of each period, `NULL` for a period
without probit), `lambda_terms` (names of the correction terms),
`selection_test` (the cluster-robust and bootstrap Wald tests), `imr`
(the inverse Mills ratio of every observation of the selection sample,
in the row order of the data), `residuals`, `nobs` (selected
observations used), `nobs_selection`, `n_groups`, `n_groups_selection`,
`df.residual`, `periods`, `model`, `mundlak`, `lambda`, `call`. Methods:
`print`, `summary`, `coef`, `vcov`, `confint`, `predict`, `nobs`, `tidy`
and `glance`.

## Estimator

Step 1. For each period `t` a probit of the selection indicator `s_it`
on the selection regressors `z_it` and, with `mundlak = TRUE`, the
individual means `zbar_i` of the time-varying selection regressors (the
Mundlak version of Chamberlain's device) is fitted on all observations
of that period (by iteratively reweighted least squares, as
[`stats::glm()`](https://rdrr.io/r/stats/glm.html), with the convergence
tolerance `1e-12`). The inverse Mills ratio is
`lambda_it = dnorm(xb_it) / pnorm(xb_it)` with `xb_it` the probit index.
The means are taken over all observations of the individual in the
selection sample, selected or not. A period in which every observation
(or none) is selected has no probit; its inverse Mills ratio is set to
zero and it gets no correction term (for instance the first wave of a
cohort, before any attrition).

Step 2. The outcome equation is fitted on the selected observations,
adding the inverse Mills ratio interacted with period dummies,
`lambda_it * d_t` (`lambda = "period"`, Wooldridge's specification, one
coefficient per period with a probit) or a single term `lambda_it`
(`lambda = "common"`). With `model = "within"` the individual effects
are removed by the within transformation on the selected sample; with
`model = "pooling"` the outcome equation is estimated by pooled OLS and,
when `mundlak = TRUE`, includes the same Mundlak means `zbar_i` to
account for the correlation between the individual effect and the
regressors. Period effects are not added automatically; include them in
`outcome` (e.g. `+ factor(wave)`) when needed.

The selection regressors should contain the outcome regressors plus at
least one excluded variable (an instrument that affects selection but
not the outcome); without one the correction is identified only by the
nonlinearity of the inverse Mills ratio, and a warning is given.

## Inference

The coefficient standard errors default to a panel (cluster) bootstrap:
whole individuals are resampled with replacement and both steps,
including the Mundlak means and the per-period probits, are re-estimated
in each replicate. The covariance is the sample covariance of the `boot`
replicate estimates; replicates that fail are dropped and counted.

The cluster-robust covariance `vcov(fit, type = "cluster")` treats the
inverse Mills ratios as known regressors. It ignores the estimation of
the first step and is therefore not valid for the coefficients when
there is selection bias. It uses the factor
`G / (G - 1) * (n - 1) / (n - k)`, as in `plm::vcovHC(type = "sss")`.

The test of selection bias is the joint Wald test that all inverse Mills
ratio coefficients are zero. Under this null hypothesis the first step
does not affect the asymptotic distribution of the second-step
estimator, so the cluster-robust Wald statistic is valid (Wooldridge
2010, section 19.9); it is chi-squared with as many degrees of freedom
as correction terms. The same statistic computed with the bootstrap
covariance is reported as well.

## Validation

No other R package implements this panel estimator. Each step, and the
recovery of the true coefficients in simulations with selection on the
idiosyncratic error, are checked in the package tests and in the
validation suite; see the validation article on the package website.

## References

Wooldridge, J. M. (1995). Selection corrections for panel data models
under conditional mean independence assumptions. *Journal of
Econometrics*, 68, 115-132.

Wooldridge, J. M. (2010). *Econometric Analysis of Cross Section and
Panel Data*, 2nd edition, section 19.9. MIT Press.

## See also

[`summary.panselect()`](https://elbadisyimad.com/panglm/reference/summary.panselect.md),
[`vcov.panselect()`](https://elbadisyimad.com/panglm/reference/vcov.panselect.md)
and
[`predict.panselect()`](https://elbadisyimad.com/panglm/reference/predict.panselect.md)
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
[`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md)

## Examples

``` r
# Simulated cohort of 300 patients followed over 4 waves. A quality of
# life score is observed only for patients who answer the questionnaire;
# non-response depends on the same health shock that lowers the score,
# and on a reminder letter that does not affect the score itself.
set.seed(1)
N <- 300; TT <- 4
d <- data.frame(id = rep(seq_len(N), each = TT), wave = rep(seq_len(TT), N))
alpha <- rep(rnorm(N), each = TT)
d$treat <- rnorm(N * TT) + 0.5 * alpha
d$reminder <- rnorm(N * TT)
v <- rnorm(N * TT)
u <- 0.7 * v + rnorm(N * TT, sd = 0.7)
d$respond <- as.numeric(0.3 + 0.5 * d$treat + d$reminder + 0.3 * alpha + v > 0)
d$qol <- ifelse(d$respond == 1, 1 + 0.5 * d$treat + alpha + u, NA)

fit <- panselect(respond ~ treat + reminder, qol ~ treat, data = d,
                 index = c("id", "wave"), model = "within", boot = 49, seed = 2)
summary(fit)
#> 
#> Call:
#> panselect(selection = respond ~ treat + reminder, outcome = qol ~ 
#>     treat, data = d, index = c("id", "wave"), model = "within", 
#>     boot = 49, seed = 2)
#> 
#> Wooldridge (1995) panel selection correction: fixed effects with Mundlak means, inverse Mills ratio by period 
#> Selection sample: N = 1200  groups = 300 
#> Outcome sample:   N = 710  groups = 286 
#> Probits fitted for periods: 1, 2, 3, 4 
#> 
#> Outcome equation (vcov: panel bootstrap, 49 replicates):
#>          Estimate Std. Error z value Pr(>|z|)    
#> treat     0.45225    0.05037   8.978  < 2e-16 ***
#> lambda:1  0.98046    0.18266   5.368 7.98e-08 ***
#> lambda:2  0.76097    0.17850   4.263 2.02e-05 ***
#> lambda:3  0.80668    0.21728   3.713 0.000205 ***
#> lambda:4  0.81105    0.17365   4.671 3.00e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Test of selection bias (H0: all inverse Mills ratio coefficients are zero):
#>   cluster-robust Wald chi2(4) = 56.2, p = 1.818e-11
#>   bootstrap Wald chi2(4) = 44.06, p = 6.232e-09
#> Residual df: 419 
#> 
```
