# Instrumental-variable estimators for linear panel models

Two-stage least squares for linear panel models with endogenous
regressors: fixed-effects 2SLS, random-effects generalized 2SLS
(Balestra and Varadharajan-Krishnakumar) or error-components 2SLS
(Baltagi), pooled 2SLS, and the Hausman and Taylor (1981) estimator,
which also identifies the effects of time-invariant regressors.

## Usage

``` r
paniv(
  formula,
  data,
  index,
  model = c("within", "random", "pooling", "hausmantaylor"),
  effect = c("individual", "twoways"),
  inst.method = c("bvk", "baltagi"),
  vcov = c("classical", "cluster"),
  na.action = stats::na.omit
)
```

## Arguments

- formula:

  a two-part formula `y ~ regressors | instruments`; see the Formula
  section. `lag(x, k)`, `lead(x, k)` and `diff(x)` are taken within each
  individual on the grid of sorted distinct time values (needs
  `index = c(id, time)`); the first waves of each individual get `NA`
  and leave the sample, and `lag(x, 0:2)` gives several lags

- data:

  a data.frame in long format

- index:

  length-1 or length-2 character vector giving the individual (and time)
  columns

- model:

  one of `"within"`, `"random"`, `"pooling"` or `"hausmantaylor"`

- effect:

  `"individual"`, or `"twoways"` (only for `model = "within"`)

- inst.method:

  instrument set for `model = "random"`: `"bvk"` (G2SLS, the default, as
  in `plm`) or `"baltagi"` (EC2SLS)

- vcov:

  `"classical"` or `"cluster"` (cluster-robust by individual); see
  [`vcov.paniv()`](https://ielbadisy.github.io/panglm/reference/vcov.paniv.md)
  to recompute after fitting

- na.action:

  missing-data action passed to
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)

## Value

an object of class `"paniv"`, a list with components `coefficients`,
`vcov`, `vcov_type`, `residuals` (structural residuals `y - X b` on the
original scale, in the panel order), `model`, `effect`, `inst.method`,
`endogenous` (names of the endogenous regressors), `instruments` (names
of the excluded instruments, or for Hausman-Taylor the four variable
groups in `varlist`), `first_stage` (a data frame with one row per
endogenous regressor, or `NULL`), `overid` (the Sargan or Hansen test,
or `NULL` when exactly identified), `ercomp` (variance components and
`theta` for the random and Hausman-Taylor estimators), `nobs`,
`n_groups`, `df.residual`, `call` and `formula`. Methods: `print`,
`summary`, `coef`, `vcov`, `confint`, `predict`, `nobs`, `tidy` and
`glance`.

## Formula

The formula has two parts separated by `|`. For `model = "pooling"`,
`"within"` and `"random"` the second part lists the complete instrument
set, as in [`plm::plm()`](https://rdrr.io/pkg/plm/man/plm.html):
exogenous regressors are repeated on both sides and excluded instruments
appear only on the right, so that `y ~ x1 + x2 | z1 + x2` treats `x1` as
endogenous, instrumented by `z1`. Regressors absent from the second part
are the endogenous ones.

For `model = "hausmantaylor"` the conventions of
[`plm::pht()`](https://rdrr.io/pkg/plm/man/pht.html) are used: the first
part lists all regressors and the second part lists those that are
exogenous (uncorrelated with the individual effect), e.g.
`y ~ x1 + x2 + z1 + z2 | x1 + z1` makes `x2` and `z2` endogenous.
Whether a regressor is time varying or time invariant is detected from
the data. Identification requires at least as many exogenous
time-varying regressors as endogenous time-invariant ones.

## Estimators

- `"within"`: 2SLS on within-demeaned regressors and instruments
  (one-way or, with `effect = "twoways"`, two-way demeaning). Time
  invariant columns are dropped. The classical covariance uses the
  degrees of freedom `n - k - N` (`N` individual effects), as in `plm`.

- `"random"`: the variance components are estimated in the Swamy-Arora
  way from within and between 2SLS residuals (the unbalanced-panel
  degrees of freedom of `plm` are used, which reduce to the usual ones
  for balanced panels). The data are then quasi-demeaned with
  `theta_i = 1 - sqrt(s2_nu / (s2_nu + T_i s2_mu))`. With
  `inst.method = "bvk"` (the default, G2SLS) the instruments are
  quasi-demeaned as well; with `inst.method = "baltagi"` (EC2SLS) the
  instrument matrix is the within and between transforms of the
  instruments.

- `"pooling"`: ordinary 2SLS on the stacked data.

- `"hausmantaylor"`: the variance components follow
  [`plm::pht()`](https://rdrr.io/pkg/plm/man/pht.html): `s2_nu` is the
  within residual variance, and `s2_1` comes from a 2SLS regression of
  the estimated individual effects on the time-invariant regressors,
  instrumented by the exogenous time-invariant and time-varying
  regressors. The GLS-transformed regression is then estimated by 2SLS
  with instruments: within transforms of all time-varying regressors,
  individual means of the exogenous time-varying regressors, the
  exogenous time-invariant regressors and a constant. For unbalanced
  panels `s2_mu = (s2_1 - s2_nu) / T_h` with `T_h` the harmonic mean of
  the `T_i`, as in [`plm::pht()`](https://rdrr.io/pkg/plm/man/pht.html).
  A negative `s2_mu` is set to zero (with a warning), which reduces the
  estimator to pooled 2SLS on the same instruments.

In every case the classical covariance is `s2 inverse(Xhat'Xhat)` with
`s2` the mean square of the structural residuals of the transformed
regression. The cluster covariance is the sandwich
`B (sum_g u_g u_g') B` with `u_g = Xhat_g' e_g` summed within the panel
individual (or the supplied cluster) and the factor
`G / (G - 1) * (n - 1) / (n - k)`, the convention of the rest of the
package and of `plm::vcovHC(type = "sss")`.

## Diagnostics

For each endogenous regressor, the first-stage regression of its
transformed values on the transformed instruments is compared with the
regression on the transformed included exogenous regressors only. The
result reports the partial R-squared of the excluded instruments and the
F statistic for their joint significance (with the absorbed fixed
effects counted in the residual degrees of freedom); for
`vcov = "cluster"` a cluster-robust Wald F is reported as well. These
first-stage statistics are not computed for `"hausmantaylor"`, whose
instruments are built from the regressors themselves. When the model is
overidentified, the Sargan statistic `n e'P_W e / e'e` (classical) or
the cluster-robust Hansen J statistic of the two-step efficient GMM
estimator (cluster) tests the validity of the overidentifying
restrictions; it is chi-squared with `L - k` degrees of freedom under
the null.

## Validation

Checked in the package's validation suite on GitHub: the within, random
and pooled 2SLS estimators against
[`plm::plm()`](https://rdrr.io/pkg/plm/man/plm.html) (same
`inst.method`), [`AER::ivreg()`](https://rdrr.io/pkg/AER/man/ivreg.html)
and
[`fixest::feols()`](https://lrberge.github.io/fixest/reference/feols.html)
on [`plm::Crime`](https://rdrr.io/pkg/plm/man/Crime.html), balanced and
unbalanced; the Hausman-Taylor estimator against
[`plm::pht()`](https://rdrr.io/pkg/plm/man/pht.html) and
`plm::plm(random.method = "ht")` on
[`plm::Wages`](https://rdrr.io/pkg/plm/man/Wages.html) (Baltagi 2013,
Table 7.4).

## References

Baltagi, B. H. (1981). Simultaneous equations with error components.
*Journal of Econometrics*, 17, 189-200.

Balestra, P. and Varadharajan-Krishnakumar, J. (1987). Full information
estimations of a system of simultaneous equations with error component
structure. *Econometric Theory*, 3, 223-246.

Hausman, J. A. and Taylor, W. E. (1981). Panel data and unobservable
individual effects. *Econometrica*, 49, 1377-1398.

## See also

[`summary.paniv()`](https://ielbadisy.github.io/panglm/reference/summary.paniv.md),
[`vcov.paniv()`](https://ielbadisy.github.io/panglm/reference/vcov.paniv.md)
and
[`predict.paniv()`](https://ielbadisy.github.io/panglm/reference/predict.paniv.md)
for the fitted model,
[tidy-methods](https://ielbadisy.github.io/panglm/reference/tidy-methods.md)
for `tidy()` and `glance()`.

Other panel models:
[`panfmm()`](https://ielbadisy.github.io/panglm/reference/panfmm.md),
[`panfmmselect()`](https://ielbadisy.github.io/panglm/reference/panfmmselect.md),
[`panglm()`](https://ielbadisy.github.io/panglm/reference/panglm.md),
[`pangmm()`](https://ielbadisy.github.io/panglm/reference/pangmm.md),
[`panhurdle()`](https://ielbadisy.github.io/panglm/reference/panhurdle.md),
[`panqr()`](https://ielbadisy.github.io/panglm/reference/panqr.md),
[`panselect()`](https://ielbadisy.github.io/panglm/reference/panselect.md),
[`twopart()`](https://ielbadisy.github.io/panglm/reference/twopart.md)

## Examples

``` r
# Simulated panel of 500 patients over 4 waves. Private insurance is
# endogenous for doctor visits (a time-varying health shock drives both)
# and is instrumented by an employer offer of insurance.
set.seed(1)
N <- 500; TT <- 4
d <- data.frame(id = rep(seq_len(N), each = TT), wave = rep(seq_len(TT), N))
alpha <- rep(rnorm(N), each = TT)       # unobserved frailty
s0 <- rep(rnorm(N), each = TT)          # family background
need <- rnorm(N * TT)                   # time-varying health shock
d$school <- round(12 + 2 * s0 + alpha)  # correlated with the frailty
d$female <- rep(rbinom(N, 1, 0.5), each = TT)
d$offer <- rbinom(N * TT, 1, plogis(s0))
d$income <- rnorm(N * TT) + 0.5 * alpha
d$insured <- as.numeric(d$offer + need + rnorm(N * TT) > 0.5)
d$visits <- 1 + 0.6 * d$insured + 0.2 * d$income + alpha + need + rnorm(N * TT)

# the within estimator is biased upwards; FE-2SLS is not (true value 0.6)
coef(panglm(visits ~ insured + income, data = d, index = c("id", "wave"),
            model = "within"))
#>  insured   income 
#> 1.697161 0.175453 
fe2sls <- paniv(visits ~ insured + income | offer + income, data = d,
                index = c("id", "wave"), model = "within", vcov = "cluster")
summary(fe2sls)
#> 
#> Call:
#> paniv(formula = visits ~ insured + income | offer + income, data = d, 
#>     index = c("id", "wave"), model = "within", vcov = "cluster")
#> 
#> Estimator: fixed-effects 2SLS 
#> N = 2000  groups = 500  vcov: cluster 
#> Endogenous: insured 
#> Excluded instruments: offer 
#> 
#> Coefficients:
#>         Estimate Std. Error z value Pr(>|z|)    
#> insured  0.86178    0.36261   2.377   0.0175 *  
#> income   0.17291    0.03735   4.629 3.67e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> First stage (excluded instruments):
#>         partial.R2     F df1  df2  p.value F.robust  p.robust
#> insured    0.04144 64.76   1 1498 1.71e-15    66.21 3.252e-15
#> Residual df: 1498 
#> 

# Hausman-Taylor: the effect of schooling (time invariant, correlated
# with the frailty) on health costs, identified by the exogenous
# time-varying employer offer; income is endogenous (true values 0.4,
# 0.2, 0.3, 0.1)
d$cost <- 2 + 0.4 * d$offer + 0.2 * d$income + 0.3 * d$female +
  0.1 * d$school + alpha + rnorm(N * TT)
ht <- paniv(cost ~ offer + income + female + school | offer + female,
            data = d, index = c("id", "wave"), model = "hausmantaylor")
summary(ht)
#> 
#> Call:
#> paniv(formula = cost ~ offer + income + female + school | offer + 
#>     female, data = d, index = c("id", "wave"), model = "hausmantaylor")
#> 
#> Estimator: Hausman-Taylor 
#> N = 2000  groups = 500  vcov: classical 
#> Exogenous time-varying:   offer 
#> Endogenous time-varying:  income 
#> Exogenous time-invariant: female 
#> Endogenous time-invariant: school 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  2.29707    0.43720   5.254 1.49e-07 ***
#> offer        0.47878    0.05873   8.152 3.59e-16 ***
#> income       0.24645    0.02615   9.425  < 2e-16 ***
#> female       0.24550    0.09953   2.467   0.0136 *  
#> school       0.07426    0.03699   2.008   0.0447 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Variance components: idiosyncratic = 0.9919, individual = 0.9807
#> theta: 0.5508 
#> Residual df: 1995 
#> 
```
