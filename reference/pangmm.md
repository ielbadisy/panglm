# Difference and system GMM for dynamic linear panel models

Arellano and Bond (1991) difference GMM and Blundell and Bond (1998)
system GMM for linear panel models with a lagged dependent variable and
an unobserved individual effect,
`y_it = sum_j rho_j y_i,t-j + x_it'b + mu_i + e_it`. Both one-step and
two-step estimators are available, with robust one-step and Windmeijer
(2005) corrected two-step standard errors, the Arellano-Bond tests for
first- and second-order serial correlation of the differenced residuals
and the Sargan and Hansen tests of the overidentifying restrictions.

## Usage

``` r
pangmm(
  formula,
  data,
  index,
  transformation = c("d", "ld"),
  model = c("onestep", "twostep"),
  effect = c("twoways", "individual"),
  collapse = FALSE,
  maxlag = NULL,
  vcov = c("robust", "classical", "cluster")
)
```

## Arguments

- formula:

  a formula `y ~ regressors | gmm instruments` or
  `y ~ regressors | gmm instruments | standard instruments`; see the
  Formula section

- data:

  a data.frame in long format

- index:

  length-2 character vector giving the individual and time columns

- transformation:

  `"d"` (difference GMM, Arellano-Bond) or `"ld"` (system GMM,
  Blundell-Bond)

- model:

  `"onestep"` or `"twostep"`

- effect:

  `"twoways"` (time dummies, the default) or `"individual"`

- collapse:

  logical: collapse the GMM-style instruments to one column per lag
  (and, for system GMM, one column for the level equation). A named
  logical vector applies it to selected GMM variables, named as in the
  formula, e.g. `c("log(bmi)" = TRUE)`.

- maxlag:

  optional positive integer: the largest lag of the GMM-style
  instruments, overriding larger upper lags in the formula
  (`lag(y, 2:99)` with `maxlag = 4` uses lags 2 to 4)

- vcov:

  `"robust"` (one-step sandwich or Windmeijer-corrected two-step, both
  clustered by individual) or `"classical"`; `"cluster"` is accepted as
  a synonym of `"robust"`, the name used by the other estimators of the
  package. See
  [`vcov.pangmm()`](https://elbadisyimad.com/panglm/reference/vcov.pangmm.md)
  to switch after fitting

## Value

an object of class `"pangmm"`, a list with components `coefficients`
(the final estimates; time dummies are named by their period),
`coefficients_onestep`, `vcov`, `vcov_type`, `vcov_list` (the classical
and robust covariances), `residuals` (one numeric vector per individual,
differenced rows first and, for system GMM, level rows after; rows not
used are zero), `artest` (data frame of the Arellano-Bond tests),
`sargan`, `hansen` and `wald` (lists with `statistic`, `df` and
`p.value`; `wald` tests the joint significance of the coefficients other
than time dummies and the constant), `sigma2` (`s2`), `n_instruments`,
`nobs` (number of equations used), `nobs_diff`, `nobs_level`,
`n_groups`, `n_periods`, `timedummies` (names of the time-effect
coefficients), `transformation`, `model`, `effect`, `collapse`, `call`
and `formula`. Methods: `print`, `summary`, `coef`, `vcov`, `confint`,
`predict`, `residuals`, `nobs`, `tidy` and `glance`.

## Formula

The formula follows
[`plm::pgmm()`](https://rdrr.io/pkg/plm/man/pgmm.html) and has two or
three parts separated by `|`:

`y ~ lag(y, 1) + x + lag(x, 1) | lag(y, 2:99) | z`

- The first part lists the regressors. `lag(v, k)` is the lag of order
  `k` of `v` within the individual, `k` may be a vector (`lag(v, 0:2)`
  gives `v` and its first two lags) and `lag(v)` means `lag(v, 1)`. Any
  expression of the data can be lagged, e.g. `lag(log(bmi), 1)`. The
  intercept is not estimated: the individual effects are differenced
  out, and time effects or the system-GMM constant are added through
  `effect`.

- The second part lists the GMM-style instruments: `lag(v, a:b)` uses,
  for the differenced equation of period `t`, the levels `v_(t-a)`, ...,
  `v_(t-b)` that exist in the panel (`b = 99` means all available lags).
  For an endogenous variable `a = 2`, for a predetermined one `a = 1`.

- The optional third part lists standard (IV-style) instruments, with
  the same [`lag()`](https://rdrr.io/r/stats/lag.html) syntax; each
  enters in first differences in the differenced equation (and in levels
  in the level equation of system GMM). Without a third part, every
  regressor whose variable does not appear among the GMM-style
  instruments is used as its own standard instrument, as in
  [`plm::pgmm()`](https://rdrr.io/pkg/plm/man/pgmm.html).

Lags are taken on the time grid formed by the sorted distinct values of
the time index, so missing waves are respected and waves need not be
equally spaced in calendar time (`plm` instead lags by the numeric value
of the time index, so biennial waves coded as years have no lag one
period back there); rows with a missing value in the response or a
regressor are dropped from the moment conditions and missing instruments
are set to zero, as in `plm` and Stata's `xtabond2`.

## Estimators

With `transformation = "d"` (difference GMM) the model is estimated in
first differences, with the GMM-style instruments of each period in
their own columns (block-diagonal instrument matrix) or, with
`collapse = TRUE`, one column per lag (Roodman 2009). The one-step
weight matrix is `(sum_i Z_i' H Z_i)^-1` with `H` the tridiagonal matrix
with 2 on the diagonal and -1 off it, the covariance of the differenced
errors under homoskedasticity. With `transformation = "ld"` (system GMM)
the level equations of the same periods are stacked under the
differenced ones and instrumented by the first difference of each
GMM-style variable lagged `a - 1` periods (`a` its smallest lag); the
one-step weight uses the block matrix `H = [G, D; D', I]` of
[`plm::pgmm()`](https://rdrr.io/pkg/plm/man/pgmm.html) (`fsm = "full"`).
The two-step estimator uses the weight `(sum_i Z_i' e_i e_i' Z_i)^-1`
built from the one-step residuals.

With `effect = "twoways"` (the default) time dummies are included: in
differences for difference GMM, and for system GMM a constant and time
dummies in the level equation (which instrument themselves there). With
`effect = "individual"` there are no time effects; system GMM then keeps
a constant in the level equation
([`plm::pgmm()`](https://rdrr.io/pkg/plm/man/pgmm.html) omits it).

Regressors that are not identified, because they are collinear with the
time effects (for instance age in a balanced annual panel, whose first
difference is constant) or with each other, are dropped with a warning,
as in [`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md);
they remain standard instruments when listed as such.

## Inference

`vcov = "robust"` gives the one-step sandwich covariance or, for the
two-step estimator, the Windmeijer (2005) finite-sample corrected
covariance. `vcov = "classical"` gives, for the one-step estimator,
`s2 (X'Z S1^-1 Z'X)^-1` with `S1 = sum_i Z_i' H Z_i` and
`s2 = sum(de^2) / (2 (n_d - k))` estimated from the `n_d` differenced
one-step residuals `de`; for the two-step estimator it gives the
conventional (downward biased) `(X'Z S2^-1 Z'X)^-1` with
`S2 = sum_i Z_i' e1_i e1_i' Z_i`. (The one-step
[`vcov()`](https://rdrr.io/r/stats/vcov.html) of
[`plm::pgmm()`](https://rdrr.io/pkg/plm/man/pgmm.html) is
`(X'Z S1^-1 Z'X)^-1 / N`, without the error variance.)

The Arellano-Bond statistic of order `m` tests for autocorrelation of
order `m` in the differenced residuals; it is standard normal under the
null of no autocorrelation. First-order autocorrelation is expected by
construction; second-order autocorrelation indicates that `y_(t-2)` is
not a valid instrument. It is computed with the robust covariance of the
fit, as in `summary(plm::pgmm())`. The Sargan statistic
`(Z'e1)' (sum Z_i' H Z_i)^-1 (Z'e1) / s2` uses the one-step residuals
and is valid only under homoskedasticity; the Hansen J statistic
`(Z'e2)' (sum Z_i' e1_i e1_i' Z_i)^-1 (Z'e2)` uses the two-step
residuals and is robust. Both are chi-squared with `L - k` degrees of
freedom (`L` instruments, `k` coefficients);
[`plm::sargan()`](https://rdrr.io/pkg/plm/man/sargan.html) returns the
Hansen statistic by default and `N s2` times the Sargan statistic with
`weights = "onestep"`. Many instruments weaken the Hansen test (Roodman
2009); collapsing or limiting the lag depth keeps the instrument count
well below the number of individuals.

## Finite-sample bias

The estimators are consistent as the number of individuals grows, but in
moderate samples they can be badly biased, and matching another package
does not protect against it. Two mechanisms matter. Too many instruments
overfit the endogenous regressors and pull the estimate toward pooled
OLS, while making the Hansen test look too good (Roodman 2009). And when
the variance of the individual effect is large relative to the
idiosyncratic variance, the level moments of system GMM become weak and
the estimate again drifts toward OLS (Bun and Windmeijer 2010). In a
simulated BMI panel (true persistence 0.6, effect SD 1.5, error SD 1, 8
waves, 63 instruments) uncollapsed system GMM averaged 0.78 one-step and
0.72 two-step with 400 individuals, 0.65 and 0.61 with 2000 individuals;
collapsing the instruments gave 0.62 already with 400. Practical rules:
prefer the two-step estimator, use `collapse = TRUE` and a small
`maxlag`, keep the instrument count well below the number of individuals
(`pangmm()` warns when it is not), treat a Hansen p-value close to 1 as
a warning sign rather than as support, and check that the coefficient of
the lagged outcome lies between the fixed-effects (downward biased) and
pooled OLS (upward biased) estimates (Bond 2002).

## Validation

Checked against [`plm::pgmm()`](https://rdrr.io/pkg/plm/man/pgmm.html)
and its `vcovHC()`, `mtest()`, `sargan()` and `pwaldtest()` methods
(difference and system GMM, one and two steps, with the Windmeijer
correction), and [`predict()`](https://rdrr.io/r/stats/predict.html)
against the transformed equations stored by
[`plm::pgmm()`](https://rdrr.io/pkg/plm/man/pgmm.html), on
[`plm::EmplUK`](https://rdrr.io/pkg/plm/man/EmplUK.html) and a simulated
unbalanced panel, in the package's validation suite on GitHub.

## References

Bond, S. (2002). Dynamic panel data models: a guide to micro data
methods and practice. *Portuguese Economic Journal*, 1, 141-162.

Bun, M. J. G. and Windmeijer, F. (2010). The weak instrument problem of
the system GMM estimator in dynamic panel data models. *The Econometrics
Journal*, 13, 95-126.

Arellano, M. and Bond, S. (1991). Some tests of specification for panel
data: Monte Carlo evidence and an application to employment equations.
*Review of Economic Studies*, 58, 277-297.

Blundell, R. and Bond, S. (1998). Initial conditions and moment
restrictions in dynamic panel data models. *Journal of Econometrics*,
87, 115-143.

Windmeijer, F. (2005). A finite sample correction for the variance of
linear efficient two-step GMM estimators. *Journal of Econometrics*,
126, 25-51.

Roodman, D. (2009). How to do xtabond2: An introduction to difference
and system GMM in Stata. *Stata Journal*, 9, 86-136.

## See also

[`summary.pangmm()`](https://elbadisyimad.com/panglm/reference/summary.pangmm.md),
[`vcov.pangmm()`](https://elbadisyimad.com/panglm/reference/vcov.pangmm.md)
and
[`predict.pangmm()`](https://elbadisyimad.com/panglm/reference/predict.pangmm.md)
for the fitted model,
[tidy-methods](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
for `tidy()` and `glance()`.

Other panel models:
[`panfmm()`](https://elbadisyimad.com/panglm/reference/panfmm.md),
[`panfmmselect()`](https://elbadisyimad.com/panglm/reference/panfmmselect.md),
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md),
[`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md),
[`paniv()`](https://elbadisyimad.com/panglm/reference/paniv.md),
[`panqr()`](https://elbadisyimad.com/panglm/reference/panqr.md),
[`panselect()`](https://elbadisyimad.com/panglm/reference/panselect.md),
[`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md)

## Examples

``` r
# Habit persistence in smoking: cigarettes per day depend on last
# wave's consumption (state dependence) and on the price of tobacco,
# with a permanent individual propensity to smoke. The lagged outcome is
# correlated with that propensity, so pooled OLS overstates persistence
# and the within estimator understates it (Nickell bias); GMM is
# consistent. True values: persistence 0.5, price -0.3.
set.seed(1)
N <- 300; TT <- 7
d <- data.frame(id = rep(seq_len(N), each = TT), wave = rep(seq_len(TT), N))
propensity <- rnorm(N)
d$price <- rnorm(N * TT)
d$cigs <- NA
for (i in seq_len(N)) {
  r <- (i - 1) * TT + seq_len(TT)
  y <- propensity[i] / 0.5 + rnorm(1)
  for (t in seq_len(TT)) {
    y <- 0.5 * y - 0.3 * d$price[r[t]] + propensity[i] + rnorm(1)
    d$cigs[r[t]] <- y
  }
}

# Arellano-Bond difference GMM, two-step with Windmeijer standard errors
ab <- pangmm(cigs ~ lag(cigs, 1) + price | lag(cigs, 2:99), data = d,
             index = c("id", "wave"), model = "twostep")
summary(ab)
#> 
#> Call:
#> pangmm(formula = cigs ~ lag(cigs, 1) + price | lag(cigs, 2:99), 
#>     data = d, index = c("id", "wave"), model = "twostep")
#> 
#> Estimator: two-way two-step difference GMM (Arellano-Bond) 
#> N = 300  T = 7  equations used = 1500  instruments = 21  vcov: robust (Windmeijer) 
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> lag(cigs, 1)  0.55025    0.04883  11.268   <2e-16 ***
#> price        -0.27276    0.03055  -8.927   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Arellano-Bond AR(1) test: z = -10.24, p = < 2.2e-16
#> Arellano-Bond AR(2) test: z = -0.9451, p = 0.3446
#> Sargan test (not robust): chi2(14) = 7.391, p = 0.9186
#> Hansen J test (robust): chi2(14) = 6.259, p = 0.9595
#> Wald test for coefficients: chi2(2) = 194.5, p = < 2.2e-16
#> 

# Blundell-Bond system GMM with collapsed instruments and at most three
# lags of the outcome as instruments
bb <- pangmm(cigs ~ lag(cigs, 1) + price | lag(cigs, 2:99), data = d,
             index = c("id", "wave"), transformation = "ld",
             model = "twostep", collapse = TRUE, maxlag = 4)
summary(bb)
#> 
#> Call:
#> pangmm(formula = cigs ~ lag(cigs, 1) + price | lag(cigs, 2:99), 
#>     data = d, index = c("id", "wave"), transformation = "ld", 
#>     model = "twostep", collapse = TRUE, maxlag = 4)
#> 
#> Estimator: two-way two-step system GMM (Blundell-Bond) 
#> N = 300  T = 7  equations used = 3300  instruments = 12  vcov: robust (Windmeijer) 
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> lag(cigs, 1)  0.50920    0.03825  13.313   <2e-16 ***
#> price        -0.26263    0.02786  -9.428   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Arellano-Bond AR(1) test: z = -8.576, p = < 2.2e-16
#> Arellano-Bond AR(2) test: z = -1.041, p = 0.2978
#> Sargan test (not robust): chi2(4) = 5.905, p = 0.2064
#> Hansen J test (robust): chi2(4) = 4.635, p = 0.3269
#> Wald test for coefficients: chi2(2) = 268.8, p = < 2.2e-16
#> 

# Employment equation of Arellano and Bond (1991)
if (requireNamespace("plm", quietly = TRUE)) {
  data("EmplUK", package = "plm")
  emp <- pangmm(log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 0:1) +
                  log(capital) + lag(log(output), 0:1) | lag(log(emp), 2:99),
                data = EmplUK, index = c("firm", "year"), model = "twostep")
  summary(emp)
}
#> 
#> Call:
#> pangmm(formula = log(emp) ~ lag(log(emp), 1:2) + lag(log(wage), 
#>     0:1) + log(capital) + lag(log(output), 0:1) | lag(log(emp), 
#>     2:99), data = EmplUK, index = c("firm", "year"), model = "twostep")
#> 
#> Estimator: two-way two-step difference GMM (Arellano-Bond) 
#> N = 140  T = 9  equations used = 611  instruments = 38  vcov: robust (Windmeijer) 
#> 
#> Coefficients:
#>                     Estimate Std. Error z value Pr(>|z|)    
#> lag(log(emp), 1)     0.47415    0.18540   2.557 0.010544 *  
#> lag(log(emp), 2)    -0.05297    0.05175  -1.024 0.306051    
#> log(wage)           -0.51320    0.14557  -3.526 0.000423 ***
#> lag(log(wage), 1)    0.22464    0.14195   1.583 0.113528    
#> log(capital)         0.29272    0.06263   4.674 2.95e-06 ***
#> log(output)          0.60977    0.15626   3.902 9.53e-05 ***
#> lag(log(output), 1) -0.44637    0.21730  -2.054 0.039960 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Arellano-Bond AR(1) test: z = -1.538, p = 0.1239
#> Arellano-Bond AR(2) test: z = -0.2797, p = 0.7797
#> Sargan test (not robust): chi2(25) = 73.86, p = 1.013e-06
#> Hansen J test (robust): chi2(25) = 30.11, p = 0.2201
#> Wald test for coefficients: chi2(7) = 142, p = < 2.2e-16
#> 
```
