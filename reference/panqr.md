# Quantile regression for panel data with individual effects

Fixed-effects quantile regression for linear panel models: the two-step
estimator of Canay (2011) and the method of moments quantile regression
(MM-QR) of Machado and Santos Silva (2019), the estimator of the Stata
command `xtqreg`. A pooled quantile regression that ignores the
individual effects is available as a baseline. Several quantiles can be
estimated in one call.

## Usage

``` r
panqr(
  formula,
  data,
  index,
  tau = 0.5,
  model = c("canay", "mmqr", "pooling"),
  vcov = NULL,
  boot = 200,
  na.action = stats::na.omit,
  seed = NULL,
  R,
  cluster = NULL
)
```

## Arguments

- formula:

  a formula `y ~ regressors`. `lag(x, k)`, `lead(x, k)` and `diff(x)`
  are taken within each individual on the grid of sorted distinct time
  values (see "Lags and differences" in
  [`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md);
  needs `index = c(id, time)`); the first waves of each individual get
  `NA` and leave the sample, and `lag(x, 0:2)` gives several lags

- data:

  a data.frame in long format

- index:

  length-1 or length-2 character vector giving the individual (and time)
  columns

- tau:

  numeric vector of quantiles, each strictly between 0 and 1

- model:

  `"canay"`, `"mmqr"` or `"pooling"`

- vcov:

  `"analytic"` (cluster-robust), `"iid"` (MM-QR only), `"bootstrap"` or
  `"none"`; the default is `"bootstrap"` for `"canay"` and `"analytic"`
  otherwise

- boot:

  number of bootstrap replications for `vcov = "bootstrap"`

- na.action:

  missing-data action passed to
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)

- seed:

  optional random seed for the bootstrap; the global random number state
  is restored afterwards. With `seed = NULL` (the default) the draws
  come from the global random number stream.

- R:

  deprecated; use `boot`

- cluster:

  optional name of a column of `data` to cluster the analytic covariance
  on (`vcov = "analytic"` only, including the MM-QR location and scale
  covariances); the default clusters by individual

## Value

an object of class `"panqr"`, a list with components `coefficients` (a
matrix with one column per quantile), `vcov` (a list of covariance
matrices, one per quantile, or `NULL` for `vcov = "none"`), `vcov_type`,
`tau`, `model`, `boot` (the bootstrap replicates, an array
`boot x k x length(tau)`, for `vcov = "bootstrap"`), `alpha` (one per
individual: for `"canay"` the normalized effects subtracted from the
outcome, for `"mmqr"` the location effects including the constant), for
`"mmqr"` also `location`, `scale` (lists with `coefficients` and
`vcov`), `q` (the quantiles of the standardized errors) and `delta` (the
individual scale effects including the constant, so that
`alpha_i(tau) = alpha_i + delta_i q(tau)`), `cluster` and `n_clusters`
(for the analytic covariances), `dropped` (time-invariant regressors
dropped), `n_singletons` (individuals observed once, dropped for
`"canay"` and `"mmqr"`), `nobs`, `n_groups`, `call` and `formula`.
Methods: `print`, `summary`, `coef`, `vcov`, `confint`, `nobs`, `plot`,
`tidy` and `glance`.

## Estimators

- `"canay"`: the individual effects are estimated from the within (mean)
  regression, `alpha_i = ybar_i - xbar_i'b - c` with `c = ybar - xbar'b`
  the overall constant (the normalization of Stata's `xtreg, fe`), and
  the quantile regression of the adjusted outcome `y - alpha_i` on the
  regressors (with an intercept) is fitted by
  [`quantreg::rq.fit()`](https://rdrr.io/pkg/quantreg/man/rq.fit.html)
  with `method = "br"`. The intercept is the `tau`-quantile of
  `c + error`. The model assumes that the individual effects are pure
  location shifts, and the estimator is consistent when both `N` and `T`
  grow.

- `"mmqr"`: the location-scale model
  `y_it = alpha_i + x_it'b + (delta_i + x_it'g) U_it`. The location
  parameters `b` come from the within regression of `y` on `x`; with `R`
  its residuals, the scale parameters `g` come from the within
  regression of `2 R (I(R >= 0) - mean(I(R >= 0)))`, an estimate of
  `|R|` that is robust to the first step; `q(tau)` is the sample
  `tau`-quantile of the standardized residuals `R / sigma_it`, with
  `sigma_it` the fitted scale including the individual effects. The
  quantile coefficients are `b(tau) = b + g q(tau)`. The scale variables
  are the regressors, as in `xtqreg`. Only slopes are reported: the
  intercept is absorbed in the individual effects
  `alpha_i(tau) = alpha_i + delta_i q(tau)`.

- `"pooling"`:
  [`quantreg::rq.fit()`](https://rdrr.io/pkg/quantreg/man/rq.fit.html)
  of `y` on the regressors, ignoring the panel structure (except for the
  clustered covariance).

Time-invariant regressors are not identified with individual effects and
are dropped (with a warning) for `"canay"` and `"mmqr"`, as are
individuals observed only once (with a message).

## Covariance

`vcov = "bootstrap"` resamples whole individuals with replacement (a
cluster bootstrap) `boot` times, repeats every step of the estimator and
returns the covariance of the replicates. `vcov = "analytic"` gives:

- `"mmqr"`: a cluster-robust (by individual, or by `cluster`) sandwich
  built from the influence function of `b(tau) = b + g q(tau)`. The
  influence of observation `it` adds the three steps of the estimator:
  the location within regression, `Q^{-1} x~_it R_it`; the scale within
  regression, `q(tau) Q^{-1} x~_it sigma_it V_it`; and the quantile
  `q(tau)`, `g s_i w_it / sum(s)`, with `Q = sum x~ x~'`, `s_i` the
  scale of individual `i` at the regressor means and
  `w_it = (tau - I(U_it <= q(tau))) / f(q(tau)) - U_it - q(tau) V_it`.
  These are summed within clusters into `psi_g`, and
  `V = G / (G - 1) * sum_g psi_g psi_g'` for `G` clusters. The density
  of the standardized error at `q(tau)` is the difference quotient
  `2h / (q(tau + h) - q(tau - h))` (as Stata's `qreg, vce(iid)`) with
  the Hall and Sheather bandwidth `h`. The scale influence uses the
  centered `V = 2 U (I(U >= 0) - mean(I(U >= 0))) - 1` of the paper
  (also used by the Stata command `mmqreg`). Unlike the iid formula, it
  allows the errors of an individual to be dependent (serial
  correlation, or an individual component in `U`) and the errors to be
  heteroskedastic.

`vcov = "iid"` (MM-QR only) gives the asymptotic covariance of Machado
and Santos Silva (2019) under independent observations, with the matrix
layout of `xtqreg` (version 1.5) and the same density estimate and
centered `V`; `xtqreg` 1.5 omits the `- 1` in `V`, which inflates its
standard errors away from the median. This was the analytic covariance
of `"mmqr"` before panglm 1.6.1. It replaces observation-level products
by their sample means and ignores the dependence within individuals; it
can be far off with short panels (see the package's validation suite).

The other models with `vcov = "analytic"`:

- `"canay"`: a cluster-robust (by individual) sandwich built from the
  influence function of the two-step estimator, which corrects the
  quantile regression score `(tau - I(e < 0)) x` for the estimation of
  the individual effects and of the within slopes. The conditional
  density is estimated with the Gaussian kernel and bandwidth of
  `quantreg::summary.rq(se = "ker")`. This is a large-`N`, large-`T`
  approximation; Canay (2011) recommends the bootstrap. It treats the
  individual effects as fixed, so its intercept SE omits the sampling
  variation of the mean individual effect that the cluster bootstrap
  includes; the slope SEs of the two agree.

- `"pooling"`: the same kernel sandwich without the first-step terms,
  clustered by individual.

The `"canay"` and `"pooling"` cluster sandwiches carry no finite-sample
factor; the `"mmqr"` one carries `G / (G - 1)`. For `"mmqr"` the
location and scale coefficients are returned with the cluster-robust
within covariance of the rest of the package
(`G / (G - 1) * (n - 1) / (n - k)`, as `plm::vcovHC(type = "sss")`).

## Validation

Checked in the package's validation suite on GitHub: the pooled model
against [`quantreg::rq()`](https://rdrr.io/pkg/quantreg/man/rq.html) and
its kernel sandwich; the Canay and MM-QR estimators against
implementations built from
[`plm::plm()`](https://rdrr.io/pkg/plm/man/plm.html),
[`stats::lm()`](https://rdrr.io/r/stats/lm.html) and
[`quantreg::rq()`](https://rdrr.io/pkg/quantreg/man/rq.html), with the
MM-QR iid covariance checked against the formulas of the Stata command
`xtqreg`; the MM-QR cluster covariance against the cluster bootstrap and
by Monte Carlo coverage with serially correlated errors; and Monte Carlo
checks that the analytic standard errors match the sampling standard
deviation.

## References

Canay, I. A. (2011). A simple approach to quantile regression for panel
data. *The Econometrics Journal*, 14(3), 368-386.
[doi:10.1111/j.1368-423X.2011.00349.x](https://doi.org/10.1111/j.1368-423X.2011.00349.x)

Machado, J. A. F. and Santos Silva, J. M. C. (2019). Quantiles via
moments. *Journal of Econometrics*, 213(1), 145-173.
[doi:10.1016/j.jeconom.2019.04.009](https://doi.org/10.1016/j.jeconom.2019.04.009)

## See also

[`summary.panqr()`](https://elbadisyimad.com/panglm/reference/summary.panqr.md),
[`vcov.panqr()`](https://elbadisyimad.com/panglm/reference/vcov.panqr.md)
and
[`plot.panqr()`](https://elbadisyimad.com/panglm/reference/plot.panqr.md)
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
[`panselect()`](https://elbadisyimad.com/panglm/reference/panselect.md),
[`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md)

## Examples

``` r
# Simulated panel of 300 patients over 6 years. Annual health costs have
# a patient effect and a heavy upper tail that widens with chronic
# conditions: a location-scale model in which the effect of chronic
# conditions grows across the cost distribution.
set.seed(1)
N <- 300; TT <- 6
d <- data.frame(id = rep(seq_len(N), each = TT), year = rep(seq_len(TT), N))
alpha <- rep(rnorm(N), each = TT)
d$chronic <- rpois(N * TT, 1 + 0.5 * (alpha > 0))
d$age <- rep(sample(40:70, N, TRUE), each = TT) + d$year
d$cost <- 2 + alpha + 0.5 * d$chronic + 0.02 * d$age +
  (1 + 0.3 * d$chronic) * rexp(N * TT)

if (requireNamespace("quantreg", quietly = TRUE)) {
  # median and upper tail of costs, MM-QR with analytic standard errors
  mm <- panqr(cost ~ chronic + age, data = d, index = c("id", "year"),
              tau = c(0.5, 0.9), model = "mmqr")
  summary(mm)

  # Canay two-step estimator with a small cluster bootstrap
  cn <- panqr(cost ~ chronic + age, data = d, index = c("id", "year"),
              tau = 0.9, model = "canay", boot = 50, seed = 1)
  coef(cn)
}
#> Warning: 0.1% of the fitted values of the scale function are not positive
#> Warning: Solution may be nonunique
#> Warning: Solution may be nonunique
#> (Intercept)     chronic         age 
#>  4.11072454  1.19437845  0.02214631 
```
