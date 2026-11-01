# Fixed-effects hurdle model for panel count data

A hurdle decomposition for panels where a large fraction of groups
report *structural* zeros (never observed to occur), rather than zeros
arising from ordinary Poisson/NB sampling variation at a low mean, a
pattern that shows up as extreme, family-invariant overdispersion under
either `family = "poisson"` or `family = "negbin"` in
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md).

## Usage

``` r
panhurdle(
  formula,
  data,
  index,
  family = c("poisson", "negbin"),
  maxit = 100,
  tol = 1e-10
)
```

## Arguments

- formula:

  a model formula, e.g. `y ~ x1 + x2`

- data:

  a data.frame containing the variables in `formula` and `index`

- index:

  length-2 character vector identifying the individual and time columns,
  e.g. `c("id", "time")`

- family:

  count distribution for positive observations, either `"poisson"` or
  `"negbin"` (NB2 parameterization). It is stored on the fit as
  `count_family`.

- maxit:

  maximum IRLS/Newton iterations

- tol:

  convergence tolerance

## Value

an object of class `"panglm_hurdle"`: a list with `zero` (the
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md)
binomial fit for the zero-vs-positive part) and `count` (the
zero-truncated fixed-effects fit for the count part). Methods: `print`,
`summary`, `coef`, `vcov`, `confint`, `nobs`, `plot`, `tidy` and
`glance`.

## Details

The model is fit in two independent parts, both using fixed effects, and
reported together rather than folded into
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md)'s
single-family dispatch, since they operate on different sample sizes
(the count part drops the zeros):

1.  **Zero vs. positive.** `1(y > 0)` fit via the exact conditional
    logistic regression already used for
    `panglm(..., family = "binomial")` (Chamberlain 1980): no new
    estimator, just a new outcome vector.

2.  **Count given `y > 0`.** A zero-truncated Poisson or NB2
    fixed-effects model (Allison-Waterman-style dummy-variable
    intercepts, the same unconditional approach used for FE-NB2) on the
    subsample where `y > 0`. The log-density is the ordinary count
    log-density minus `log(1 - P(Y = 0))`; see
    [`pscl::hurdle()`](https://rdrr.io/pkg/pscl/man/hurdle.html) for the
    reference form of the truncated-count likelihood.

Groups with no remaining within-group variation after truncation (a
single surviving observation) contribute no information about the shared
count-part slope, the same profiling fact that applies to the other
fixed-effects estimators in this package; they are not dropped
explicitly since a perfectly-fit singleton doesn't diverge, it simply
carries zero weight. Groups whose positive counts are all equal to 1 are
dropped from the count part with a message: their intercept has no
finite estimate (the truncated likelihood approaches its maximum as the
group mean goes to 0, whatever the slopes), so they carry no information
either.

Both parts remove one effect per individual, so regressors that do not
vary within individuals (sex, baseline covariates) are not identified:
they are dropped with a warning, as in `panglm(model = "within")`. A
regressor that varies only in individuals a part does not use (for
example individuals with a single positive count, for the count part) is
dropped from that part alone, with a warning.

The count part is fit by Newton-Raphson with the analytic gradient and
Hessian of the joint likelihood. The group intercepts enter through a
diagonal Hessian block that is eliminated at each step, so the cost is
linear in the number of groups and the intercepts are never stored as
dummy columns. Standard errors are the slope block of the inverse
observed information of the joint likelihood (slopes, intercepts and,
for NB2, `log(theta)`). When the NB2 shape parameter runs off to
infinity (no overdispersion), it is stopped at `1e8` with a warning and
the fit equals the zero-truncated Poisson one.

## Validation

Numerical accuracy is checked in the package tests and in the validation
suite; see the validation article on the package website.

## See also

[hurdle-methods](https://elbadisyimad.com/panglm/reference/hurdle-methods.md)
and
[`plot.panglm_hurdle()`](https://elbadisyimad.com/panglm/reference/plot.panglm_hurdle.md)
for the fitted model,
[`pandispersiontest()`](https://elbadisyimad.com/panglm/reference/pandispersiontest.md)
to check the count part,
[`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md) for
continuous costs.

Other panel models:
[`panfmm()`](https://elbadisyimad.com/panglm/reference/panfmm.md),
[`panfmmselect()`](https://elbadisyimad.com/panglm/reference/panfmmselect.md),
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md),
[`pangmm()`](https://elbadisyimad.com/panglm/reference/pangmm.md),
[`paniv()`](https://elbadisyimad.com/panglm/reference/paniv.md),
[`panqr()`](https://elbadisyimad.com/panglm/reference/panqr.md),
[`panselect()`](https://elbadisyimad.com/panglm/reference/panselect.md),
[`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md)

## Examples

``` r
data(copd)
fit <- panhurdle(exacerbations ~ crp, data = copd, index = c("id", "visit"))
#> 26 group(s) with no within-group outcome variation dropped (all-0 or all-1) -- they carry no information for the conditional likelihood.
#> 6 group(s) whose positive counts are all 1 dropped from the count part -- their intercept is not finite and they carry no information.
fit
#> Fixed-effects hurdle model (zero part: conditional logit; count part: zero-truncated Poisson)
#> 
#> Zero part, 1(y > 0):
#>    crp 
#> 0.1089 
#> 
#> Count part, y | y > 0:
#>     crp 
#> 0.04246 
#> 
#> n = 300  positive = 238 
```
