# Tidy and glance methods for the panglm model classes

Implement the
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`generics::glance()`](https://generics.r-lib.org/reference/glance.html)
generics for the fits of
[`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md),
[`pangmm()`](https://elbadisyimad.com/panglm/reference/pangmm.md),
[`paniv()`](https://elbadisyimad.com/panglm/reference/paniv.md),
[`panqr()`](https://elbadisyimad.com/panglm/reference/panqr.md),
[`panselect()`](https://elbadisyimad.com/panglm/reference/panselect.md),
[`panfmm()`](https://elbadisyimad.com/panglm/reference/panfmm.md),
[`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md),
[`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md) and
[`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md),
the hurdle model, and for the test objects of the Hausman, joint Wald
and dispersion tests, so that they work with `broom` and `modelsummary`.

## Usage

``` r
# S3 method for class 'panglm_hurdle'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'panglm_hurdle'
glance(x, ...)

# S3 method for class 'twopart'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'twopart'
glance(x, ...)

# S3 method for class 'pangmm'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'pangmm'
glance(x, ...)

# S3 method for class 'paniv'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'paniv'
glance(x, ...)

# S3 method for class 'panqr'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'panqr'
glance(x, ...)

# S3 method for class 'panselect'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'panselect'
glance(x, ...)

# S3 method for class 'panfmm'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'panfmm'
glance(x, ...)

# S3 method for class 'panglm_attgt'
tidy(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  type = c("attgt", "simple", "dynamic", "group", "calendar"),
  ...
)

# S3 method for class 'panglm_attgt'
glance(x, ...)

# S3 method for class 'panglm_etwfe'
tidy(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  type = c("simple", "group", "calendar", "event", "coefficients"),
  ...
)

# S3 method for class 'panglm_etwfe'
glance(x, ...)

# S3 method for class 'panglm_eventstudy'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)

# S3 method for class 'panglm_eventstudy'
glance(x, ...)

# S3 method for class 'panglm_hausman'
tidy(x, ...)

# S3 method for class 'panglm_hausman'
glance(x, ...)

# S3 method for class 'panglm_jointtest'
tidy(x, ...)

# S3 method for class 'panglm_jointtest'
glance(x, ...)

# S3 method for class 'panglm_dispersiontest'
tidy(x, ...)

# S3 method for class 'panglm_dispersiontest'
glance(x, ...)
```

## Arguments

- x:

  a fitted model or test object

- conf.int:

  if `TRUE`, add `conf.low`/`conf.high` columns

- conf.level:

  confidence level for the interval

- ...:

  further arguments passed to the
  [`summary()`](https://rdrr.io/r/base/summary.html) aggregation (attgt,
  etwfe), otherwise unused

- type:

  for [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md)
  and [`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md)
  fits, which effects to tidy (see Details)

## Value

a `data.frame`: one row per coefficient (or effect) for `tidy()` of a
fit, one row for `glance()` and for the test objects (see Details).

## Details

`tidy()` returns one row per coefficient with columns `term`,
`estimate`, `std.error`, `statistic` (the Wald z), `p.value` and, with
`conf.int = TRUE`, the Wald interval `conf.low`, `conf.high`. The
estimates and standard errors are those of
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html), so the interval equals
[`confint()`](https://rdrr.io/r/stats/confint.html). Extra columns
identify the rows where a fit has several blocks of coefficients:
`component` (twopart: `"binary"` or `"positive"`; panfmm: `"class"` or
`"prior"`), `class` (panfmm), `tau` (panqr: one row per term and
quantile), `group`, `time`, `event` (attgt, etwfe, eventstudy).

For [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md),
`type = "attgt"` (the default) gives one row per group-time effect
`ATT(g, t)`, named as in [`coef()`](https://rdrr.io/r/stats/coef.html);
the other types are the aggregations of
[`summary.panglm_attgt()`](https://elbadisyimad.com/panglm/reference/summary.panglm_attgt.md).
For [`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md),
`type = "simple"` (the default) gives the overall ATT, `"group"`,
`"calendar"` and `"event"` the aggregations of
[`summary.panglm_etwfe()`](https://elbadisyimad.com/panglm/reference/summary.panglm_etwfe.md),
and `"coefficients"` the coefficients of the underlying regression. For
[`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md)
the reference period is a row with estimate and standard error zero and
a missing statistic. A statistic is missing whenever its standard error
is zero or missing.

For the test objects, `tidy()` and `glance()` both return a one-row data
frame in the style of `broom`'s `htest` methods: `statistic`, `p.value`,
`parameter` (the degrees of freedom; for an F test also `den.df`),
`distribution` and `method`.

`glance()` returns one row with the estimator, the variance type, the
numbers of observations and individuals, the log-likelihood, AIC and BIC
where the model is fitted by maximum likelihood (twopart, panfmm), and
the specification tests reported by
[`summary()`](https://rdrr.io/r/base/summary.html) (for example the
Sargan, Hansen and Arellano-Bond tests of
[`pangmm()`](https://elbadisyimad.com/panglm/reference/pangmm.md), the
pre-trend test of
[`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md) and
[`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md)).

## Validation

Estimates and standard errors are checked against the reference fits of
each model ([`plm::plm()`](https://rdrr.io/pkg/plm/man/plm.html),
[`plm::pgmm()`](https://rdrr.io/pkg/plm/man/pgmm.html),
[`quantreg::rq()`](https://rdrr.io/pkg/quantreg/man/rq.html),
[`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html)
and
[`etwfe::emfx()`](http://grantmcdermott.com/etwfe/reference/emfx.md)) in
the package's validation suite on GitHub, and the intervals against
[`confint()`](https://rdrr.io/r/stats/confint.html) of the same fit.

## See also

[`tidy.panglm()`](https://elbadisyimad.com/panglm/reference/tidy.panglm.md)
and
[`glance.panglm()`](https://elbadisyimad.com/panglm/reference/glance.panglm.md)
for [`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md)
fits.

## Examples

``` r
data(copd)
if (requireNamespace("generics", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(id = rep(1:100, each = 5), t = rep(1:5, 100), x = rnorm(500))
  d$y <- d$x + rep(rnorm(100), each = 5) + rnorm(500)
  fit <- panqr(y ~ x, data = d, index = c("id", "t"), tau = c(0.25, 0.75))
  generics::tidy(fit, conf.int = TRUE)
  generics::glance(fit)
}
#>   model vcov.type n.tau nobs n.groups
#> 1 canay bootstrap     2  500      100
```
