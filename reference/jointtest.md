# Joint Wald test of coefficients

Tests `H0: R beta = r` with the statistic
`(R b - r)' (R V R')^(-1) (R b - r)`, referred to a chi-squared
distribution with `nrow(R)` degrees of freedom. The coefficients are
chosen by name (`terms`), by a regular expression (`pattern`) or through
an explicit restriction matrix `R`. Typical uses are the joint test of
the spline columns of a regressor (does `x` matter at all), of the
Mundlak means of a correlated random effects fit, or of the leads of an
event study.

## Usage

``` r
jointtest(object, terms = NULL, pattern = NULL, R = NULL, r = 0, vcov = NULL)
```

## Arguments

- object:

  a `"panglm"` fit

- terms:

  character vector of coefficient names

- pattern:

  regular expression selecting coefficient names (used when `terms` and
  `R` are `NULL`)

- R:

  optional restriction matrix with one column per coefficient

- r:

  right-hand side of the restrictions (recycled), default 0

- vcov:

  `NULL` (the covariance stored on the fit) or a type passed to
  [`vcov.panglm()`](https://ielbadisy.github.io/panglm/reference/vcov.panglm.md),
  e.g. `"cluster"`

## Value

an object of class `"panglm_jointtest"`: a list with `statistic`, `df`,
`p.value`, `terms` and `vcov_type`

## Validation

Checked against a hand computation of the Wald statistic, in both the
`terms` and the `R` matrix forms, in the package tests.

## See also

Other tests:
[`lineartest()`](https://ielbadisy.github.io/panglm/reference/lineartest.md),
[`pandispersiontest()`](https://ielbadisy.github.io/panglm/reference/pandispersiontest.md),
[`panhausman()`](https://ielbadisy.github.io/panglm/reference/panhausman.md)

## Examples

``` r
data(copd)
fit <- panglm(fev1 ~ splines::ns(crp, 3), data = copd,
              index = c("id", "visit"), model = "within", vcov = "cluster")
jointtest(fit, pattern = "crp")
#> 
#> Joint Wald test (vcov: cluster)
#> H0: all of splines::ns(crp, 3)1, splines::ns(crp, 3)2, splines::ns(crp, 3)3 are zero
#> Chi-squared = 37.27, df = 3, p-value = 4.031e-08
#> 
```
