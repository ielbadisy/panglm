# Test whether the effect of a regressor is linear

For a regressor entered through a spline or polynomial basis whose span
contains the linear function (true for
[`splines::ns()`](https://rdrr.io/r/splines/ns.html),
[`splines::bs()`](https://rdrr.io/r/splines/bs.html) and
[`poly()`](https://rdrr.io/r/stats/poly.html)), linearity means the
basis coefficients are proportional to the coefficients `a` that
reproduce `x` itself (`x = c + B a`). The test is the joint Wald test of
the `q - 1` restrictions that the coefficient vector is orthogonal to
the complement of `a`, so it accepts any covariance type, including
cluster-robust ones. For a test of whether the regressor matters at all,
use
[`jointtest()`](https://ielbadisy.github.io/panglm/reference/jointtest.md)
on its columns.

## Usage

``` r
lineartest(object, variable, vcov = NULL)
```

## Arguments

- object:

  a `"panglm"` fit

- variable:

  name of the regressor

- vcov:

  `NULL` (the covariance stored on the fit) or a type passed to
  [`vcov.panglm()`](https://ielbadisy.github.io/panglm/reference/vcov.panglm.md)

## Value

an object of class `"panglm_jointtest"`, as returned by
[`jointtest()`](https://ielbadisy.github.io/panglm/reference/jointtest.md).

## Validation

The package tests check that the test rejects a nonlinear curve and
keeps its size under linearity, and the fixed-effects likelihood-ratio
statistic against the Wald statistic.

## See also

Other tests:
[`jointtest()`](https://ielbadisy.github.io/panglm/reference/jointtest.md),
[`pandispersiontest()`](https://ielbadisy.github.io/panglm/reference/pandispersiontest.md),
[`panhausman()`](https://ielbadisy.github.io/panglm/reference/panhausman.md)

## Examples

``` r
data(copd)
fit <- panglm(fev1 ~ splines::ns(crp, 3), data = copd,
              index = c("id", "visit"), model = "within", vcov = "cluster")
lineartest(fit, "crp")
#> 
#> Joint Wald test (vcov: cluster)
#> H0: the effect of crp is linear 
#> Chi-squared = 1.338, df = 2, p-value = 0.5122
#> 
```
