# Hausman specification test (fixed vs. random effects)

Tests the null that the random-effects estimator is consistent (i.e.
individual effects are uncorrelated with the regressors) against the
alternative that only the fixed-effects estimator is consistent.

## Usage

``` r
panhausman(
  fe,
  re = NULL,
  method = c("classical", "mundlak"),
  data = NULL,
  vcov = c("cluster", "HC1", "classical")
)
```

## Arguments

- fe:

  a `"panglm"` fit with `model = "within"`

- re:

  a `"panglm"` fit with `model = "random"`, same formula/data/family;
  required for `method = "classical"`, ignored for `"mundlak"`

- method:

  `"classical"` (default) or `"mundlak"`

- data:

  for `method = "mundlak"`: the data `fe` was fit on

- vcov:

  for `method = "mundlak"`: covariance of the pooled fit, one of
  `"cluster"` (default), `"HC1"`, `"classical"`

## Value

a list with class `"panglm_hausman"`: `statistic`, `parameter` (df; two
values for an F reference), `p.value`, `method`, `distribution`
(`"chisq"` or `"F"`), and the compared coefficient vectors (`coef_fe`,
`coef_re` for `"classical"`; `coef_mean` for `"mundlak"`)

## Details

`method = "classical"` computes the Hausman (1978) statistic
`t(b_fe - b_re) %*% solve(V_fe - V_re) %*% (b_fe - b_re)` over the
coefficients common to both models (the random-effects intercept is
dropped), with a chi-square reference. It always uses the model-based
(classical) covariances, also when `fe` was fit with a robust `vcov`:
the contrast variance `V_fe - V_re` relies on the random-effects
estimator being efficient under the null, which holds only for
model-based covariances and independent, homoskedastic idiosyncratic
errors.

`method = "mundlak"` is the regression-based (correlated random effects)
version of the same test (Mundlak 1978; Wooldridge 2010, section
10.7.3), valid under heteroskedasticity and serial correlation. The
pooled model `y ~ x + xbar` is fit with the individual means `xbar` of
the time-varying regressors of `fe`, and the null `coef(xbar) = 0` is
tested with a Wald statistic using the `vcov` covariance. With
`vcov = "cluster"` the reference is `F(k, G - 1)` (`G` individuals),
otherwise chi-square(k). Gaussian family only.

## Validation

The classical test is checked against
[`plm::phtest()`](https://rdrr.io/pkg/plm/man/phtest.html) in the
package tests. The Mundlak test is the Wald test on the individual means
of the correlated random effects fit of
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md), whose
pooled and random versions are checked to reproduce the within slope.

## See also

[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md) with
`cre = TRUE`, whose summary reports the Mundlak test.

Other tests:
[`jointtest()`](https://elbadisyimad.com/panglm/reference/jointtest.md),
[`lineartest()`](https://elbadisyimad.com/panglm/reference/lineartest.md),
[`pandispersiontest()`](https://elbadisyimad.com/panglm/reference/pandispersiontest.md)

## Examples

``` r
data(copd)
fit_fe <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
                  model = "within", family = "gaussian")
fit_re <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
                  model = "random", family = "gaussian")
panhausman(fit_fe, fit_re)
#> 
#>  Hausman specification test (fixed vs. random effects)
#> 
#> chisq = 0.08488 , df = 1 , p-value = 0.7708 
#> alternative hypothesis: one model is inconsistent
#> 
panhausman(fit_fe, method = "mundlak", data = copd)
#> 
#>  Mundlak specification test (fixed vs. random effects), vcov: cluster 
#> 
#> F = 0.1647 , df = 1, 59 , p-value = 0.6864 
#> alternative hypothesis: one model is inconsistent
#> 
```
