# Pearson chi-squared/df overdispersion test

A quick diagnostic for count models: the Pearson goodness-of-fit
statistic `sum((y - mu)^2 / mu)` divided by its residual degrees of
freedom should be close to 1 under a correctly-specified Poisson model;
values well above 1 indicate overdispersion (evidence for a negative
binomial specification instead). The diagnostic is available for pooled
and fixed-effects Poisson or negative binomial models. It is not applied
to random-effects models because their marginal variance includes the
fitted mixing distribution and is not equal to the observation-level
Poisson variance used by this statistic.

## Usage

``` r
pandispersiontest(object)
```

## Arguments

- object:

  a `"panglm"` fit with `family` `"poisson"` or `"negbin"`

## Value

a list of class `"panglm_dispersiontest"` with the Pearson statistic,
its degrees of freedom, the ratio (statistic/df), a one-sided p-value
against `chisq(df)`, and `n_excluded`: the number of observations
dropped from the statistic because their fitted value was exactly zero
or unavailable (e.g. rows in a structurally-zero panel individual
screened out of a fixed-effects negbin fit)

## Validation

No reference implementation exists for panel fits; the package tests
check that it flags simulated overdispersed panels and not equidispersed
ones.

## See also

[`negbinfamily()`](https://ielbadisy.github.io/panglm/reference/panglm-families.md)
for the negative binomial alternative.

Other tests:
[`jointtest()`](https://ielbadisy.github.io/panglm/reference/jointtest.md),
[`lineartest()`](https://ielbadisy.github.io/panglm/reference/lineartest.md),
[`panhausman()`](https://ielbadisy.github.io/panglm/reference/panhausman.md)

## Examples

``` r
data(copd)
fit <- panglm(exacerbations ~ treatment + age + smoker + crp, data = copd,
              index = c("id", "visit"), model = "pooling", family = "poisson")
pandispersiontest(fit)
#> 
#>  Pearson chi-squared/df overdispersion test
#> 
#> chisq = 378.5 , df = 295 , ratio = 1.283 , p-value = 0.0007215 
#> alternative hypothesis: overdispersion (ratio > 1)
#> 
```
