# Tidy a panglm fit into a data frame of coefficients

Implements the
[`generics::tidy()`](https://generics.r-lib.org/reference/tidy.html)
generic so `panglm` results work with `broom`/`modelsummary`
table-building workflows.

## Usage

``` r
# S3 method for class 'panglm'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)
```

## Arguments

- x:

  a `"panglm"` fit

- conf.int:

  if `TRUE`, add `conf.low`/`conf.high` columns

- conf.level:

  confidence level for the interval

- ...:

  unused

## Value

a `data.frame` with one row per coefficient

## See also

[tidy-methods](https://ielbadisy.github.io/panglm/reference/tidy-methods.md)
for the other model classes.

Other panglm methods:
[`confint.panglm()`](https://ielbadisy.github.io/panglm/reference/confint.panglm.md),
[`effectcurve()`](https://ielbadisy.github.io/panglm/reference/effectcurve.md),
[`fitted.panglm()`](https://ielbadisy.github.io/panglm/reference/fitted.panglm.md),
[`glance.panglm()`](https://ielbadisy.github.io/panglm/reference/glance.panglm.md),
[`logLik.panglm()`](https://ielbadisy.github.io/panglm/reference/logLik.panglm.md),
[`model.matrix.panglm()`](https://ielbadisy.github.io/panglm/reference/model.matrix.panglm.md),
[`partialeffects.panordered()`](https://ielbadisy.github.io/panglm/reference/partialeffects.md),
[`plot.panglm()`](https://ielbadisy.github.io/panglm/reference/plot.panglm.md),
[`plot.panglm_effectcurve()`](https://ielbadisy.github.io/panglm/reference/plot.panglm_effectcurve.md),
[`predict.panglm()`](https://ielbadisy.github.io/panglm/reference/predict.panglm.md),
[`predict.panordered()`](https://ielbadisy.github.io/panglm/reference/predict.panordered.md),
[`summary.panglm()`](https://ielbadisy.github.io/panglm/reference/summary.panglm.md),
[`vcov.panglm()`](https://ielbadisy.github.io/panglm/reference/vcov.panglm.md)

## Examples

``` r
data(copd)
fit <- panglm(fev1 ~ treatment + age + smoker + crp, data = copd,
              index = c("id", "visit"), model = "pooling", family = "gaussian")
if (requireNamespace("generics", quietly = TRUE)) {
  generics::tidy(fit, conf.int = TRUE)
}
#>          term    estimate   std.error statistic      p.value    conf.low
#> 1 (Intercept)  4.09147857 0.225597468 18.136190 1.650863e-73  3.64931566
#> 2   treatment  0.20855645 0.056577512  3.686207 2.276208e-04  0.09766656
#> 3         age -0.01423926 0.003133068 -4.544830 5.497949e-06 -0.02037996
#> 4      smoker -0.35520702 0.055904799 -6.353784 2.100818e-10 -0.46477841
#> 5         crp -0.03289479 0.013135302 -2.504304 1.226925e-02 -0.05863951
#>      conf.high
#> 1  4.533641481
#> 2  0.319446334
#> 3 -0.008098562
#> 4 -0.245635625
#> 5 -0.007150074
```
