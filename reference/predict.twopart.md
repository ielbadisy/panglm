# Predict from a two-part model

Predict from a two-part model

## Usage

``` r
# S3 method for class 'twopart'
predict(object, newdata = NULL, type = c("response", "prob", "positive"), ...)
```

## Arguments

- object:

  a `"twopart"` fit

- newdata:

  optional data frame; defaults to the (complete-case) data used for
  fitting. For `cre = TRUE` fits, the individual means are recomputed
  within `newdata` unless its `mean_*` columns are supplied.

- type:

  `"response"` for `E[y] = P(y > 0) E[y | y > 0]`, `"prob"` for
  `P(y > 0)`, or `"positive"` for `E[y | y > 0]` (retransformed for the
  lognormal part)

- ...:

  unused

## Value

a numeric vector with one prediction per row of `newdata`

## See also

Other twopart methods:
[`confint.twopart()`](https://elbadisyimad.com/panglm/reference/confint.twopart.md)

## Examples

``` r
data(copd)
set.seed(1)
copd$cost <- ifelse(copd$exacerbations > 0,
                    exp(6 + 0.1 * copd$crp + rnorm(nrow(copd), 0, 0.5)), 0)
tp <- twopart(cost ~ treatment + crp, data = copd, index = c("id", "visit"))
head(predict(tp))                    # E[y]
#> [1] 657.3327 376.3459 552.9015 608.6868 511.9711 558.2423
head(predict(tp, type = "prob"))     # P(y > 0)
#> [1] 0.7720320 0.7609261 0.7686167 0.7705176 0.7670897 0.7688072
head(predict(tp, type = "positive")) # E[y | y > 0]
#> [1] 851.4320 494.5893 719.3462 789.9713 667.4202 726.1148
```
