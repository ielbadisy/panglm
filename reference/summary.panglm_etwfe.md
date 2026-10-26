# Average treatment effects on the treated from an ETWFE fit

Averages, over the treated observations, the difference between the
predicted mean with and without the observation's cohort by period
effect (on the response scale). The averages are taken over all treated
observations (`"simple"`), by cohort (`"group"`), by calendar period
(`"calendar"`) or by time since treatment (`"event"`), with delta-method
standard errors. With sampling weights in
[`etwfe()`](https://ielbadisy.github.io/panglm/reference/etwfe.md), the
averages are weighted. They correspond to
[`etwfe::emfx()`](http://grantmcdermott.com/etwfe/reference/emfx.md).
With `control = "never"`, the event study also reports the pre-treatment
placebo effects, with event time `-1` as the reference (zero).

## Usage

``` r
# S3 method for class 'panglm_etwfe'
summary(
  object,
  type = c("simple", "group", "calendar", "event"),
  level = object$level,
  byxvar = !is.null(object$xvar),
  ...
)
```

## Arguments

- object:

  a `"panglm_etwfe"` fit from
  [`etwfe()`](https://ielbadisy.github.io/panglm/reference/etwfe.md)

- type:

  `"simple"`, `"group"`, `"calendar"` or `"event"`

- level:

  confidence level; defaults to that of the fit

- byxvar:

  logical; for a fit with `xvar`, report the ATTs separately for each
  value of the covariate (by aggregation level and value). The default
  is `TRUE` for such fits; `overall` is always the ATT over all treated
  observations.

- ...:

  unused

## Value

an object of class `"summary.panglm_etwfe"`: a list with `overall` (the
simple ATT: `estimate`, `se`, `statistic`, `p.value`, `lower`, `upper`),
`estimates` (the same columns by `level` of the aggregation, and by
value of the `xvar` covariate in a column named after it when
`byxvar = TRUE`; `NULL` for `"simple"` without `byxvar`), `cells` and
the fit information.

## See also

Other etwfe methods:
[`plot.panglm_etwfe()`](https://ielbadisy.github.io/panglm/reference/plot.panglm_etwfe.md)

## Examples

``` r
set.seed(2)
d <- expand.grid(id = 1:150, t = 1:8)
d$first <- rep(c(4, 6, NA), each = 50)[d$id]
d$y <- rbinom(nrow(d), 1, plogis(-1 + 0.1 * d$t +
  0.8 * (!is.na(d$first) & d$t >= d$first)))
fit <- etwfe(y ~ 1, d, index = c("id", "t"), treat = "first", family = "binomial")
summary(fit, type = "event")
#> 
#> Extended TWFE (Wooldridge), logit model, control group: not yet treated
#> Observations: 1200  Individuals: 150  vcov: cluster 
#> Fixed effects: cohort and period
#> 
#> ATT (response scale): 0.2108 (SE 0.06119, 95% CI 0.09088 to 0.3308, p = 0.0005711)
#> 
#>  event estimate      se statistic  p.value    lower  upper
#>      0   0.2275 0.07474     3.044 0.002335  0.08102 0.3740
#>      1   0.1772 0.07631     2.322 0.020212  0.02766 0.3268
#>      2   0.2018 0.08722     2.314 0.020657  0.03090 0.3728
#>      3   0.1862 0.12379     1.504 0.132495 -0.05640 0.4288
#>      4   0.2871 0.10781     2.663 0.007735  0.07584 0.4985
#> 
```
