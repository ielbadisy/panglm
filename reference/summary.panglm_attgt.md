# Aggregate group-time average treatment effects

Combines the `ATT(g, t)` of an
[`attgt()`](https://ielbadisy.github.io/panglm/reference/attgt.md) fit
into summary effects, as in Callaway and Sant'Anna (2021, Section 3) and
[`did::aggte()`](https://bcallaway11.github.io/did/reference/aggte.html):

- `"simple"`: the average of all post-treatment `ATT(g, t)` weighted by
  cohort size;

- `"dynamic"`: the event study, `ATT(e)` for each event time
  `e = t - g`, averaging across cohorts weighted by cohort size; the
  overall effect is the mean of the post-treatment `ATT(e)`, `e >= 0`;

- `"group"`: the average effect of each cohort over its post-treatment
  periods; the overall effect weights cohorts by size;

- `"calendar"`: the average effect in each calendar period across the
  cohorts treated by then (weighted by size); the overall effect is the
  mean over periods.

## Usage

``` r
# S3 method for class 'panglm_attgt'
summary(
  object,
  type = c("simple", "dynamic", "group", "calendar"),
  min_e = -Inf,
  max_e = Inf,
  balance_e = NULL,
  level = object$level,
  ...
)
```

## Arguments

- object:

  a `"panglm_attgt"` fit from
  [`attgt()`](https://ielbadisy.github.io/panglm/reference/attgt.md)

- type:

  `"simple"`, `"dynamic"`, `"group"` or `"calendar"`

- min_e, max_e:

  for `"dynamic"`: event-time window

- balance_e:

  for `"dynamic"`: if not `NULL`, only cohorts observed for at least
  `balance_e` periods after treatment are used, so that the composition
  of cohorts is the same at every event time up to `balance_e`

- level:

  confidence level; defaults to that of the fit

- ...:

  unused

## Value

an object of class `"summary.panglm_attgt"`: a list with `overall` (a
one-row data frame with `estimate`, `se`, `lower`, `upper`), `estimates`
(the effects by event time, cohort or period: column `level` plus the
same columns; `NULL` for `"simple"`), `type`, `pretrend`, `inf` (the
influence functions of the overall effect), `crit` (the simultaneous
critical value of the uniform bands, or `NULL`), `boot` and the fit
settings. The influence functions of the effects by level are the
attribute `"inf"` of `estimates`.

## Details

Standard errors come from the influence functions and account for the
estimation of the cohort-size weights. When the fit was bootstrapped
([`attgt()`](https://ielbadisy.github.io/panglm/reference/attgt.md) with
`bstrap = TRUE`), the standard errors of the aggregated effects come
from the same multiplier draws and, with `cband = TRUE`, the effects by
level get uniform confidence bands (`band_lower`, `band_upper`) from the
simultaneous critical value `crit`, as in
[`did::aggte()`](https://bcallaway11.github.io/did/reference/aggte.html).
Clustering follows the fit.

## See also

Other attgt methods:
[`plot.panglm_attgt()`](https://ielbadisy.github.io/panglm/reference/plot.panglm_attgt.md)

## Examples

``` r
set.seed(2)
d <- expand.grid(id = 1:150, t = 1:8)
d$first <- rep(c(4, 6, NA), each = 50)[d$id]
d$y <- rnorm(150)[d$id] + 0.2 * d$t +
  ifelse(!is.na(d$first) & d$t >= d$first, 1, 0) + rnorm(nrow(d))
fit <- attgt(y ~ 1, d, index = c("id", "t"), treat = "first")
summary(fit, type = "dynamic")
#> 
#> Callaway and Sant'Anna group-time ATT, aggregated event study (by time since treatment)
#> Control group: never treated  Estimator: unconditional  Individuals: 150 
#> 
#> ATT (mean of e >= 0): 1.176 (SE 0.163, 95% CI 0.8569 to 1.496)
#> 
#>  event estimate     se    lower    upper
#>     -4   0.4697 0.2704 -0.06039  0.99972
#>     -3  -0.5765 0.2717 -1.10908 -0.04401
#>     -2   0.1442 0.1984 -0.24472  0.53316
#>     -1  -0.2597 0.2103 -0.67194  0.15259
#>      0   0.8400 0.2041  0.44001  1.24005
#>      1   1.1887 0.1602  0.87475  1.50275
#>      2   1.3308 0.1882  0.96196  1.69970
#>      3   1.2700 0.2685  0.74374  1.79633
#>      4   1.2525 0.2656  0.73181  1.77311
#> 
#> Pre-trend test (pre-treatment ATT(g, t) = 0): Chi-squared = 22.19, df = 6, p-value = 0.001119
#> 
summary(fit, type = "calendar")
#> 
#> Callaway and Sant'Anna group-time ATT, aggregated by calendar period
#> Control group: never treated  Estimator: unconditional  Individuals: 150 
#> 
#> ATT (mean over periods): 1.141 (SE 0.1613, 95% CI 0.825 to 1.457)
#> 
#>  time estimate     se  lower upper
#>     4   0.8908 0.2821 0.3378 1.444
#>     5   1.2792 0.2631 0.7635 1.795
#>     6   1.0762 0.2131 0.6585 1.494
#>     7   1.1842 0.1931 0.8056 1.563
#>     8   1.2755 0.2099 0.8641 1.687
#> 
#> Pre-trend test (pre-treatment ATT(g, t) = 0): Chi-squared = 22.19, df = 6, p-value = 0.001119
#> 
```
