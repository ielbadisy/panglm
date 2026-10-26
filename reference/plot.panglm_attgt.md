# Plot an event study from group-time effects

Plots the dynamic aggregation of
[`attgt()`](https://ielbadisy.github.io/panglm/reference/attgt.md)
(`ATT(e)` by event time with pointwise confidence intervals), or the
cohort or calendar aggregation. When the fit has uniform confidence
bands (`cband = TRUE`), they are drawn as thin grey whiskers behind the
pointwise intervals.

## Usage

``` r
# S3 method for class 'panglm_attgt'
plot(
  x,
  type = c("dynamic", "group", "calendar"),
  xlab = NULL,
  ylab = "ATT",
  ...
)

# S3 method for class 'summary.panglm_attgt'
plot(x, type = NULL, xlab = NULL, ylab = "ATT", ...)
```

## Arguments

- x:

  a `"panglm_attgt"` fit from
  [`attgt()`](https://ielbadisy.github.io/panglm/reference/attgt.md), or
  its summary

- type:

  aggregation to plot (for a fit): `"dynamic"`, `"group"` or
  `"calendar"`

- xlab, ylab:

  axis labels

- ...:

  passed to
  [`summary.panglm_attgt()`](https://ielbadisy.github.io/panglm/reference/summary.panglm_attgt.md)
  (for example `min_e`, `max_e`) and then to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html)

## Value

the plotted summary, invisibly

## See also

Other attgt methods:
[`summary.panglm_attgt()`](https://ielbadisy.github.io/panglm/reference/summary.panglm_attgt.md)

## Examples

``` r
set.seed(2)
d <- expand.grid(id = 1:150, t = 1:8)
d$first <- rep(c(4, 6, NA), each = 50)[d$id]
d$y <- rnorm(150)[d$id] + 0.2 * d$t +
  ifelse(!is.na(d$first) & d$t >= d$first, 1, 0) + rnorm(nrow(d))
plot(attgt(y ~ 1, d, index = c("id", "t"), treat = "first"))
```
