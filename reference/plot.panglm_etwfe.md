# Plot ETWFE treatment effects

Plots the event-study aggregation of an
[`etwfe()`](https://ielbadisy.github.io/panglm/reference/etwfe.md) fit
(or the cohort or calendar aggregation) with pointwise confidence
intervals.

## Usage

``` r
# S3 method for class 'panglm_etwfe'
plot(x, type = c("event", "group", "calendar"), xlab = NULL, ylab = "ATT", ...)

# S3 method for class 'summary.panglm_etwfe'
plot(x, type = NULL, xlab = NULL, ylab = "ATT", ...)
```

## Arguments

- x:

  a `"panglm_etwfe"` fit from
  [`etwfe()`](https://ielbadisy.github.io/panglm/reference/etwfe.md), or
  its summary

- type:

  aggregation to plot (for a fit): `"event"`, `"group"` or `"calendar"`

- xlab, ylab:

  axis labels

- ...:

  passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html)

## Value

the plotted summary, invisibly

## See also

Other etwfe methods:
[`summary.panglm_etwfe()`](https://ielbadisy.github.io/panglm/reference/summary.panglm_etwfe.md)

## Examples

``` r
set.seed(2)
d <- expand.grid(id = 1:150, t = 1:8)
d$first <- rep(c(4, 6, NA), each = 50)[d$id]
d$y <- rnorm(150)[d$id] + 0.2 * d$t +
  ifelse(!is.na(d$first) & d$t >= d$first, 1, 0) + rnorm(nrow(d))
plot(etwfe(y ~ 1, d, index = c("id", "t"), treat = "first", control = "never"))
```
