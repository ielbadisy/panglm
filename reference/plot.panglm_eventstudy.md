# Plot event-study coefficients

Plots the event-time coefficients of
[`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md)
with pointwise confidence intervals, in the styles of
[`plot.panglm_attgt()`](https://elbadisyimad.com/panglm/reference/plot.panglm_attgt.md)
(see its Styles section). The omitted reference period (zero by
construction) is drawn as an open symbol at zero, without an interval.
With `bin = TRUE` the end points pool all earlier (later) event times.

## Usage

``` r
# S3 method for class 'panglm_eventstudy'
plot(
  x,
  min_e = -Inf,
  max_e = Inf,
  xlab = "Event time",
  ylab = "Estimate",
  main = NULL,
  ylim = NULL,
  legend = FALSE,
  style = c("classic", "points", "line"),
  ref = x$ref,
  onset = NULL,
  onset_label = NA,
  caps = NULL,
  ...
)
```

## Arguments

- x:

  a `"panglm_eventstudy"` object from
  [`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md)

- min_e, max_e:

  event window to plot

- xlab, ylab, main:

  axis labels and title

- ylim:

  y-axis limits

- legend:

  logical: add a compact legend of the symbols

- style:

  `"classic"` (default), `"points"` or `"line"` (see Styles); for the
  event-time, cohort and calendar plots

- ref:

  event time of the reference period, marked when its estimate is
  exactly zero; the default is the period omitted by
  [`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md)
  (`x$ref`)

- onset:

  position of the dashed treatment-onset line on the event time axis:
  `NULL` gives `0` for `style = "classic"` and `-0.5` (between the last
  pre-treatment and the first treated period) otherwise

- onset_label:

  text written above the onset line; the default (`NA`) is
  "Intervention" for `style = "classic"` and none otherwise, `NULL`
  suppresses it

- caps:

  logical: draw caps on the confidence-interval whiskers; the default
  (`NULL`) draws them for `style = "classic"` only

- ...:

  passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html)

## Value

`x`, invisibly

## See also

Other eventstudy methods:
[`eventstudy-methods`](https://elbadisyimad.com/panglm/reference/eventstudy-methods.md)

## Examples

``` r
set.seed(1)
d <- expand.grid(id = 1:100, year = 2015:2020)
d$first <- ifelse(d$id <= 50, 2018, NA)
e <- d$year - d$first
d$visits <- 3 + rnorm(100)[d$id] + ifelse(!is.na(e) & e >= 0, 1, 0) +
  rnorm(nrow(d))
es <- eventstudy(visits ~ 1, data = d, index = c("id", "year"),
                 treat = "first", leads = 3, lags = 2)
plot(es, main = "Effect of the policy by event time")

plot(es, min_e = -2, legend = TRUE)

plot(es, style = "points")

plot(es, style = "line")
```
