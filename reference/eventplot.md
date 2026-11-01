# Compare event studies on one plot

Draws the event-time effects of several fits on one axis, offset
horizontally, with pointwise confidence intervals and a legend: for
example the two-way fixed-effects event study of
[`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md),
the Callaway and Sant'Anna estimator of
[`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md) and the
extended two-way fixed-effects estimator of
[`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md). Uniform
bands are not drawn here; use the
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method of a
single fit for them.

## Usage

``` r
eventplot(
  ...,
  labels = NULL,
  min_e = -Inf,
  max_e = Inf,
  col = c("black", "#D55E00", "#0072B2", "#009E73", "#CC79A7"),
  pch = c(19, 17, 15, 18, 25),
  offset = 0.2,
  legend = TRUE,
  xlab = "Event time",
  ylab = "Estimate",
  main = NULL,
  ylim = NULL,
  style = c("classic", "points", "line"),
  ref = NULL,
  onset = NULL,
  onset_label = NA,
  caps = NULL
)
```

## Arguments

- ...:

  fits from
  [`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md),
  [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md) or
  [`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md) (or
  the event-time summaries of the last two). Names of the arguments are
  used as labels.

- labels:

  optional labels, one per fit (default: the argument names, else the
  estimator)

- min_e, max_e:

  event window to plot

- col, pch:

  colours and (filled) plotting symbols, recycled over the fits; the
  hollow version of each symbol is used before treatment

- offset:

  horizontal distance between the fits, in event-time units

- legend:

  logical: add a legend

- xlab, ylab, main, ylim:

  axis labels, title and y-axis limits

- style:

  `"classic"` (default), `"points"` or `"line"` (see Styles); for the
  event-time, cohort and calendar plots

- ref:

  event time of the reference period, marked when the estimate of a fit
  is exactly zero there; the default (`NULL`) uses the omitted period of
  each fit (`-1`, or the `ref` of an
  [`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md)
  fit)

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

## Value

a data frame of the plotted estimates (`model`, `event`, `estimate`,
`lower`, `upper`), invisibly

## Details

The styles are those of
[`plot.panglm_attgt()`](https://elbadisyimad.com/panglm/reference/plot.panglm_attgt.md).
With `style = "classic"` (the default) the points are filled, the
whiskers capped, each reference period is a hollow symbol at zero and
the dashed onset line is at event time `0`, labelled "Intervention".
With `style = "points"` the points are hollow before treatment and
filled after, and the onset line is at `-0.5`. With `style = "line"`
each fit is a coloured line with its pointwise intervals as a
semi-transparent band (the fits are not offset).

## See also

Other policy evaluation:
[`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md),
[`did()`](https://elbadisyimad.com/panglm/reference/did.md),
[`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md),
[`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md),
[`panelplot()`](https://elbadisyimad.com/panglm/reference/panelplot.md),
[`treatmentplot()`](https://elbadisyimad.com/panglm/reference/treatmentplot.md),
[`trendplot()`](https://elbadisyimad.com/panglm/reference/trendplot.md)

## Examples

``` r
set.seed(2)
d <- expand.grid(id = 1:150, t = 1:8)
d$first <- rep(c(4, 6, NA), each = 50)[d$id]
d$y <- rnorm(150)[d$id] + 0.2 * d$t +
  ifelse(!is.na(d$first) & d$t >= d$first, 1, 0) + rnorm(nrow(d))
es <- eventstudy(y ~ 1, d, index = c("id", "t"), treat = "first",
                 leads = 3, lags = 3)
cs <- attgt(y ~ 1, d, index = c("id", "t"), treat = "first")
ew <- etwfe(y ~ 1, d, index = c("id", "t"), treat = "first", control = "never")
eventplot(TWFE = es, CS = cs, ETWFE = ew, min_e = -3, max_e = 3)

eventplot(TWFE = es, ETWFE = ew, min_e = -3, max_e = 3, style = "points")

eventplot(TWFE = es, CS = cs, min_e = -3, max_e = 3, style = "line")
```
