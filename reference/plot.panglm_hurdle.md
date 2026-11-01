# Plot a fitted panglm hurdle model

Displays coefficient estimates and Wald confidence intervals for the
binary hurdle component, the positive-count component, or both, in the
dot-and-whisker layout of
[`plot.panglm()`](https://elbadisyimad.com/panglm/reference/plot.panglm.md).
`which = "correlation"` draws the heatmap of the correlation matrix of
all the estimates (see
[`plot.panglm()`](https://elbadisyimad.com/panglm/reference/plot.panglm.md));
the two parts are estimated separately, so the between-part block is
zero.

## Usage

``` r
# S3 method for class 'panglm_hurdle'
plot(
  x,
  which = c("both", "zero", "count", "correlation"),
  level = 0.95,
  parm = NULL,
  cells = FALSE,
  ...
)
```

## Arguments

- x:

  a fitted `"panglm_hurdle"` object

- which:

  one of `"both"`, `"zero"`, `"count"` or `"correlation"`

- level:

  confidence level for coefficient intervals

- parm:

  names or positions of the terms to plot; by default every term except
  `"(Intercept)"`. For `"correlation"`, names as in `coef(x)` (with the
  `zero:` and `count:` prefixes).

- cells:

  for `"correlation"`, print the correlation in each cell

- ...:

  additional arguments passed to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html)

## Value

`x`, invisibly

## See also

Other panhurdle methods:
[`hurdle-methods`](https://elbadisyimad.com/panglm/reference/hurdle-methods.md)

## Examples

``` r
data(copd)
fit <- panhurdle(exacerbations ~ crp, data = copd,
                     index = c("id", "visit"))
#> 26 group(s) with no within-group outcome variation dropped (all-0 or all-1) -- they carry no information for the conditional likelihood.
plot(fit)

plot(fit, which = "correlation")
```
