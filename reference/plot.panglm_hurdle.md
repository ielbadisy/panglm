# Plot a fitted panglm hurdle model

Displays coefficient estimates and Wald confidence intervals for the
binary hurdle component, the positive-count component, or both.

## Usage

``` r
# S3 method for class 'panglm_hurdle'
plot(x, which = c("both", "zero", "count"), level = 0.95, ...)
```

## Arguments

- x:

  a fitted `"panglm_hurdle"` object

- which:

  one of `"both"`, `"zero"`, or `"count"`

- level:

  confidence level for coefficient intervals

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
```
