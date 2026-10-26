# Plot event-study coefficients

Plot event-study coefficients

## Usage

``` r
# S3 method for class 'panglm_eventstudy'
plot(x, xlab = "Event time", ylab = "Estimate", ...)
```

## Arguments

- x:

  a `"panglm_eventstudy"` object from
  [`eventstudy()`](https://ielbadisy.github.io/panglm/reference/eventstudy.md)

- xlab, ylab:

  axis labels

- ...:

  passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html)

## Value

`x`, invisibly

## See also

Other eventstudy methods:
[`eventstudy-methods`](https://ielbadisy.github.io/panglm/reference/eventstudy-methods.md)

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
```
