# Plot a latent class (finite mixture) model

`which = "coefficients"` compares the coefficients of the latent classes
in one dot-and-whisker panel (the classes dodged within each term, as in
[`coefplot()`](https://elbadisyimad.com/panglm/reference/coefplot.md));
the coefficients of the class-membership (concomitant) model are not
shown. `which = "posterior"` is a classification-quality diagnostic: the
histogram of each individual's largest posterior class probability. Mass
close to one means well-separated classes; the subtitle reports its mean
and the relative entropy `1 - sum(-p log p) / (N log k)` (one for a
perfect classification). `which = "correlation"` draws the heatmap of
the correlation matrix of all the estimates, `cov2cor(vcov(x))` (see
[`plot.panglm()`](https://elbadisyimad.com/panglm/reference/plot.panglm.md)).

## Usage

``` r
# S3 method for class 'panfmm'
plot(
  x,
  which = c("coefficients", "posterior", "correlation"),
  parm = NULL,
  level = 0.95,
  cells = FALSE,
  main = NULL,
  xlab = NULL,
  ...
)
```

## Arguments

- x:

  a fitted `"panfmm"` object

- which:

  `"coefficients"` (default), `"posterior"` or `"correlation"`

- parm:

  names or positions of the terms to plot (without the class prefix); by
  default every term except `"(Intercept)"`. For `"correlation"`, names
  as in `coef(x)`.

- level:

  confidence level of the intervals

- cells:

  for `"correlation"`, print the correlation in each cell

- main:

  optional plot title

- xlab:

  optional horizontal-axis label

- ...:

  additional arguments passed to
  [`graphics::plot.default()`](https://rdrr.io/r/graphics/plot.default.html)

## Value

`x`, invisibly

## See also

Other panfmm methods:
[`predict.panfmm()`](https://elbadisyimad.com/panglm/reference/predict.panfmm.md),
[`summary.panfmm()`](https://elbadisyimad.com/panglm/reference/summary.panfmm.md),
[`vcov.panfmm()`](https://elbadisyimad.com/panglm/reference/vcov.panfmm.md)

## Examples

``` r
data(copd)
fit <- panfmm(exacerbations ~ crp + treatment, data = copd,
              index = c("id", "visit"), k = 2, family = "poisson",
              nstart = 2, seed = 1)
plot(fit)

plot(fit, which = "posterior")

plot(fit, which = "correlation", cells = TRUE)
```
