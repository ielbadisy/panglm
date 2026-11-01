# Compare coefficient estimates of several fitted models

Dot-and-whisker plot of the coefficients of one or more fitted models,
terms on the vertical axis and Wald confidence intervals as horizontal
segments. Several models are drawn side by side for each term (dodged
vertically, with distinct plotting symbols and colours and a legend), so
that, for example, pooled, fixed-effects, random-effects and correlated
random effects fits of the same specification can be compared in one
figure. Every term present in at least one model is shown; a model that
does not contain a term has no mark on that row. The left margin is
widened to the longest term label so that no label is clipped.

## Usage

``` r
coefplot(
  ...,
  parm = NULL,
  level = 0.95,
  labels = NULL,
  main = "Coefficient estimates",
  xlab = "Estimate",
  col = NULL,
  pch = NULL,
  ref = 0
)
```

## Arguments

- ...:

  fitted models, or a single unclassed list of fitted models. Argument
  names are used as legend labels.

- parm:

  names or positions (in the order the terms first appear across the
  models) of the terms to plot; by default every term except
  `"(Intercept)"`

- level:

  confidence level of the intervals

- labels:

  legend labels of the models; by default the argument names, or the
  deparsed arguments

- main:

  plot title

- xlab:

  horizontal-axis label

- col:

  colours of the models, recycled

- pch:

  plotting symbols of the models, recycled

- ref:

  position of the dashed reference line; `NULL` for none

## Value

the models, invisibly, as a named list

## Details

Any fitted object with [`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) methods can be used, e.g.
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md),
[`paniv()`](https://elbadisyimad.com/panglm/reference/paniv.md),
[`pangmm()`](https://elbadisyimad.com/panglm/reference/pangmm.md),
[`panfmm()`](https://elbadisyimad.com/panglm/reference/panfmm.md) or
[`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md)
fits, and also [`stats::glm()`](https://rdrr.io/r/stats/glm.html) fits.
When [`vcov()`](https://rdrr.io/r/stats/vcov.html) fails or has no entry
for a term, only the estimate is drawn.

## See also

Other panglm methods:
[`confint.panglm()`](https://elbadisyimad.com/panglm/reference/confint.panglm.md),
[`effectcurve()`](https://elbadisyimad.com/panglm/reference/effectcurve.md),
[`fitted.panglm()`](https://elbadisyimad.com/panglm/reference/fitted.panglm.md),
[`fixef()`](https://elbadisyimad.com/panglm/reference/fixef.md),
[`glance.panglm()`](https://elbadisyimad.com/panglm/reference/glance.panglm.md),
[`logLik.panglm()`](https://elbadisyimad.com/panglm/reference/logLik.panglm.md),
[`model.matrix.panglm()`](https://elbadisyimad.com/panglm/reference/model.matrix.panglm.md),
[`partialeffects.panordered()`](https://elbadisyimad.com/panglm/reference/partialeffects.md),
[`plot.panglm()`](https://elbadisyimad.com/panglm/reference/plot.panglm.md),
[`plot.panglm_ape()`](https://elbadisyimad.com/panglm/reference/plot.panglm_ape.md),
[`plot.panglm_effectcurve()`](https://elbadisyimad.com/panglm/reference/plot.panglm_effectcurve.md),
[`predict.panglm()`](https://elbadisyimad.com/panglm/reference/predict.panglm.md),
[`predict.panordered()`](https://elbadisyimad.com/panglm/reference/predict.panordered.md),
[`residuals.panglm()`](https://elbadisyimad.com/panglm/reference/residuals.panglm.md),
[`summary.panglm()`](https://elbadisyimad.com/panglm/reference/summary.panglm.md),
[`tidy.panglm()`](https://elbadisyimad.com/panglm/reference/tidy.panglm.md),
[`vcov.panglm()`](https://elbadisyimad.com/panglm/reference/vcov.panglm.md)

## Examples

``` r
data(copd)
f <- exacerbations ~ crp + treatment + age
pool <- panglm(f, data = copd, index = c("id", "visit"),
               model = "pooling", family = "poisson")
re <- panglm(f, data = copd, index = c("id", "visit"),
             model = "random", family = "poisson")
cre <- panglm(f, data = copd, index = c("id", "visit"),
              model = "random", family = "poisson", cre = TRUE)
coefplot(Pooled = pool, RE = re, CRE = cre)
```
