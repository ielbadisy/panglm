# Individual effects of a panel model

`fixef()` returns the estimated individual effects `alpha_i` of a
fixed-effects (`model = "within"`) fit, and `ranef()` the empirical
Bayes predictions of the random effects of a `model = "random"` fit.

## Usage

``` r
fixef(object, ...)

ranef(object, ...)

# S3 method for class 'panglm'
fixef(object, ...)

# S3 method for class 'panglm_hurdle'
fixef(object, ...)

# S3 method for class 'panglm'
ranef(object, type = c("mean", "mode"), ...)
```

## Arguments

- object:

  a fitted
  [`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md)
  model (or a
  [`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md)
  fit for `fixef()`)

- ...:

  unused

- type:

  for `ranef()`, `"mean"` (posterior mean, default) or `"mode"`
  (posterior mode)

## Value

`fixef()`: a named numeric vector of individual effects (names are the
individual identifiers), or for two-way fits a list with components
`individual` and `time`. `ranef()`: a matrix with one row per individual
(row names are the identifiers) and one column per random effect, with
attribute `"type"`.

## Fixed effects

- Gaussian within fits: `alpha_i = ybar_i - xbar_i' beta`, the
  intercepts of the least-squares dummy-variable regression. They are
  not normalized:
  [`fixest::fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html)
  gives the same values for a model without an intercept.

- Poisson within fits: the conditional (Poisson fixed-effects) estimator
  does not estimate `alpha_i`, but the individual intercept that
  maximizes the unconditional likelihood at the conditional `beta` has
  the closed form `alpha_i = log(sum_t y_it / sum_t exp(x_it' beta))`,
  the values reported by
  [`fixest::fepois()`](https://lrberge.github.io/fixest/reference/feglm.html).
  Individuals with no positive outcome (`alpha_i = -Inf`) are omitted,
  as in fixest.

- Negative binomial and gamma within fits: the dummy-variable intercepts
  of the joint maximum likelihood fit. The negative binomial fit drops
  the individuals whose outcome is always zero.

- Two-way fits (Gaussian, Poisson, negative binomial): a list with the
  `individual` and `time` effects, normalized so that the time effects
  average zero.

- The conditional logit (`family = "binomial"`, one-way or two-way) and
  the fixed-effects ordered logit condition the individual effects out
  of the likelihood and do not estimate them: `fixef()` stops with an
  error. The individual intercepts are not consistently estimable with a
  fixed number of periods anyway (incidental parameters).

- [`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md)
  fits: the dummy-variable intercepts of the zero-truncated count part
  (individuals with at least one positive count); the binary part is a
  conditional logit.

## Random effects

`ranef()` returns the posterior mean (`type = "mean"`, default) or the
posterior mode (`type = "mode"`) of each individual's random effects
given the individual's outcomes, at the estimated parameters (empirical
Bayes).

- Gaussian random intercept: the BLUP
  `sigma_mu^2 T_i / (sigma_mu^2 T_i + sigma_e^2) * mean_t(y_it - x_it' beta)`
  with the variance components of the fit (Swamy-Arora, or maximum
  likelihood for weighted fits). Mean and mode coincide.

- Gaussian random slopes: the BLUPs stored by the fit.

- Binomial, gamma, beta, Tobit and ordered random intercepts, which are
  normal on the link scale: the posterior mean is computed by adaptive
  Gauss-Hermite quadrature with the `R` nodes of the fit, centred at the
  posterior mode and scaled by the curvature there, the rule used for
  the likelihood. With `R = 1` (Laplace) the mean is the mode.

- Poisson and negative binomial random intercepts (Hausman, Hall and
  Griliches), whose mixing distribution is not normal: the effect is
  reported on the log scale of the multiplicative effect on the mean.
  For the Poisson-gamma model `nu_i ~ Gamma(d, d)` and the posterior is
  `Gamma(d + sum_t y_it, d + sum_t exp(x_it' beta))`, so the posterior
  mean of `log(nu_i)` is `digamma(d + Y_i) - log(d + L_i)` and its mode
  `log((d + Y_i) / (d + L_i))`. For the beta-mixed negative binomial
  model the effect is `log(delta_i)`, with `delta_i / (1 + delta_i)`
  distributed as `Beta(b, a)`, giving
  `digamma(b + Y_i) - digamma(a + L_i)` and
  `log((b + Y_i) / (a + L_i))`.

- Random slopes of the Poisson, negative binomial, gamma and binomial
  families: posterior modes are those stored by the fit; posterior means
  use the adaptive product rule of the fit (`R^q` nodes).

Empirical Bayes predictions are shrunk towards zero, the more so the
fewer observations an individual has: their spread is smaller than the
estimated standard deviation of the random effect (for the Gaussian
model, `var(b_i) = sigma_mu^2 * shrinkage_i`), so their histogram is
narrower than the assumed normal density (see
[`plot.panglm()`](https://elbadisyimad.com/panglm/reference/plot.panglm.md),
`which = "effects"`).

## Generics

The generics `fixef()` and `ranef()` are also defined by nlme, whose
generics lme4 and fixest re-export. panglm exports its own generics so
that it does not depend on nlme, and when nlme is loaded it also
registers its methods for nlme's generics: whichever package is attached
last, `fixef(fit)` and `ranef(fit)` dispatch to the panglm methods.

## See also

[`plot.panglm()`](https://elbadisyimad.com/panglm/reference/plot.panglm.md)
with `which = "effects"` to plot them.

Other panglm methods:
[`coefplot()`](https://elbadisyimad.com/panglm/reference/coefplot.md),
[`confint.panglm()`](https://elbadisyimad.com/panglm/reference/confint.panglm.md),
[`effectcurve()`](https://elbadisyimad.com/panglm/reference/effectcurve.md),
[`fitted.panglm()`](https://elbadisyimad.com/panglm/reference/fitted.panglm.md),
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
fe <- panglm(fev1 ~ crp, data = copd, index = c("id", "visit"),
             model = "within")
head(fixef(fe))
#>        1        2        3        4        5        6 
#> 3.072628 3.011334 2.915205 3.492994 3.210295 2.814040 

re <- panglm(fev1 ~ crp + treatment, data = copd,
             index = c("id", "visit"), model = "random")
head(ranef(re))
#>   (Intercept)
#> 1  -0.1462000
#> 2  -0.2056603
#> 3  -0.2981489
#> 4   0.2606700
#> 5   0.1480610
#> 6  -0.3970215

hosp <- panglm(hospitalized ~ crp + treatment, data = copd,
               index = c("id", "visit"), model = "random",
               family = "binomial")
head(ranef(hosp))
#>          [,1]
#> 1  0.02573862
#> 2  0.02580552
#> 3  0.03808138
#> 4  0.01191119
#> 5  0.56973142
#> 6 -0.43496147
```
