# Package index

## Model fitting

The main estimator for panel GLMs (pooled, within, between and random
effects) and the two models built on it for costs and structural zeros.

- [`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md) :
  Fit a generalized linear model for panel data
- [`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md) :
  Two-part model for panel cost data
- [`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md)
  : Fixed-effects hurdle model for panel count data
- [`panglm-package`](https://elbadisyimad.com/panglm/reference/panglm-package.md)
  : panglm: Generalized Linear Models for Panel Data

## Outcome families

Family objects for the `family` argument of
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md).

- [`gaussianfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  [`poissonfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  [`binomialfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  [`negbinfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  [`gammafamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  [`orderedfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  [`betafamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  [`tobitfamily()`](https://elbadisyimad.com/panglm/reference/panglm-families.md)
  : Family specifications for panglm

## Policy evaluation

Difference-in-differences, event studies and staggered-adoption
estimators.

- [`did()`](https://elbadisyimad.com/panglm/reference/did.md) : Two-way
  fixed-effects difference-in-differences
- [`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md)
  : Event-study regression with a pre-trend test
- [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md) :
  Group-time average treatment effects for staggered adoption
- [`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md) :
  Extended two-way fixed-effects difference-in-differences

## Endogeneity

Instrumental variables, Hausman-Taylor, dynamic panels by GMM and sample
selection.

- [`paniv()`](https://elbadisyimad.com/panglm/reference/paniv.md) :
  Instrumental-variable estimators for linear panel models
- [`pangmm()`](https://elbadisyimad.com/panglm/reference/pangmm.md) :
  Difference and system GMM for dynamic linear panel models
- [`panselect()`](https://elbadisyimad.com/panglm/reference/panselect.md)
  : Sample-selection correction for linear panel models

## Heterogeneity

Quantile regression and latent classes. Random slopes are fitted by
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md) with
`model = "random"` and its `random` argument.

- [`panqr()`](https://elbadisyimad.com/panglm/reference/panqr.md) :
  Quantile regression for panel data with individual effects
- [`panfmm()`](https://elbadisyimad.com/panglm/reference/panfmm.md) :
  Latent class (finite mixture) models for panel data
- [`panfmmselect()`](https://elbadisyimad.com/panglm/reference/panfmmselect.md)
  : Compare latent class panel models with different numbers of classes

## Tests

Specification and hypothesis tests on fitted models.

- [`panhausman()`](https://elbadisyimad.com/panglm/reference/panhausman.md)
  : Hausman specification test (fixed vs. random effects)
- [`pandispersiontest()`](https://elbadisyimad.com/panglm/reference/pandispersiontest.md)
  : Pearson chi-squared/df overdispersion test
- [`jointtest()`](https://elbadisyimad.com/panglm/reference/jointtest.md)
  : Joint Wald test of coefficients
- [`lineartest()`](https://elbadisyimad.com/panglm/reference/lineartest.md)
  : Test whether the effect of a regressor is linear

## Post-estimation methods

Methods for
[`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md) fits,
then for the other model classes.

### panglm fits

- [`summary(`*`<panglm>`*`)`](https://elbadisyimad.com/panglm/reference/summary.panglm.md)
  : Summarize a panglm fit
- [`vcov(`*`<panglm>`*`)`](https://elbadisyimad.com/panglm/reference/vcov.panglm.md)
  : Recompute the coefficient covariance matrix for a panglm fit
- [`confint(`*`<panglm>`*`)`](https://elbadisyimad.com/panglm/reference/confint.panglm.md)
  : Confidence intervals for panglm coefficients
- [`predict(`*`<panglm>`*`)`](https://elbadisyimad.com/panglm/reference/predict.panglm.md)
  : Predict from a fitted panglm model
- [`predict(`*`<panordered>`*`)`](https://elbadisyimad.com/panglm/reference/predict.panordered.md)
  : Predict from an ordered panel model
- [`fitted(`*`<panglm>`*`)`](https://elbadisyimad.com/panglm/reference/fitted.panglm.md)
  : Extract fitted values from a panglm model
- [`logLik(`*`<panglm>`*`)`](https://elbadisyimad.com/panglm/reference/logLik.panglm.md)
  : Extract the model log-likelihood
- [`model.matrix(`*`<panglm>`*`)`](https://elbadisyimad.com/panglm/reference/model.matrix.panglm.md)
  : Design matrix of a panglm fit
- [`plot(`*`<panglm>`*`)`](https://elbadisyimad.com/panglm/reference/plot.panglm.md)
  : Plot a fitted panglm model
- [`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md)
  : Average partial effects
- [`effectcurve()`](https://elbadisyimad.com/panglm/reference/effectcurve.md)
  : Estimated effect curve of a continuous regressor
- [`plot(`*`<panglm_effectcurve>`*`)`](https://elbadisyimad.com/panglm/reference/plot.panglm_effectcurve.md)
  : Plot an effect curve

### Other model classes

- [`predict(`*`<twopart>`*`)`](https://elbadisyimad.com/panglm/reference/predict.twopart.md)
  : Predict from a two-part model
- [`confint(`*`<twopart>`*`)`](https://elbadisyimad.com/panglm/reference/confint.twopart.md)
  : Confidence intervals for a two-part model
- [`coef(`*`<panglm_hurdle>`*`)`](https://elbadisyimad.com/panglm/reference/hurdle-methods.md)
  [`vcov(`*`<panglm_hurdle>`*`)`](https://elbadisyimad.com/panglm/reference/hurdle-methods.md)
  [`nobs(`*`<panglm_hurdle>`*`)`](https://elbadisyimad.com/panglm/reference/hurdle-methods.md)
  [`confint(`*`<panglm_hurdle>`*`)`](https://elbadisyimad.com/panglm/reference/hurdle-methods.md)
  [`summary(`*`<panglm_hurdle>`*`)`](https://elbadisyimad.com/panglm/reference/hurdle-methods.md)
  [`print(`*`<summary.panglm_hurdle>`*`)`](https://elbadisyimad.com/panglm/reference/hurdle-methods.md)
  : Coefficient methods for the hurdle model
- [`plot(`*`<panglm_hurdle>`*`)`](https://elbadisyimad.com/panglm/reference/plot.panglm_hurdle.md)
  : Plot a fitted panglm hurdle model
- [`vcov(`*`<panglm_eventstudy>`*`)`](https://elbadisyimad.com/panglm/reference/eventstudy-methods.md)
  [`nobs(`*`<panglm_eventstudy>`*`)`](https://elbadisyimad.com/panglm/reference/eventstudy-methods.md)
  [`confint(`*`<panglm_eventstudy>`*`)`](https://elbadisyimad.com/panglm/reference/eventstudy-methods.md)
  [`summary(`*`<panglm_eventstudy>`*`)`](https://elbadisyimad.com/panglm/reference/eventstudy-methods.md)
  [`print(`*`<summary.panglm_eventstudy>`*`)`](https://elbadisyimad.com/panglm/reference/eventstudy-methods.md)
  [`logLik(`*`<panglm_eventstudy>`*`)`](https://elbadisyimad.com/panglm/reference/eventstudy-methods.md)
  : Methods for event-study fits
- [`plot(`*`<panglm_eventstudy>`*`)`](https://elbadisyimad.com/panglm/reference/plot.panglm_eventstudy.md)
  : Plot event-study coefficients
- [`summary(`*`<panglm_attgt>`*`)`](https://elbadisyimad.com/panglm/reference/summary.panglm_attgt.md)
  : Aggregate group-time average treatment effects
- [`plot(`*`<panglm_attgt>`*`)`](https://elbadisyimad.com/panglm/reference/plot.panglm_attgt.md)
  [`plot(`*`<summary.panglm_attgt>`*`)`](https://elbadisyimad.com/panglm/reference/plot.panglm_attgt.md)
  : Plot an event study from group-time effects
- [`summary(`*`<panglm_etwfe>`*`)`](https://elbadisyimad.com/panglm/reference/summary.panglm_etwfe.md)
  : Average treatment effects on the treated from an ETWFE fit
- [`plot(`*`<panglm_etwfe>`*`)`](https://elbadisyimad.com/panglm/reference/plot.panglm_etwfe.md)
  [`plot(`*`<summary.panglm_etwfe>`*`)`](https://elbadisyimad.com/panglm/reference/plot.panglm_etwfe.md)
  : Plot ETWFE treatment effects
- [`summary(`*`<paniv>`*`)`](https://elbadisyimad.com/panglm/reference/summary.paniv.md)
  : Summarize a paniv fit
- [`vcov(`*`<paniv>`*`)`](https://elbadisyimad.com/panglm/reference/vcov.paniv.md)
  : Covariance of a paniv fit
- [`predict(`*`<paniv>`*`)`](https://elbadisyimad.com/panglm/reference/predict.paniv.md)
  : Predictions from a panel instrumental-variables fit
- [`summary(`*`<pangmm>`*`)`](https://elbadisyimad.com/panglm/reference/summary.pangmm.md)
  : Summarize a pangmm fit
- [`vcov(`*`<pangmm>`*`)`](https://elbadisyimad.com/panglm/reference/vcov.pangmm.md)
  : Covariance of a pangmm fit
- [`predict(`*`<pangmm>`*`)`](https://elbadisyimad.com/panglm/reference/predict.pangmm.md)
  : Predictions from a dynamic panel GMM fit
- [`summary(`*`<panselect>`*`)`](https://elbadisyimad.com/panglm/reference/summary.panselect.md)
  : Summarize a panselect fit
- [`vcov(`*`<panselect>`*`)`](https://elbadisyimad.com/panglm/reference/vcov.panselect.md)
  : Covariance of a panselect fit
- [`predict(`*`<panselect>`*`)`](https://elbadisyimad.com/panglm/reference/predict.panselect.md)
  : Predictions from a panel sample-selection fit
- [`summary(`*`<panqr>`*`)`](https://elbadisyimad.com/panglm/reference/summary.panqr.md)
  : Summarize a panqr fit
- [`vcov(`*`<panqr>`*`)`](https://elbadisyimad.com/panglm/reference/vcov.panqr.md)
  : Covariance of a panqr fit
- [`plot(`*`<panqr>`*`)`](https://elbadisyimad.com/panglm/reference/plot.panqr.md)
  : Plot a panel quantile regression fit
- [`summary(`*`<panfmm>`*`)`](https://elbadisyimad.com/panglm/reference/summary.panfmm.md)
  : Summarize a panfmm fit
- [`vcov(`*`<panfmm>`*`)`](https://elbadisyimad.com/panglm/reference/vcov.panfmm.md)
  : Covariance of a panfmm fit
- [`predict(`*`<panfmm>`*`)`](https://elbadisyimad.com/panglm/reference/predict.panfmm.md)
  : Predictions from a latent class panel model

### Tidiers

- [`tidy(`*`<panglm>`*`)`](https://elbadisyimad.com/panglm/reference/tidy.panglm.md)
  : Tidy a panglm fit into a data frame of coefficients
- [`glance(`*`<panglm>`*`)`](https://elbadisyimad.com/panglm/reference/glance.panglm.md)
  : One-row model summary of a panglm fit
- [`tidy(`*`<panglm_hurdle>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`glance(`*`<panglm_hurdle>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`tidy(`*`<twopart>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`glance(`*`<twopart>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`tidy(`*`<pangmm>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`glance(`*`<pangmm>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`tidy(`*`<paniv>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`glance(`*`<paniv>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`tidy(`*`<panqr>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`glance(`*`<panqr>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`tidy(`*`<panselect>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`glance(`*`<panselect>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`tidy(`*`<panfmm>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`glance(`*`<panfmm>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`tidy(`*`<panglm_attgt>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`glance(`*`<panglm_attgt>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`tidy(`*`<panglm_etwfe>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`glance(`*`<panglm_etwfe>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`tidy(`*`<panglm_eventstudy>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`glance(`*`<panglm_eventstudy>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`tidy(`*`<panglm_hausman>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`glance(`*`<panglm_hausman>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`tidy(`*`<panglm_jointtest>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`glance(`*`<panglm_jointtest>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`tidy(`*`<panglm_dispersiontest>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  [`glance(`*`<panglm_dispersiontest>`*`)`](https://elbadisyimad.com/panglm/reference/tidy-methods.md)
  : Tidy and glance methods for the panglm model classes

## Data

- [`copd`](https://elbadisyimad.com/panglm/reference/copd.md) :
  Synthetic COPD follow-up panel

## Deprecated

Former names, kept as wrappers that warn and call the new function.

- [`gaussian_family()`](https://elbadisyimad.com/panglm/reference/panglm-deprecated.md)
  [`binomial_family()`](https://elbadisyimad.com/panglm/reference/panglm-deprecated.md)
  [`poisson_family()`](https://elbadisyimad.com/panglm/reference/panglm-deprecated.md)
  [`negbin_family()`](https://elbadisyimad.com/panglm/reference/panglm-deprecated.md)
  [`gamma_family()`](https://elbadisyimad.com/panglm/reference/panglm-deprecated.md)
  [`beta_family()`](https://elbadisyimad.com/panglm/reference/panglm-deprecated.md)
  [`ordered_family()`](https://elbadisyimad.com/panglm/reference/panglm-deprecated.md)
  [`tobit_family()`](https://elbadisyimad.com/panglm/reference/panglm-deprecated.md)
  [`panglm_dispersiontest()`](https://elbadisyimad.com/panglm/reference/panglm-deprecated.md)
  [`panglm_hausman()`](https://elbadisyimad.com/panglm/reference/panglm-deprecated.md)
  [`panglm_hurdle()`](https://elbadisyimad.com/panglm/reference/panglm-deprecated.md)
  : Deprecated functions in panglm
