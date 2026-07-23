# panglm 0.1.0

Initial release.

## Estimators

* `model = "pooling"`: gaussian/poisson/binomial GLM via IRLS.
* `model = "within"`: gaussian (exact demeaning) and poisson (conditional
  MLE) fixed effects; `effect = "twoways"` adds a time fixed effect for
  gaussian.
* `model = "random"`: gaussian (Swamy-Arora), poisson (Poisson-Gamma), and
  binomial (Gauss-Hermite quadrature) random effects.

## Inference

* `vcov(fit, type = "HC1" | "cluster")` for pooled/within models.
* `panglm_hausman()` for FE-vs-RE specification testing.
* `tidy()`/`glance()` methods registered against `generics::tidy`/`glance`
  for `broom`/`modelsummary` compatibility.

## Known limitations

* No between/twoways effects for poisson/binomial.
* No fixed-effects (conditional) binomial.
* No ordinal/tobit families.
* Not yet on CRAN.
