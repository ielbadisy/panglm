# panglm 0.1.0

Initial release.

## Estimators

* `model = "pooling"`: gaussian/poisson/binomial GLM via IRLS.
* `model = "within"`: gaussian (exact demeaning), poisson (conditional
  MLE), and binomial (exact conditional logistic regression, Chamberlain
  1980) fixed effects; `effect = "twoways"` adds a time fixed effect for
  gaussian, parallelized in C++ (RcppParallel).
* `model = "random"`: gaussian (Swamy-Arora), poisson (Poisson-Gamma), and
  binomial (Gauss-Hermite quadrature) random effects.

Fixed-effects binomial has no equivalent in `pglm`; it's validated against
`survival::clogit(method = "exact")`.

## Inference

* `vcov(fit, type = "HC1" | "cluster")` for pooled/within models.
* `panglm_hausman()` for FE-vs-RE specification testing.
* `tidy()`/`glance()` methods registered against `generics::tidy`/`glance`
  for `broom`/`modelsummary` compatibility.

## Known limitations

* No between/twoways effects for poisson/binomial.
* No ordinal/tobit families.
* Not yet on CRAN.
