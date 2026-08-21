## Submission

This is the first submission of `panglm` to CRAN.

`panglm` provides fast estimation of pooled, fixed-effects (within), and
random-effects generalized linear models for panel data (gaussian,
poisson, binomial, and negative-binomial families), with a compiled
(Rcpp/RcppArmadillo/RcppParallel) numerical core. Covers two-way fixed
effects, exact conditional logistic FE (Chamberlain 1980), a fixed-effects
hurdle model for structural-zero counts, Hausman specification testing,
HC1/cluster-robust sandwich variance, and `broom` (`tidy()`/`glance()`)
compatibility. Every estimator is validated in tests against a reference
implementation (`plm`, `pglm`, `fixest`, `MASS::glm.nb`, `survival::clogit`,
`pscl::hurdle`, or `glmmTMB`, as applicable).

## Test environments

* local: Ubuntu 24.04, R 4.5.1 (via `R CMD check --as-cran`)
* win-builder / R-hub: to be run before submission

## R CMD check results

See below for the latest local `R CMD check --as-cran` result.

## Downstream dependencies

None (new package).
