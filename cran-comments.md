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

0 errors | 0 warnings | 1 note

* This is a new submission.
* One NOTE on the local machine only: the Ubuntu R toolchain supplies
  `-mno-omit-leaf-frame-pointer` through its system `Makeconf`. The package
  does not set this flag in `src/Makevars`; it does not appear on other
  platforms.

## Downstream dependencies

None (new package).
