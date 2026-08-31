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
* win-builder: R-devel and R-release

## R CMD check results

0 errors | 0 warnings | 3 notes

* This is a new submission.
* NOTE: "unable to verify current time" / future file timestamps. This is a
  clock artifact of the local check environment (no NTP access); it does not
  reflect anything in the package.
* NOTE: compilation used `-mno-omit-leaf-frame-pointer`. This flag is
  injected by the local Ubuntu R toolchain's system `Makeconf`. The package
  does not set it in `src/Makevars`; it does not appear on other platforms.
* The installed size is ~12 Mb, almost entirely the compiled
  Rcpp/RcppArmadillo/RcppParallel shared object in `libs/`.

## Downstream dependencies

None (new package).

## Downstream dependencies

None (new package).
