## Submission

This is an update of panglm from 1.1.4 to 1.6.0. It adds new panel
estimators (between, Tobit, ordered, beta, gamma, random slopes, panel IV
and Hausman-Taylor, dynamic GMM, sample selection, quantile regression,
latent class models, and staggered difference-in-differences), tidy methods
for every model class, and task-based vignettes. Exported names with
underscores are deprecated through wrappers that keep working and warn. See
NEWS.md for details.

Every estimator is covered by unit tests and by validation scripts kept in
the GitHub repository.

## Test environments

* local: Ubuntu 24.04, R 4.5.1, `R CMD check --as-cran` on a fresh clone

## R CMD check results

0 errors | 0 warnings | 1 note

* NOTE: compilation used `-mno-omit-leaf-frame-pointer`. This flag comes
  from the local Ubuntu R toolchain (`Makeconf`), not from the package's
  `src/Makevars`.

Parallel code respects `OMP_THREAD_LIMIT` and uses at most 2 threads under
`R CMD check`; examples, tests and vignettes run with a CPU/elapsed ratio
close to 1.

## Downstream dependencies

None.
