## Update submission (1.1.4 -> 1.7.0)

This update adds new panel estimators (between, Tobit, ordered, beta,
gamma, random slopes, panel IV and Hausman-Taylor, dynamic GMM, sample
selection, quantile regression, latent class models, and staggered
difference-in-differences), tidy methods for every model class, plot
methods for every fitted object, five health economics teaching panels and
task-based vignettes. Exported names with underscores are deprecated
through wrappers that keep working and warn. Estimates of existing models
are unchanged. See NEWS.md for details.

Every estimator is covered by unit tests and by validation scripts kept in
the GitHub repository.

## Test environments

* Local: Ubuntu 24.04, R 4.5.1

## R CMD check results

0 errors | 0 warnings | 1 note

* "Compilation used the following non-portable flag(s):
  '-mno-omit-leaf-frame-pointer'". This flag comes from the local Ubuntu R
  toolchain (`Makeconf`), not from the package's `src/Makevars`.

Parallel code respects `OMP_THREAD_LIMIT` and uses at most 2 threads under
`R CMD check`.

## Reverse dependencies

None on CRAN.
