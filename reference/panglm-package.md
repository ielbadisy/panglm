# panglm: Generalized Linear Models for Panel Data

Provides generalized linear models for panel data, including pooled,
fixed-effects, and random-effects estimators for continuous, binary,
count, ordinal, positive (gamma), and bounded (beta) outcomes, with the
fixed-effects ordered logit of Baetschmann, Staub and Winkelmann (2015)
[doi:10.1111/rssa.12090](https://doi.org/10.1111/rssa.12090) and beta
regression (Ferrari and Cribari-Neto 2004)
[doi:10.1080/0266476042000214501](https://doi.org/10.1080/0266476042000214501)
. Also offers correlated random effects (Mundlak 1978)
[doi:10.2307/1913646](https://doi.org/10.2307/1913646) , dynamic models
with the Wooldridge (2005)
[doi:10.1002/jae.770](https://doi.org/10.1002/jae.770) initial
conditions solution, observation weights, average partial effects,
regression-spline curves and linearity tests, difference-in-differences
and event studies, and two-part models for costs (Duan et al. 1983)
[doi:10.1080/07350015.1983.10509330](https://doi.org/10.1080/07350015.1983.10509330)
, including correlated random intercepts (Van Ourti 2004)
[doi:10.1002/hec.869](https://doi.org/10.1002/hec.869) . Estimators and
tests follow standard panel-data references, including Hausman (1978)
[doi:10.2307/1913827](https://doi.org/10.2307/1913827) , Chamberlain
(1980) [doi:10.2307/2297110](https://doi.org/10.2307/2297110) , Allison
and Waterman (2002)
[doi:10.1111/1467-9531.00117](https://doi.org/10.1111/1467-9531.00117) ,
and Croissant and Millo (2008)
[doi:10.18637/jss.v027.i02](https://doi.org/10.18637/jss.v027.i02) .
Core numerical routines are implemented using 'Rcpp', 'RcppArmadillo',
and 'RcppParallel'.

## Threads

The within and random-effects estimators for Gaussian, binomial and
Poisson models compute per-group quantities in parallel with
RcppParallel. By default they use all available cores. Set
`options(panglm.threads = k)` to use at most `k` threads; the
environment variable `RCPP_PARALLEL_NUM_THREADS` is also honoured. When
`OMP_THREAD_LIMIT` is set, it caps the number of threads, and under
`R CMD check` (when `_R_CHECK_LIMIT_CORES_` is set) at most 2 threads
are used. Matrix algebra done in R goes through the BLAS library R is
linked to, whose threads are controlled outside panglm (for example with
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` or `MKL_NUM_THREADS`).

## See also

Useful links:

- <https://elbadisyimad.com/panglm/>

- <https://CRAN.R-project.org/package=panglm>

- <https://github.com/ielbadisy/panglm>

- Report bugs at <https://github.com/ielbadisy/panglm/issues>

## Author

**Maintainer**: Imad El Badisy <elbadisyimad@gmail.com>
