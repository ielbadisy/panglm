# Validation against established packages

Every estimator in panglm is checked against an established R package,
or against a hand-coded reference when no package fits the same model.
The checks live in the [validation
folder](https://github.com/ielbadisy/panglm/tree/main/validation) of the
repository, one script per model family. Each script fits the same model
on the same data with panglm and with the reference, then compares
coefficients, standard errors, log-likelihoods, test statistics and
predictions value by value.

The last full run used panglm 1.6.0 and R 4.5.1. It made 109,739
comparisons and all of them passed. Each comparison has its own
tolerance, stated in the script: from `1e-12` to `1e-6` for closed-form
estimators, and looser where the reference uses numerical integration, a
different optimizer or the bootstrap.

| Script | What it checks | Reference | Comparisons |
|----|----|----|---:|
| `validate-pooled.R` | Pooled Gaussian, Poisson, logit, probit, Gamma and negative binomial | stats, MASS, plm, fixest, sandwich | 255 |
| `validate-within.R` | Fixed effects, one-way and two-way, linear and nonlinear | plm, fixest, pglm, survival, MASS, sandwich | 142 |
| `validate-between.R` | Between estimator | plm, sandwich | 38 |
| `validate-twowaylogit.R` | Two-way fixed effects logit (conditional) | survival, sandwich, brute force | 85 |
| `validate-slopes.R` | Random intercepts and slopes (Laplace, quadrature) | lme4, nlme, glmmTMB, GLMMadaptive, clubSandwich | 5,748 |
| `validate-ordered.R` | Ordered logit and probit: pooled, random and fixed effects | MASS, ordinal, survival | 87 |
| `validate-tobit.R` | Pooled and random effects Tobit | AER, censReg, sandwich | 58 |
| `validate-paniv.R` | Panel instrumental variables | AER, fixest, plm | 2,926 |
| `validate-hausmantaylor.R` | Hausman-Taylor | plm | 84 |
| `validate-gmm.R` | Difference and system GMM | plm | 990 |
| `validate-select.R` | Sample selection with Mundlak terms | sampleSelection, plm, sandwich | 36,451 |
| `validate-qr.R` | Quantile regression (Canay, pooled, MMQR) | quantreg, plm | 174 |
| `validate-fmm.R` | Finite mixtures and selection of the number of classes | flexmix, MASS, sampleSelection, sandwich | 40,224 |
| `validate-did.R` | Staggered adoption: [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md) and [`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md) | did, etwfe | 2,119 |
| `validate-did-weights.R` | The same with sampling weights | did, etwfe | 2,054 |
| `validate-did-boot.R` | Multiplier bootstrap for [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md) | did, analytic standard errors | 333 |
| `validate-did-gaps.R` | DiD and event studies with gaps and nonlinear outcomes | etwfe, fixest, survival, sandwich | 2,327 |
| `validate-methods.R` | Post-estimation methods across model classes | AER, did, etwfe, plm, quantreg, sampleSelection, generics | 15,644 |

## Running the suite

Clone the repository, install panglm and the reference packages, then
run

``` sh
Rscript validation/run-all.R
```

from the package root. Each script writes a CSV of its comparisons to
`validation/results/` and stops with an error if any comparison fails.
The difference-in-differences scripts take the longest, because the
reference implementations compute derivatives numerically.
