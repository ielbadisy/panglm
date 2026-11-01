# Synthetic COPD follow-up panel with costs and a staggered programme

A simulated longitudinal (panel) dataset of chronic obstructive
pulmonary disease (COPD) patients followed over five clinic visits, used
throughout this package's examples and vignettes. It is entirely
synthetic (not real patient data), generated to have a realistic panel
structure (time-invariant baseline covariates plus individual-level
frailties or random intercepts). Because every data-generating value is
known (see Details), it is the dataset on which students can check that
an estimator recovers the truth before using it on the real panels
shipped with the package
([docvisits](https://elbadisyimad.com/panglm/reference/docvisits.md),
[meps](https://elbadisyimad.com/panglm/reference/meps.md),
[medicaid](https://elbadisyimad.com/panglm/reference/medicaid.md),
[healthspend](https://elbadisyimad.com/panglm/reference/healthspend.md),
[cigar](https://elbadisyimad.com/panglm/reference/cigar.md)). The
generating code is in `data-raw/copd.R` in the package's GitHub
repository.

## Usage

``` r
copd
```

## Format

A data frame with 300 rows (60 patients x 5 visits) and 16 columns; the
panel index is `c("id", "visit")`:

- id:

  integer patient identifier (1-60), the panel individual index

- visit:

  integer visit number (1-5), the panel time index

- treatment:

  integer treatment arm, time-invariant (0 = control, 1 = active)

- age:

  numeric patient age in whole years at baseline, time-invariant

- smoker:

  integer current-smoker indicator, time-invariant (0/1)

- bmi:

  numeric body mass index (kg/m^2) at baseline, one decimal,
  time-invariant; not used to generate the outcomes

- crp:

  numeric C-reactive protein (mg/L) measured at each visit, one decimal
  and at least 0.1; the one time-varying covariate, needed for
  fixed-effects ("within") fits to have something left to estimate after
  demeaning

- fev1:

  numeric forced expiratory volume in 1 second (litres), two decimals
  and at least 0.5; a continuous lung-function outcome, one value per
  visit

- exacerbations:

  integer number of COPD exacerbations since the previous visit; a count
  outcome, one value per visit

- hospitalized:

  integer indicator (0/1) of a hospital admission in the interval since
  the previous visit; a binary outcome, one value per visit

- clinic:

  integer clinic (1-6) treating the patient, time-invariant; 10 patients
  per clinic

- adopt:

  integer visit at which the patient's clinic adopted the
  care-management programme (3 or 4), 0 if the clinic never adopted; the
  `treat` column for
  [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md),
  [`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md) and
  [`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md)

- programme:

  integer indicator (0/1) that the programme is in place at this visit
  (`adopt > 0` and `visit >= adopt`)

- cost:

  numeric outpatient plus inpatient cost (euros, two decimals) since the
  previous visit; 28.7 percent zeros and a right-skewed positive part

- adherence:

  numeric proportion of prescribed inhaler doses taken since the
  previous visit, in (0, 1), three decimals

- symptoms:

  ordered factor of symptom severity at the visit: none \< mild \<
  moderate \< severe

## Source

Simulated; see `data-raw/copd.R` in the package's GitHub repository.

## Details

The first ten columns (`id` to `hospitalized`) are unchanged from
earlier releases. The last six (`clinic` to `symptoms`) add health
economics outcomes and were generated afterwards in a separate random
number stream, so the clinical columns did not change.

Data-generating values (patient random intercepts are normal with mean
zero and are independent of the covariates; `visit` enters linearly).

- `exacerbations`: Poisson with mean
  `f_i exp(0.4 - 0.35 treatment + 0.02 (age - 65) + 0.3 smoker + 0.05 crp)`,
  with a gamma frailty `f_i` of mean 1 and variance 0.25.

- `hospitalized`: logit with index
  `-1.5 - 0.6 treatment + 0.03 (age - 65) + 0.5 smoker + 0.08 crp + r_i`,
  `sd(r_i) = 0.4`.

- `fev1`:
  `3.2 - 0.02 (age - 65) - 0.25 smoker + 0.15 treatment - 0.03 crp + r_i + e`,
  `sd(e) = 0.2` (same `r_i`).

- `cost`, any cost: probit with index
  `-0.6 + 0.4 exacerbations + 1.5 hospitalized + 0.3 smoker - 0.2 treatment + 0.05 crp + a_i`,
  `sd(a_i) = 0.5`.

- `cost`, positive part: lognormal, \`log(cost) = 5.8 + 0.1
  exacerbations

  - 0.7 hospitalized + 0.25 smoker - 0.2 treatment + 0.04 crp + 0.05
    visit

  &nbsp;

  - 0.5 programme + c_k + b_i + e`, with `sd(b_i) = 0.4`, `corr(a_i,
    b_i) = 0.5`, `sd(e) = 0.4`and clinic effects`c_k = (0.20, 0.10,
    0.15, 0.00, -0.10,
    -0.15)`for clinics 1 to 6. The programme therefore lowers expected cost among users by`1 -
    exp(-0.5)\`, about 39 percent, without changing the probability of
    any cost. Clinics 1 and 2 adopt at visit 3, clinics 3 and 4 at visit
    4, clinics 5 and 6 never; adopting clinics have higher baseline
    costs (selection on levels) but parallel trends hold. Because both
    potential costs were simulated from the same draws, the sample
    average effect on the treated rows is known: about -406.0 euros per
    visit interval (100 treated patient-visits).

- `adherence`: beta with precision 25 and logit mean
  `0.8 + 0.6 treatment - 0.5 smoker - 0.08 visit + v_i`,
  `sd(v_i) = 0.3`.

- `symptoms`: ordered logit with latent index
  `0.5 smoker - 0.6 treatment + 0.25 crp + w_i`, `sd(w_i) = 0.8`, and
  cut points 0.5, 2.0 and 3.5.

Random-effects (conditional) estimators target these values; pooled
nonlinear fits estimate population-averaged coefficients, which are
attenuated for the binary, ordered and beta outcomes. With only 60
patients the estimates are noisy: `data-raw/copd.R` also refits the
correctly specified models on a 6,000-patient panel drawn from the same
design and checks that they land close to these values.

## Teaching use

- Lung function and inflammation (`fev1` on `crp`): pooled, within,
  between and random effects, CRE and
  [`panhausman()`](https://elbadisyimad.com/panglm/reference/panhausman.md).

- Exacerbation counts: `family = "poisson"` or `"negbin"`,
  [`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md),
  [`panfmm()`](https://elbadisyimad.com/panglm/reference/panfmm.md).

- Hospital admission: `family = "binomial"` or `"probit"`, with
  `dynamic = TRUE` for state dependence.

- What drives outpatient and inpatient cost?
  [`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md)
  (probit plus lognormal or gamma, `correlated = TRUE`),
  `family = "gamma"` on positive costs, `family = "tobit"`,
  [`panqr()`](https://elbadisyimad.com/panglm/reference/panqr.md) for
  cost quantiles.

- Did the care-management programme lower costs?
  [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md),
  [`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md),
  [`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md)
  with `treat = "adopt"`, against the known effect above.

- Adherence (`family = "beta"`) and symptom severity
  (`family = "ologit"`, within = fixed-effects ordered logit).

## Examples

``` r
data(copd)
str(copd)
#> 'data.frame':    300 obs. of  16 variables:
#>  $ id           : int  1 1 1 1 1 2 2 2 2 2 ...
#>  $ visit        : int  1 2 3 4 5 1 2 3 4 5 ...
#>  $ treatment    : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ age          : num  70 70 70 70 70 64 64 64 64 64 ...
#>  $ smoker       : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ bmi          : num  28.5 28.5 28.5 28.5 28.5 24 24 24 24 24 ...
#>  $ crp          : num  5.9 0.1 4.1 5.1 3.3 4.2 4.4 1.4 4.9 3.7 ...
#>  $ fev1         : num  2.61 3.19 2.88 3.12 2.91 2.81 3.01 3.04 2.74 2.8 ...
#>  $ exacerbations: int  2 3 5 1 3 0 1 1 2 1 ...
#>  $ hospitalized : int  0 0 1 0 0 0 0 0 0 1 ...
#>  $ clinic       : int  6 6 6 6 6 4 4 4 4 4 ...
#>  $ adopt        : int  0 0 0 0 0 4 4 4 4 4 ...
#>  $ programme    : int  0 0 0 0 0 0 0 0 1 1 ...
#>  $ cost         : num  386 0 717 0 455 ...
#>  $ adherence    : num  0.891 0.889 0.716 0.834 0.708 0.762 0.787 0.856 0.785 0.769 ...
#>  $ symptoms     : Ord.factor w/ 4 levels "none"<"mild"<..: 2 2 1 2 2 4 3 1 2 3 ...
mean(copd$cost == 0)
#> [1] 0.2866667
# two-part model of cost, programme effect on the positive part
# (clinic dummies absorb the clinic cost levels that differ between
# adopting and never-adopting clinics)
tp <- twopart(cost ~ exacerbations + hospitalized + smoker + treatment +
                crp + visit + programme + factor(clinic),
              binary_formula = ~ exacerbations + hospitalized + smoker +
                treatment + crp,
              data = copd, index = c("id", "visit"), positive = "lognormal")
summary(tp)
#> 
#> Call:
#> twopart(formula = cost ~ exacerbations + hospitalized + smoker + 
#>     treatment + crp + visit + programme + factor(clinic), data = copd, 
#>     index = c("id", "visit"), binary_formula = ~exacerbations + 
#>         hospitalized + smoker + treatment + crp, positive = "lognormal")
#> 
#> Estimator: two-part model (pooling)
#> N = 300  groups = 60  positive = 214  vcov: classical 
#> 
#> Part 1: P(y > 0), binomial probit 
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   -0.54237    0.31625  -1.715   0.0863 .  
#> exacerbations  0.36359    0.07030   5.172 2.32e-07 ***
#> hospitalized   1.77472    0.35753   4.964 6.91e-07 ***
#> smoker         0.27834    0.19551   1.424   0.1546    
#> treatment      0.13107    0.19405   0.675   0.4994    
#> crp            0.01649    0.04424   0.373   0.7093    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Part 2: E[y | y > 0], linear model for log(y) 
#>                 Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)      6.21068    0.17872  34.750  < 2e-16 ***
#> exacerbations    0.07996    0.02057   3.888 0.000101 ***
#> hospitalized     0.76694    0.08127   9.437  < 2e-16 ***
#> smoker           0.26218    0.08430   3.110 0.001871 ** 
#> treatment       -0.32499    0.08399  -3.870 0.000109 ***
#> crp              0.04715    0.01932   2.440 0.014672 *  
#> visit            0.03582    0.03874   0.925 0.355133    
#> programme       -0.51788    0.13450  -3.850 0.000118 ***
#> factor(clinic)2 -0.34519    0.13484  -2.560 0.010466 *  
#> factor(clinic)3  0.01610    0.14289   0.113 0.910264    
#> factor(clinic)4 -0.22703    0.14696  -1.545 0.122401    
#> factor(clinic)5 -0.21778    0.16107  -1.352 0.176358    
#> factor(clinic)6 -0.23223    0.15280  -1.520 0.128539    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Retransformation (smearing): 1.158
#> Log-likelihood: -1734 
#> 
# staggered adoption: group-time effects of the programme on cost
att <- attgt(cost ~ 1, data = copd, index = c("id", "visit"),
             treat = "adopt", control = "notyet")
summary(att)
#> 
#> Callaway and Sant'Anna group-time ATT, aggregated overall (simple weighted average)
#> Control group: not yet treated  Estimator: unconditional  Individuals: 60 
#> 
#> ATT: -184 (SE 199.1, 95% CI -574.2 to 206.2)
#> 
#> Pre-trend test (pre-treatment ATT(g, t) = 0): Chi-squared = 2.116, df = 3, p-value = 0.5487
#> 
```
