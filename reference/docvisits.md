# German health-care demand panel (doctor and hospital visits, GSOEP)

A real, public, unbalanced panel of 7,293 adults (aged 25 to 64) from
the German Socio-Economic Panel (GSOEP), observed in up to seven waves
(1984 to 1988, 1991 and 1994), 27,326 person-years in total. It is the
data set of Riphahn, Wambach and Million (2003), who studied doctor and
hospital visits in relation to public and add-on private health
insurance, and it is the "German health care data" used in Greene's
Econometric Analysis and his discrete choice lecture notes. The outcomes
cover the main types of health-econometrics panel models: a count with
many zeros and a long right tail (`docvis`, `hospvis`), binary use
indicators (`doctor`, `hospital`) and an 11-point ordinal
self-assessment (`hsat`). The unit is the person and the time index is
the calendar year; the panel index is `c("id", "year")`. The build
script, including the checks against published numbers, is
`data-raw/docvisits.R` in the package's GitHub repository.

## Usage

``` r
docvisits
```

## Format

A data frame with 27326 rows (7,293 persons, 1 to 7 waves each) and 21
columns. Panel index: `id` (individual) and `year` (time).

- id:

  integer person identifier (1 to 7293), the panel individual index

- year:

  integer calendar year of the interview (1984, 1985, 1986, 1987, 1988,
  1991, 1994), the panel time index

- docvis:

  integer number of doctor visits in the last three months (0 to 121);
  count outcome

- hospvis:

  integer number of hospital visits in the last calendar year (0 to 51);
  count outcome

- doctor:

  integer indicator (0/1) of at least one doctor visit, `docvis > 0`;
  binary outcome

- hospital:

  integer indicator (0/1) of at least one hospital visit, `hospvis > 0`;
  binary outcome

- hsat:

  integer health satisfaction, 0 (low) to 10 (high); ordinal outcome

- female:

  integer indicator, 1 = female, 0 = male; time-invariant

- age:

  integer age in years (25 to 64)

- educ:

  numeric years of schooling (7 to 18); constant within person for
  almost everyone

- married:

  integer indicator, 1 = married

- hhkids:

  integer indicator, 1 = children under 16 in the household

- hhninc:

  numeric household nominal monthly net income, in 1000 German marks
  (DM); 4 person-years report 0

- working:

  integer indicator, 1 = employed

- bluecollar:

  integer indicator, 1 = blue-collar employee

- whitecollar:

  integer indicator, 1 = white-collar employee

- self:

  integer indicator, 1 = self-employed

- civil:

  integer indicator, 1 = civil servant (Beamter)

- handicap:

  integer indicator, 1 = officially handicapped (1987 values recoded,
  see Description)

- public:

  integer indicator, 1 = insured in the public (statutory) health
  insurance

- addon:

  integer indicator, 1 = holds add-on (supplementary private) health
  insurance

## Source

Journal of Applied Econometrics Data Archive, Riphahn, Wambach and
Million (2003), file `rwm-data.zip`:
<http://qed.econ.queensu.ca/jae/2003-v18.4/riphahn-wambach-million/>
(now at <https://journaldata.zbw.eu/>, DOI
10.15456/jae.2022314.1312212303). The original survey is the German
Socio-Economic Panel (SOEP, DIW Berlin). The same 27,326 x 25 file is
distributed on CRAN under the GPL (\>= 2) as `Health` in package Rchoice
(citing Riphahn et al. 2003), and subsets are in package COUNT (GPL-2)
as `rwm`, `rwm5yr` and `rwm1984`; the build script checks that the JAE
file and the Rchoice copy are identical. The JAE readme states that "the
distributors of the data require that the data are only used for the
purpose of scientific studies": use these data for teaching and research
only, and cite Riphahn, Wambach and Million (2003).

## Details

Cleaning relative to the original file: the 40 imputed non-integer
values of `hsat` (between 6.5 and 7) are set to 7, as in Greene (2008);
the handicap indicator, which the data authors report as swapped (0
and 1) in 1987, is recoded for that year; the degree of handicap
(erroneously imputed in 1987) and five overlapping schooling-degree
dummies are dropped; income is rescaled to 1000 DM. Missing values of
`hsat`, `educ` and the handicap indicator were imputed by the data
authors with gender-specific sample means (so a few `educ` values are
not whole or half years). Note that the waves are not equally spaced
(gaps after 1988 and 1991): lags and `dynamic = TRUE` fits only use
consecutive years, so dynamic examples restrict the data to 1984 to
1988.

## Teaching use

- Is income or add-on insurance associated with more doctor visits once
  stable individual differences are removed? Count models:
  `panglm(docvis ~ ..., family = "poisson", model = "within")`
  (conditional fixed-effects Poisson) against `model = "random"` or
  `family = "negbin"`, compared with the classical
  [`panhausman()`](https://elbadisyimad.com/panglm/reference/panhausman.md)
  test or a correlated random effects fit (`cre = TRUE`, Mundlak).

- Is a doctor visit this year partly caused by a visit last year (state
  dependence) or only by persistent health? Binary model with Wooldridge
  initial conditions on the consecutive waves 1984 to 1988:
  `panglm(doctor ~ ..., family = "probit", model = "random", dynamic = TRUE)`;
  [`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md)
  for the effect on the probability.

- Which life events change health satisfaction within a person? Ordered
  outcome: `panglm(hsat ~ ..., family = "ologit", model = "within")`
  (blow-up-and-cluster fixed-effects ordered logit) against
  `model = "pooling"` or `model = "random"`.

- Do covariates act on the decision to see a doctor or on how often one
  goes? Fixed-effects hurdle count model
  `panhurdle(docvis ~ ..., family = "poisson")` (time-varying regressors
  only: time-invariant ones such as `female` are dropped with a
  warning).

- Are there latent "low user" and "high user" classes, as in Deb and
  Trivedi? Finite mixture
  `panfmm(docvis ~ ..., k = 2, family = "negbin")`, with partial pooling
  by person.

- Effects of time-invariant characteristics (gender, schooling) when
  time-varying regressors may be correlated with the individual effect:
  `paniv(..., model = "hausmantaylor")`, or `cre = TRUE` in
  [`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md).
  (`educ` varies slightly within a few persons because of the
  imputation, so
  [`paniv()`](https://elbadisyimad.com/panglm/reference/paniv.md) treats
  it as time-varying.)

## References

Riphahn, R. T., Wambach, A. and Million, A. (2003). Incentive effects in
the demand for health care: a bivariate panel count data estimation.
*Journal of Applied Econometrics*, 18(4), 387-405.
[doi:10.1002/jae.680](https://doi.org/10.1002/jae.680)

Greene, W. H. (2018). *Econometric Analysis*, 8th ed. Pearson (uses
these German health care data in several examples).

Greene, W. H. (2014). Discrete choice modeling, Part 3: panel data
binary choice models. Lecture notes, Stern School of Business, New York
University (pooled and random-effects probit of `doctor`, reproduced by
the build script).

Greene, W. H. (2008). Modeling ordered choices. Working paper, Stern
School of Business, New York University (Table 3.1, ordered logit of
collapsed health satisfaction in 1988).

Hilbe, J. M. (2011). *Negative Binomial Regression*, 2nd ed. Cambridge
University Press.

## Examples

``` r
data(docvisits)
str(docvisits)
#> 'data.frame':    27326 obs. of  21 variables:
#>  $ id         : int  1 1 1 2 2 2 2 3 3 3 ...
#>  $ year       : int  1984 1985 1986 1984 1985 1986 1988 1984 1986 1987 ...
#>  $ docvis     : int  1 0 0 0 1 2 1 0 0 10 ...
#>  $ hospvis    : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ doctor     : int  1 0 0 0 1 1 1 0 0 1 ...
#>  $ hospital   : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ hsat       : int  8 8 7 7 8 7 8 10 9 10 ...
#>  $ female     : int  0 0 0 1 1 1 1 1 1 1 ...
#>  $ age        : int  54 55 56 44 45 46 48 58 60 61 ...
#>  $ educ       : num  15 15 15 9 9 9 9 11 11 11 ...
#>  $ married    : int  1 1 1 1 1 1 1 0 0 0 ...
#>  $ hhkids     : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ hhninc     : num  3.05 4.51 3.5 3.05 3.18 ...
#>  $ working    : int  1 1 1 0 0 0 0 0 0 0 ...
#>  $ bluecollar : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ whitecollar: int  1 1 1 0 0 0 0 0 0 0 ...
#>  $ self       : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ civil      : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ handicap   : int  0 0 0 0 0 0 0 0 1 1 ...
#>  $ public     : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ addon      : int  0 0 0 0 0 0 0 0 0 0 ...
table(docvisits$year)
#> 
#> 1984 1985 1986 1987 1988 1991 1994 
#> 3874 3794 3792 3666 4483 4340 3377 
# sample means by gender
aggregate(cbind(docvis, hospvis, doctor, addon, public, hsat) ~ female,
          data = docvisits, FUN = mean)
#>   female   docvis   hospvis    doctor      addon    public     hsat
#> 1      0 2.625711 0.1277821 0.5595029 0.01755248 0.8610546 6.924595
#> 2      1 3.790797 0.1496599 0.7048842 0.02017886 0.9125583 6.634411

idx <- c("id", "year")

# Count outcome: fixed-effects (conditional) Poisson versus random
# effects, with a Hausman test
fe <- panglm(docvis ~ age + hhninc + hhkids + married + working + addon,
             data = docvisits, index = idx, model = "within",
             family = "poisson")
re <- panglm(docvis ~ age + hhninc + hhkids + married + working + addon,
             data = docvisits, index = idx, model = "random",
             family = "poisson")
summary(fe)
#> 
#> Call:
#> panglm(formula = docvis ~ age + hhninc + hhkids + married + working + 
#>     addon, data = docvisits, index = idx, model = "within", family = "poisson")
#> 
#> Family: poisson  Link: log  Model: within  Effect: individual 
#> n = 27326  groups = 7293  vcov: classical 
#> 
#> Coefficients:
#>          Estimate Std. Error z value Pr(>|z|)    
#> age      0.029550   0.001453  20.335  < 2e-16 ***
#> hhninc  -0.024972   0.004211  -5.930 3.03e-09 ***
#> hhkids   0.002555   0.015733   0.162   0.8710    
#> married -0.044662   0.021756  -2.053   0.0401 *  
#> working -0.082671   0.013492  -6.128 8.93e-10 ***
#> addon   -0.078030   0.035683  -2.187   0.0288 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -47689 
#> Dispersion parameter: 1 
#> Residual df: 20027  Iterations: 3
#> 
panhausman(fe, re)
#> 
#>  Hausman specification test (fixed vs. random effects)
#> 
#> chisq = 131.3, df = 6, p-value = < 2.2e-16
#> alternative hypothesis: one model is inconsistent
#> 

# Ordinal outcome: fixed-effects ordered logit (blow-up and cluster)
# of health satisfaction
fo <- panglm(hsat ~ age + hhninc + hhkids + married + working + handicap,
             data = docvisits, index = idx, model = "within",
             family = "ologit")
summary(fo)
#> 
#> Call:
#> panglm(formula = hsat ~ age + hhninc + hhkids + married + working + 
#>     handicap, data = docvisits, index = idx, model = "within", 
#>     family = "ologit")
#> 
#> Family: ordered  Link: logit  Model: within  Effect: individual 
#> n = 27326  groups = 7293  vcov: cluster 
#> 
#> Coefficients:
#>           Estimate Std. Error z value Pr(>|z|)    
#> age      -0.101061   0.005848 -17.281  < 2e-16 ***
#> hhninc    0.043260   0.014430   2.998  0.00272 ** 
#> hhkids   -0.108355   0.057332  -1.890  0.05877 .  
#> married   0.078653   0.081018   0.971  0.33164    
#> working  -0.033533   0.057502  -0.583  0.55979    
#> handicap -0.302447   0.095196  -3.177  0.00149 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -31886 
#> Iterations: 4
#> 

# Binary outcome with state dependence, consecutive waves 1984-1988
early <- subset(docvisits, year <= 1988)
dp <- panglm(doctor ~ age + hhninc + hhkids + married + working,
             data = early, index = idx, model = "pooling",
             family = "probit", dynamic = TRUE, vcov = "cluster")
coef(dp)[c("lag_doctor", "initial_doctor")]
#>     lag_doctor initial_doctor 
#>      0.6453130      0.5457217 

# \donttest{
# Dynamic random-effects probit (Wooldridge initial conditions)
dr <- panglm(doctor ~ age + hhninc + hhkids + married + working,
             data = early, index = idx, model = "random",
             family = "probit", dynamic = TRUE, R = 8)
summary(dr)
#> 
#> Call:
#> panglm(formula = doctor ~ age + hhninc + hhkids + married + working, 
#>     data = early, index = idx, model = "random", family = "probit", 
#>     dynamic = TRUE, R = 8)
#> 
#> Family: binomial  Link: probit  Model: random  Effect: individual 
#> n = 11862  groups = 4578  vcov: classical 
#> 
#> Coefficients:
#>                Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)    -0.46859    0.10435  -4.491 7.10e-06 ***
#> age             0.03564    0.01465   2.433   0.0150 *  
#> hhninc          0.00838    0.02040   0.411   0.6812    
#> hhkids          0.01990    0.08823   0.226   0.8215    
#> married         0.08094    0.12941   0.625   0.5317    
#> working        -0.05942    0.08069  -0.736   0.4615    
#> lag_doctor      0.28250    0.04546   6.214 5.16e-10 ***
#> initial_doctor  0.91035    0.05176  17.588  < 2e-16 ***
#> mean_age       -0.02611    0.01474  -1.772   0.0765 .  
#> mean_hhninc    -0.04014    0.02407  -1.668   0.0954 .  
#> mean_hhkids    -0.14491    0.09799  -1.479   0.1392    
#> mean_married    0.04751    0.13791   0.345   0.7305    
#> mean_working   -0.18651    0.09119  -2.045   0.0408 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -6666 
#> Dispersion parameter: 1 
#> Random-intercept SD: 0.6453 
#> Mundlak test (individual means jointly zero): chi2(5) = 12.73, p = 0.026
#> Residual df: 11848  Iterations: 8
#> 

# Latent "low user" / "high user" classes (Deb-Trivedi style),
# women 1984-1988
women <- subset(docvisits, female == 1 & year <= 1988)
lc <- panfmm(docvis ~ age + educ + hhninc + hhkids + working + addon,
             data = women, index = idx, k = 2, family = "negbin",
             nstart = 3, seed = 1)
summary(lc)
#> 
#> Call:
#> panfmm(formula = docvis ~ age + educ + hhninc + hhkids + working + 
#>     addon, data = women, index = idx, k = 2, family = "negbin", 
#>     nstart = 3, seed = 1)
#> 
#> Model: 2-class latent class NB2 model 
#> n = 9422  groups = 3058  vcov: classical 
#> 
#> Class 1 (prior 0.616, 2005 individuals by modal posterior, mean outcome 1.589):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.046588   0.233423   0.200    0.842    
#> age          0.012500   0.002849   4.388 1.14e-05 ***
#> educ         0.002881   0.013182   0.219    0.827    
#> hhninc      -0.024563   0.015359  -1.599    0.110    
#> hhkids      -0.281618   0.063043  -4.467 7.93e-06 ***
#> working     -0.041179   0.054439  -0.756    0.449    
#> addon        0.247310   0.189783   1.303    0.193    
#> log(theta)  -0.333372   0.043434  -7.675 1.65e-14 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Class 2 (prior 0.384, 1053 individuals by modal posterior, mean outcome 7.162):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.915087   0.176266  10.865  < 2e-16 ***
#> age          0.010880   0.002151   5.057 4.25e-07 ***
#> educ        -0.016707   0.010656  -1.568 0.116917    
#> hhninc      -0.060582   0.012302  -4.924 8.46e-07 ***
#> hhkids      -0.179046   0.046069  -3.886 0.000102 ***
#> working     -0.034213   0.041000  -0.834 0.404016    
#> addon        0.018952   0.134168   0.141 0.887668    
#> log(theta)   0.205423   0.035269   5.824 5.73e-09 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> theta (NB2, variance mu + mu^2 / theta):
#> class1 class2 
#> 0.7165 1.2280 
#> 
#> Class probabilities (multinomial logit, class 1 = reference):
#>                    Estimate Std. Error z value Pr(>|z|)    
#> prior2:(Intercept) -0.47183    0.09796  -4.816 1.46e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> logLik = -21814.97  df = 17  AIC = 43663.94  BIC = 43785.5 
#> 

# Hausman-Taylor: effect of the time-invariant female indicator when
# income may be correlated with the individual effect
ht <- paniv(log(1 + docvis) ~ age + hhninc + working + female |
              age + working + female,
            data = docvisits, index = idx, model = "hausmantaylor")
summary(ht)
#> 
#> Call:
#> paniv(formula = log(1 + docvis) ~ age + hhninc + working + female | 
#>     age + working + female, data = docvisits, index = idx, model = "hausmantaylor")
#> 
#> Estimator: Hausman-Taylor 
#> n = 27326  groups = 7293  vcov: classical 
#> Exogenous time-varying:   age, working 
#> Endogenous time-varying:  hhninc 
#> Exogenous time-invariant: female 
#> Endogenous time-invariant: (none) 
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.2619675  0.0369980   7.081 1.44e-12 ***
#> age          0.0147995  0.0007294  20.291  < 2e-16 ***
#> hhninc       0.0009869  0.0041519   0.238    0.812    
#> working     -0.0875829  0.0147979  -5.919 3.25e-09 ***
#> female       0.2276753  0.0193597  11.760  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Variance components: idiosyncratic = 0.4976, individual = 0.4844
#> theta: 0.2882 to 0.6423 
#> 
#> Sargan overidentification test: chi2(2) = 61.51, p = 4.396e-14
#> Residual df: 27321 
#> 

# Replicates Greene's random-effects probit (DOCTOR on AGE, EDUC,
# HHNINC; log L = -16290.72); slow (about 1.5 minutes) on 27,326 rows:
# panglm(doctor ~ age + educ + I(hhninc / 10), data = docvisits,
#        index = idx, model = "random", family = "probit")
# }
```
