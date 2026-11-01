# Health care expenditure and use of US adults, MEPS panel 23 (2018-2019)

A real public two-year person panel from the US Medical Expenditure
Panel Survey (MEPS) Household Component, longitudinal file HC-217 (panel
23). The unit is a person, the time index is the calendar year (2018,
2019): each of 3,000 adults contributes exactly two rows. Expenditure
outcomes have the classic health-economics shape (a mass of zeros plus a
long right tail), utilisation outcomes are overdispersed counts with
many zeros, and insurance, income and chronic conditions vary over time
within person, which makes the data a natural testbed for two-part,
hurdle, tobit, gamma, quantile, selection and correlated random effects
(Mundlak) panel models.

## Usage

``` r
meps
```

## Format

A data frame with 6,000 rows (3,000 persons x 2 years) and 20 columns.
The panel index pair is `c("id", "year")`.

- id:

  integer person identifier (1-3000), the panel individual index; a new
  sequential id, not the MEPS `DUPERSID`

- year:

  integer calendar year (2018 or 2019), the panel time index (MEPS
  suffix Y1, Y2)

- totexp:

  numeric total health care expenditure in the year, US dollars, all
  payers (`TOTEXPY1`, `TOTEXPY2`); 14 percent zeros

- oopexp:

  numeric out-of-pocket (self or family) expenditure, US dollars
  (`TOTSLFY1`, `TOTSLFY2`); never larger than `totexp`

- anyexp:

  integer indicator (0/1) of any health care expenditure in the year,
  `as.integer(totexp > 0)`

- officevis:

  integer number of office-based provider visits (`OBTOTVY1`,
  `OBTOTVY2`)

- ervis:

  integer number of emergency room visits (`ERTOTY1`, `ERTOTY2`)

- ipstays:

  integer number of hospital inpatient discharges (`IPDISY1`, `IPDISY2`)

- insurance:

  factor health insurance coverage in the year (`INSCOVY1`, `INSCOVY2`):
  `private` (any private coverage), `public` (public coverage only, e.g.
  Medicare, Medicaid), `uninsured` (uninsured all year); reference level
  `private`

- age:

  integer age in years at the end of the year (or last age in scope),
  top-coded at 85 (`AGEY1X`, `AGEY2X`). Age rises by one year for
  (almost) everyone, so within person it is collinear with the year: use
  it in pooled and random effects models but not as a time-varying
  regressor in within or Mundlak (`cre = TRUE`) fits

- sex:

  factor `male`, `female` (`SEX`), time-invariant

- race:

  factor race and ethnicity (`RACETHX`): `hispanic`, `white`
  (non-Hispanic white only), `black` (non-Hispanic black only), `asian`
  (non-Hispanic Asian only), `other` (non-Hispanic other or multiple
  race); time-invariant

- region:

  factor Census region at the end of the year (`REGIONY1`, `REGIONY2`):
  `northeast`, `midwest`, `south`, `west`

- povcat:

  factor family income as a percent of the federal poverty line, in
  categories (`POVCATY1`, `POVCATY2`): `poor` (below 100 percent,
  including negative income), `nearpoor` (100 to below 125), `low` (125
  to below 200), `middle` (200 to below 400), `high` (400 or more)

- povlev:

  numeric family income as a percent of the poverty line, one decimal
  (`POVLEVY1`, `POVLEVY2`); can be negative

- health:

  ordered factor self-rated (perceived) health status, `excellent` \<
  `verygood` \< `good` \< `fair` \< `poor`, from the interview round
  closing the year (2018: round 3 `RTHLTH3`; 2019: round 5 `RTHLTH5`),
  as in the MEPS annual files; 4 values `NA` (refused, don't know or not
  ascertained)

- nchronic:

  integer number (0-8) of eight MEPS priority conditions ever diagnosed
  as reported in that year: high blood pressure, coronary heart disease,
  stroke, emphysema, high cholesterol, diabetes, arthritis, asthma
  (`HIBPDXYt`, `CHDDXYt`, `STRKDXYt`, `EMPHDXYt`, `CHOLDXYt`,
  `DIABDXYt_M18`, `ARTHDXYt`, `ASTHDXYt`); refused, don't know and
  inapplicable answers count as not diagnosed. The longitudinal file
  carries the same "ever diagnosed" answers in both years, so `nchronic`
  is constant within person: a baseline morbidity measure, dropped by
  within (fixed effects) estimators

- longwt:

  numeric MEPS panel 23 longitudinal person weight (`LONGWT`), constant
  within person

- varstr:

  integer variance estimation stratum (`VARSTR`)

- varpsu:

  integer variance estimation primary sampling unit within stratum
  (`VARPSU`)

## Source

Agency for Healthcare Research and Quality (AHRQ), Medical Expenditure
Panel Survey, Household Component, public use file HC-217: MEPS Panel 23
Longitudinal Data File (released December 2021). Rockville, MD: AHRQ.
<https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-217>.
Data: <https://meps.ahrq.gov/mepsweb/data_files/pufs/h217/h217dta.zip>.

MEPS is produced by AHRQ, an agency of the US federal government; the
public use files carry no copyright notice. AHRQ "requests that users
cite AHRQ and the Medical Expenditure Panel Survey as the data source in
any publications or research based upon these data." The files are
released under the following Data Use Agreement (section A of the HC-217
documentation,
<https://meps.ahrq.gov/mepsweb/data_stats/download_data/pufs/h217/h217doc.shtml>),
which binds every user of this extract as well:

"Individual identifiers have been removed from the micro-data contained
in these files. Nevertheless, under sections 308 (d) and 903 (c) of the
Public Health Service Act (42 U.S.C. 242m and 42 U.S.C. 299 a-1), data
collected by the Agency for Healthcare Research and Quality (AHRQ)
and/or the National Center for Health Statistics (NCHS) may not be used
for any purpose other than for the purpose for which they were supplied;
any effort to determine the identity of any reported cases is prohibited
by law. Therefore in accordance with the above referenced Federal
Statute, it is understood that: (1) No one is to use the data in this
data set in any way except for statistical reporting and analysis; and
(2) If the identity of any person or establishment should be discovered
inadvertently, then (a) no use will be made of this knowledge, (b) the
Director Office of Management AHRQ will be advised of this incident, (c)
the information that would identify any individual or establishment will
be safeguarded or destroyed, as requested by AHRQ, and (d) no one else
will be informed of the discovered identity; and (3) No one will attempt
to link this data set with individually identifiable records from any
data sets other than the Medical Expenditure Panel Survey or the
National Health Interview Survey. By using these data you signify your
agreement to comply with the above stated statutorily based requirements
with the knowledge that deliberately making a false statement in any
matter within the jurisdiction of any department or agency of the
Federal Government violates Title 18 part 1 Chapter 47 Section 1001 and
is punishable by a fine of up to \$10,000 or up to 5 years in prison."

Validation (in `data-raw/meps.R`): on the full HC-217 file the
unweighted codebook frequencies of `INSCOVY1`, `RACETHX`, `POVCATY2`,
`TOTEXPY1` and `OBTOTVY2` are reproduced exactly, and `LONGWT` sums to
305.7 million over the 13,044 `ALL5RDS = 1` persons as in the
documentation. The `LONGWT`-weighted mean total expenditure per person
is 5,676 dollars in 2018 and 6,249 dollars in 2019, against 6,063 and
6,252 dollars from the annual files in AHRQ Statistical Briefs \#533 and
\#540; the share with any expense is 87.1 and 84.4 percent against about
87 and 86 percent. The 2018 gap (-6 percent) is expected: the
longitudinal weight covers only people in scope for the whole two years,
which excludes 2018 decedents, whose end-of-life spending is high.

## Details

The extract keeps persons with data for all five interview rounds
(`ALL5RDS = 1`, in scope for the whole of 2018 and 2019) who were aged
18 or older in 2018, and draws a stratified random sample of 3,000 of
the 10,084 eligible persons (strata: 2018 age 18-34, 35-49, 50-64, 65+;
proportional allocation; `set.seed(20261101)`). The original person
identifier `DUPERSID` is replaced by a new sequential `id` so that rows
cannot be matched back to the public file or linked to other sources.
The build script, with its validation against the codebook and against
official AHRQ estimates, is `data-raw/meps.R` in the package's GitHub
repository.

Because the sample is drawn with equal probability within age strata,
`longwt` remains approximately valid for weighted descriptive statistics
of the adult population (the 2018 and 2019 weighted means of `totexp` in
the sample are within 2 percent of those of all eligible adults in the
full file), but the design variables `varstr` and `varpsu` only describe
the original MEPS design: design-based standard errors computed from
this subsample are approximate. The panglm examples below are unweighted
model-based fits, as usual for teaching.

## Teaching use

All outcomes are observed for every person in both years, so the data
are a balanced short panel (T = 2): within and Mundlak (CRE) estimators
are identified only from people whose insurance or income changed
between 2018 and 2019 (223 insurance changes, 1,114 income-category
changes), and `age` and `nchronic` should enter only pooled or random
effects fits.

- **Does insurance status change health spending, and through which
  margin?** `totexp` is zero for 14 percent of person-years and highly
  skewed:
  [`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md)
  (probit + gamma, `model = "pooling"` or `"random"`) with
  [`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md)
  for the dollar effect of being uninsured; compare with
  `panglm(family = "tobit")` on `log1p(totexp)`.

- **Is the cross-sectional insurance gap causal or selection?**
  [`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md)
  `family = "gamma"` on positive spending, `model = "pooling"` or
  `"random"` versus `model = "within"` (identified from people who
  gained or lost coverage).

- **Who uses any care?** `anyexp` with `family = "probit"`,
  `model = "random"`, `cre = TRUE` (Mundlak test of the random effects
  assumption) and
  [`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md).

- **Utilisation counts.** `officevis` (overdispersed, 24 percent zeros)
  or `ervis`, `ipstays` (mostly zero) with `family = "negbin"` or
  `"poisson"`, `model = "random"` or `"within"`;
  [`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md)
  for the zero versus positive split with fixed effects.

- **Distributional effects.**
  [`panqr()`](https://elbadisyimad.com/panglm/reference/panqr.md)
  (`model = "canay"` or `"pooling"`) on `log1p(totexp)` or
  `log(oopexp + 1)` at several quantiles; `model = "mmqr"` needs T \>= 3
  and does not apply here.

- **Selection on spending.**
  [`panselect()`](https://elbadisyimad.com/panglm/reference/panselect.md)
  with `anyexp` as the selection equation and log spending as the
  outcome (the alternative to the two-part model); self-rated `health`
  as an ordered outcome with `family = "oprobit"`, `model = "random"`.

## References

Mitchell, E. M. (2021). Concentration of healthcare expenditures and
selected characteristics of persons with high expenses, U.S. civilian
noninstitutionalized population, 2018. Statistical Brief \#533. Agency
for Healthcare Research and Quality, Rockville, MD.
<https://meps.ahrq.gov/data_files/publications/st533/stat533.shtml>

Mitchell, E. M. (2022). Concentration of healthcare expenditures and
selected characteristics of persons with high expenses, U.S. civilian
noninstitutionalized population, 2019. Statistical Brief \#540. Agency
for Healthcare Research and Quality, Rockville, MD.
<https://meps.ahrq.gov/data_files/publications/st540/stat540.shtml>

Cohen, J. W., Cohen, S. B. and Banthin, J. S. (2009). The Medical
Expenditure Panel Survey: a national information resource to support
healthcare cost research and inform policy and practice. *Medical Care*,
47(7 Suppl 1), S44-S50.
[doi:10.1097/MLR.0b013e3181a23e3a](https://doi.org/10.1097/MLR.0b013e3181a23e3a)

Deb, P., Norton, E. C. and Manning, W. G. (2017). *Health Econometrics
Using Stata*. College Station, TX: Stata Press. (Two-part, generalized
linear and count models for MEPS expenditure and use data.)

Duan, N., Manning, W. G., Morris, C. N. and Newhouse, J. P. (1983). A
comparison of alternative models for the demand for medical care.
*Journal of Business and Economic Statistics*, 1, 115-126.
[doi:10.1080/07350015.1983.10509330](https://doi.org/10.1080/07350015.1983.10509330)

## Examples

``` r
data(meps)
str(meps)
#> 'data.frame':    6000 obs. of  20 variables:
#>  $ id       : int  1 1 2 2 3 3 4 4 5 5 ...
#>  $ year     : int  2018 2019 2018 2019 2018 2019 2018 2019 2018 2019 ...
#>  $ totexp   : num  0 54 945 2452 923 ...
#>  $ oopexp   : num  0 0 247 459 139 ...
#>  $ anyexp   : int  0 1 1 1 1 1 1 1 1 1 ...
#>  $ officevis: int  0 1 3 1 5 0 14 20 1 3 ...
#>  $ ervis    : int  0 0 0 0 5 1 0 0 0 0 ...
#>  $ ipstays  : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ insurance: Factor w/ 3 levels "private","public",..: 1 1 2 2 2 3 2 2 1 1 ...
#>  $ age      : int  41 42 82 83 45 46 70 71 45 46 ...
#>  $ sex      : Factor w/ 2 levels "male","female": 2 2 1 1 2 2 1 1 2 2 ...
#>  $ race     : Factor w/ 5 levels "hispanic","white",..: 3 3 2 2 1 1 2 2 2 2 ...
#>  $ region   : Factor w/ 4 levels "northeast","midwest",..: 2 2 1 1 3 3 4 4 4 4 ...
#>  $ povcat   : Factor w/ 5 levels "poor","nearpoor",..: 4 4 4 1 4 1 2 3 5 5 ...
#>  $ povlev   : num  347 227 276 0 302 ...
#>  $ health   : Ord.factor w/ 5 levels "excellent"<"verygood"<..: 2 1 4 3 3 3 1 3 3 2 ...
#>  $ nchronic : int  0 0 3 3 0 0 1 1 0 0 ...
#>  $ longwt   : num  55960 55960 12860 12860 8498 ...
#>  $ varstr   : int  2019 2019 2079 2079 2028 2028 2077 2077 2018 2018 ...
#>  $ varpsu   : int  1 1 1 1 1 1 1 1 2 2 ...
# survey-weighted mean spending per adult and share with any expense
sapply(split(meps, meps$year), function(d)
  c(mean = weighted.mean(d$totexp, d$longwt),
    anyexp = weighted.mean(d$anyexp, d$longwt)))
#>               2018         2019
#> mean   6566.098685 7490.9550582
#> anyexp    0.866528    0.8375669

# two-part model: P(any spending) and E(spending | spending > 0)
tp <- twopart(totexp ~ insurance + I(age / 10) + sex + nchronic + povcat,
              data = meps, index = c("id", "year"), positive = "gamma")
summary(tp)
#> 
#> Call:
#> twopart(formula = totexp ~ insurance + I(age/10) + sex + nchronic + 
#>     povcat, data = meps, index = c("id", "year"), positive = "gamma")
#> 
#> Estimator: two-part model (pooling)
#> N = 6000  groups = 3000  positive = 5151  vcov: classical 
#> 
#> Part 1: P(y > 0), binomial probit 
#>                     Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)        -0.042247   0.095502  -0.442  0.65822    
#> insurancepublic    -0.020466   0.065467  -0.313  0.75457    
#> insuranceuninsured -1.128560   0.069907 -16.144  < 2e-16 ***
#> I(age/10)           0.117951   0.016309   7.232 4.76e-13 ***
#> sexfemale           0.407118   0.047051   8.653  < 2e-16 ***
#> nchronic            0.414266   0.030872  13.419  < 2e-16 ***
#> povcatnearpoor     -0.001501   0.126960  -0.012  0.99057    
#> povcatlow           0.135083   0.085579   1.578  0.11446    
#> povcatmiddle        0.206181   0.078269   2.634  0.00843 ** 
#> povcathigh          0.344645   0.080947   4.258 2.07e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Part 2: E[y | y > 0], gamma GLM with log link 
#>                    Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)         8.26642    0.15425  53.591  < 2e-16 ***
#> insurancepublic    -0.06043    0.09105  -0.664 0.506874    
#> insuranceuninsured -1.13808    0.18171  -6.263 3.77e-10 ***
#> I(age/10)           0.08908    0.02398   3.714 0.000204 ***
#> sexfemale           0.06290    0.07200   0.874 0.382289    
#> nchronic            0.25455    0.02907   8.755  < 2e-16 ***
#> povcatnearpoor     -0.05778    0.20985  -0.275 0.783041    
#> povcatlow          -0.18631    0.13909  -1.339 0.180411    
#> povcatmiddle       -0.21901    0.12604  -1.738 0.082280 .  
#> povcathigh         -0.08807    0.12818  -0.687 0.492029    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Log-likelihood: -52455 
#> 
partialeffects(tp, variables = "insurance")   # dollars per person-year
#> 
#> Average partial effects (total mean)
#> Delta-method standard errors, vcov: classical 
#> 
#>       term            contrast estimate std.error statistic   p.value conf.low
#>  insurance    public - private   -503.2     728.1    -0.691 4.896e-01    -1930
#>  insurance uninsured - private  -6194.0     563.4   -10.990 4.086e-28    -7299
#>  conf.high
#>        924
#>      -5090
#> 

# positive spending: pooled versus within (coverage switchers)
pos <- subset(meps, totexp > 0)
fe <- panglm(totexp ~ insurance + povcat, data = pos,
             index = c("id", "year"), family = "gamma", model = "within")
summary(fe)
#> 
#> Call:
#> panglm(formula = totexp ~ insurance + povcat, data = pos, index = c("id", 
#>     "year"), model = "within", family = "gamma")
#> 
#> Family: gamma  Link: log  Model: within  Effect: individual 
#> N = 5151  groups = 2757  vcov: classical 
#> 
#> Coefficients:
#>                    Estimate Std. Error z value Pr(>|z|)  
#> insurancepublic    -0.11058    0.09615  -1.150   0.2501  
#> insuranceuninsured -0.37052    0.15898  -2.331   0.0198 *
#> povcatnearpoor     -0.02014    0.09962  -0.202   0.8398  
#> povcatlow          -0.08441    0.07736  -1.091   0.2752  
#> povcatmiddle       -0.10368    0.07751  -1.338   0.1810  
#> povcathigh          0.02023    0.08660   0.234   0.8152  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -45869 
#> Dispersion parameter: 0.5693 
#> Residual df: 2388  Iterations: 15
#> 

# office-based visits: random effects negative binomial
nb <- panglm(officevis ~ insurance + I(age / 10) + sex + nchronic,
             data = meps, index = c("id", "year"), family = "negbin",
             model = "random")
summary(nb)
#> 
#> Call:
#> panglm(formula = officevis ~ insurance + I(age/10) + sex + nchronic, 
#>     data = meps, index = c("id", "year"), model = "random", family = "negbin")
#> 
#> Family: negbin  Link: log  Model: random  Effect: individual 
#> N = 6000  groups = 3000  vcov: classical 
#> 
#> Coefficients:
#>                    Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)        -0.60205    0.07190  -8.374   <2e-16 ***
#> insurancepublic    -0.10219    0.03979  -2.568   0.0102 *  
#> insuranceuninsured -1.31614    0.09217 -14.280   <2e-16 ***
#> I(age/10)           0.17873    0.01279  13.976   <2e-16 ***
#> sexfemale           0.36359    0.03888   9.353   <2e-16 ***
#> nchronic            0.24152    0.01540  15.683   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -16138 
#> Beta mixing parameters (Hausman, Hall and Griliches): a = 1.787  b = 2.793 
#> Residual df: 5992  Iterations: 5
#> 

# \donttest{
# any use: correlated random effects probit with Mundlak test
cre <- panglm(anyexp ~ insurance + povcat + sex + nchronic, data = meps,
              index = c("id", "year"), family = "probit",
              model = "random", cre = TRUE)
summary(cre)
#> 
#> Call:
#> panglm(formula = anyexp ~ insurance + povcat + sex + nchronic, 
#>     data = meps, index = c("id", "year"), model = "random", family = "probit", 
#>     cre = TRUE)
#> 
#> Family: binomial  Link: probit  Model: random  Effect: individual 
#> N = 6000  groups = 3000  vcov: classical 
#> 
#> Coefficients:
#>                         Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)              0.47111    0.16868   2.793  0.00522 ** 
#> insurancepublic          0.12784    0.26614   0.480  0.63097    
#> insuranceuninsured      -0.82562    0.29465  -2.802  0.00508 ** 
#> povcatnearpoor          -0.29140    0.26586  -1.096  0.27305    
#> povcatlow                0.09130    0.20985   0.435  0.66350    
#> povcatmiddle             0.14648    0.20923   0.700  0.48387    
#> povcathigh               0.07002    0.23738   0.295  0.76800    
#> sexfemale                0.68685    0.09246   7.429  1.1e-13 ***
#> nchronic                 0.81274    0.06029  13.479  < 2e-16 ***
#> mean_insurancepublic    -0.01803    0.29486  -0.061  0.95123    
#> mean_insuranceuninsured -1.04420    0.33235  -3.142  0.00168 ** 
#> mean_povcatnearpoor      0.56404    0.42942   1.314  0.18901    
#> mean_povcatlow           0.19728    0.28353   0.696  0.48656    
#> mean_povcatmiddle        0.30485    0.26921   1.132  0.25747    
#> mean_povcathigh          0.76924    0.29410   2.616  0.00891 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1756 
#> Dispersion parameter: 1 
#> Random-intercept SD: 1.275 
#> Mundlak test (individual means jointly zero): chi2(6) = 21.37, p = 0.001572
#> Residual df: 5984  Iterations: 98
#> 

# censored log spending: tobit
summary(panglm(log1p(totexp) ~ insurance + I(age / 10) + sex + nchronic,
               data = meps, index = c("id", "year"), family = "tobit",
               model = "pooling"))
#> 
#> Call:
#> panglm(formula = log1p(totexp) ~ insurance + I(age/10) + sex + 
#>     nchronic, data = meps, index = c("id", "year"), model = "pooling", 
#>     family = "tobit")
#> 
#> Family: tobit  Link: identity  Model: pooling  Effect: individual 
#> N = 6000  groups = 3000  vcov: classical 
#> 
#> Coefficients:
#>                    Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)         3.92516    0.12973  30.256  < 2e-16 ***
#> insurancepublic    -0.37754    0.09132  -4.134 3.56e-05 ***
#> insuranceuninsured -4.27667    0.15789 -27.087  < 2e-16 ***
#> I(age/10)           0.34521    0.02629  13.133  < 2e-16 ***
#> sexfemale           0.80404    0.07960  10.101  < 2e-16 ***
#> nchronic            0.67195    0.03296  20.387  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -14031 
#> Dispersion parameter: 9.096 
#> Residual df: 5993  Iterations: 87
#> 

# quantiles of log spending with individual effects (the many tied
# values of log1p(totexp) make quantreg warn that the solution may be
# nonunique)
qr <- suppressWarnings(panqr(log1p(totexp) ~ insurance + povlev,
  data = meps, index = c("id", "year"), tau = c(0.25, 0.5, 0.75),
  model = "canay", boot = 50, seed = 1))
summary(qr)
#> 
#> Call:
#> panqr(formula = log1p(totexp) ~ insurance + povlev, data = meps, 
#>     index = c("id", "year"), tau = c(0.25, 0.5, 0.75), model = "canay", 
#>     boot = 50, seed = 1)
#> 
#> Estimator: Canay (2011) two-step fixed-effects quantile regression 
#> N = 6000  groups = 3000  vcov: bootstrap (cluster, boot = 50) 
#> 
#> Quantile tau = 0.25:
#>                      Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)         6.3883695  0.1057939  60.385   <2e-16 ***
#> insurancepublic    -0.0229213  0.2011713  -0.114   0.9093    
#> insuranceuninsured -1.0552424  0.4272558  -2.470   0.0135 *  
#> povlev             -0.0001311  0.0001501  -0.873   0.3827    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Quantile tau = 0.5:
#>                      Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)         6.9475052  0.1087264  63.899   <2e-16 ***
#> insurancepublic    -0.0949217  0.2007992  -0.473   0.6364    
#> insuranceuninsured -1.3055085  0.3999694  -3.264   0.0011 ** 
#> povlev             -0.0002027  0.0001572  -1.289   0.1972    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Quantile tau = 0.75:
#>                      Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)         7.5060244  0.1130061  66.421  < 2e-16 ***
#> insurancepublic    -0.1976571  0.1983995  -0.996 0.319125    
#> insuranceuninsured -1.3097238  0.3417253  -3.833 0.000127 ***
#> povlev             -0.0002764  0.0001675  -1.650 0.098989 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 

# selection correction for log spending observed only if positive
meps$logexp <- ifelse(meps$totexp > 0, log(meps$totexp), NA)
summary(panselect(anyexp ~ insurance + povcat, logexp ~ insurance,
                  data = meps, index = c("id", "year"),
                  model = "within", boot = 49, seed = 1))
#> 
#> Call:
#> panselect(selection = anyexp ~ insurance + povcat, outcome = logexp ~ 
#>     insurance, data = meps, index = c("id", "year"), model = "within", 
#>     boot = 49, seed = 1)
#> 
#> Wooldridge (1995) panel selection correction: fixed effects with Mundlak means, inverse Mills ratio by period 
#> Selection sample: N = 6000  groups = 3000 
#> Outcome sample:   N = 5151  groups = 2757 
#> Probits fitted for periods: 2018, 2019 
#> 
#> Outcome equation (vcov: panel bootstrap, 49 replicates):
#>                    Estimate Std. Error z value Pr(>|z|)  
#> insurancepublic     -0.1047     0.1253  -0.835   0.4035  
#> insuranceuninsured  -0.5953     0.3134  -1.900   0.0575 .
#> lambda:2018          0.7878     1.1659   0.676   0.4992  
#> lambda:2019          0.8570     1.0147   0.845   0.3983  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Test of selection bias (H0: all inverse Mills ratio coefficients are zero):
#>   cluster-robust Wald chi2(2) = 2.757, p = 0.252
#>   bootstrap Wald chi2(2) = 2.032, p = 0.3621
#> Residual df: 2390 
#> 

# self-rated health: random effects ordered probit
summary(panglm(health ~ nchronic + I(age / 10) + povcat, data = meps,
               index = c("id", "year"), family = "oprobit",
               model = "random"))
#> 
#> Call:
#> panglm(formula = health ~ nchronic + I(age/10) + povcat, data = meps, 
#>     index = c("id", "year"), model = "random", family = "oprobit")
#> 
#> Family: ordered  Link: probit  Model: random  Effect: individual 
#> N = 5996  groups = 3000  vcov: classical 
#> 
#> Coefficients:
#>                    Estimate Std. Error z value Pr(>|z|)    
#> nchronic            0.47591    0.02314  20.565  < 2e-16 ***
#> I(age/10)          -0.01577    0.01734  -0.909 0.363251    
#> povcatnearpoor     -0.20687    0.10846  -1.907 0.056487 .  
#> povcatlow          -0.27450    0.07652  -3.587 0.000334 ***
#> povcatmiddle       -0.45216    0.06936  -6.519 7.08e-11 ***
#> povcathigh         -0.77657    0.07086 -10.959  < 2e-16 ***
#> excellent|verygood -1.31156    0.09513 -13.786  < 2e-16 ***
#> verygood|good       0.23781    0.09392   2.532 0.011338 *  
#> good|fair           1.94846    0.10036  19.415  < 2e-16 ***
#> fair|poor           3.37760    0.11527  29.302  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -7344 
#> Random-intercept SD: 1.16 
#> Residual df: 5985  Iterations: 75
#> 
# }
```
