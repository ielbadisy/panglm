# Health insurance coverage and the ACA Medicaid expansion, US states 2009-2019

A real public state-year panel of health insurance coverage for the 50
US states and the District of Columbia, 2009 to 2019, built from the
American Community Survey (ACS) 1-year estimates of the US Census
Bureau, with the date on which each state implemented the Affordable
Care Act (ACA) Medicaid expansion to adults with incomes up to 138
percent of the federal poverty level. Adoption was staggered (2014,
2015, 2016, 2017, 2019, with ten states never adopting and several
adopting after 2019), which makes it the standard teaching example for
staggered difference-in-differences. The panel unit is the state and the
time unit is the calendar year; the panel is balanced (51 x 11).

## Usage

``` r
medicaid
```

## Format

A data frame with 561 rows (51 states x 11 years) and 20 columns. The
panel index pair is `c("state", "year")`.

- state:

  factor, state name (51 levels, District of Columbia included); the
  panel individual index

- year:

  integer calendar year of the ACS estimate (2009-2019); the panel time
  index

- abb:

  factor, two-letter postal abbreviation

- fips:

  integer state FIPS code

- region:

  factor, Census region (Northeast, Midwest, South, West)

- pop:

  integer civilian noninstitutionalized population (all ages), the
  universe of the health insurance tables (B27001)

- uninsured:

  numeric percent of the civilian noninstitutionalized population
  without health insurance coverage, all ages (B27001)

- adultpop:

  integer adults aged 18 to 64 (2009-2016) or 19 to 64 (2017-2019) for
  whom poverty status is determined (B27016); the ACS moved its age
  break from 18 to 19 in 2017, for all states at once

- uninsadult:

  numeric percent of `adultpop` without health insurance coverage

- medicaidadult:

  numeric percent of `adultpop` with Medicaid or other means-tested
  public coverage (alone or with other coverage)

- povpop:

  integer adults (same ages) with family income below 100 percent of the
  poverty threshold

- povbase:

  integer `povpop` of the state in 2013, the last year before the main
  expansion wave: a time-invariant weight, as needed by
  [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md)
  (weights constant within individuals)

- uninspov:

  numeric percent of `povpop` without health insurance coverage; the
  low-income outcome available in every year

- medicaidpov:

  numeric percent of `povpop` with Medicaid or other means-tested public
  coverage

- lowpop:

  integer adults (same ages) with family income below 138 percent of the
  poverty threshold, the ACA eligibility cut; `NA` in 2009-2010, when
  the ACS tables cut at 100 and 150 percent

- uninslow:

  numeric percent of `lowpop` without health insurance coverage; `NA` in
  2009-2010

- adopt:

  integer expansion year: the first calendar year coded as treated, `NA`
  for the 10 states that had not adopted by October 2026. Convention of
  Miller, Johnson and Wherry (2021): the year of implementation if the
  expansion took effect before July 1, otherwise the following year (so
  Michigan, April 2014, is 2014; New Hampshire, August 2014, is 2015;
  Louisiana, July 2016, is 2017). Values after 2019 (Idaho and Utah
  2020, Nebraska 2021, Missouri and Oklahoma 2022, North Carolina and
  South Dakota 2024) are kept: in this panel those states are never
  treated, and
  [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md) treats
  cohorts first treated after the last period as never treated

- expdate:

  Date the expansion took effect (KFF); `NA` if not adopted

- expanded:

  integer 0/1 indicator that `year >= adopt`

- prior:

  factor coverage of low-income adults before 2014: `"full"` for DE, DC,
  MA, NY and VT, which already covered adults up to the poverty line or
  more (Miller, Johnson and Wherry 2021, who drop them); `"early"` for
  CA, CT, MN, NJ and WA, which made limited early expansions in
  2010-2012 under the ACA option or waivers (Sommers et al. 2013);
  `"none"` otherwise. A standard robustness check drops
  `prior != "none"`. Wisconsin is coded as a non-adopter but has covered
  adults up to 100 percent of poverty since 2014 through a waiver

## Source

US Census Bureau, American Community Survey 1-year estimates, 2009-2019,
detailed tables B27001 (health insurance coverage status by sex by age)
and B27016 (health insurance coverage status and type by ratio of income
to poverty level by age), Summary File,
<https://www2.census.gov/programs-surveys/acs/summary_file/>. Works of
the US federal government are in the public domain (17 U.S.C. 105); the
Census Bureau states "You do not need permission to use American
Community Survey data"
(<https://www.census.gov/help/topics/faq.how-do-i-cite-or-source-acs-data.html>).
Attribution: "Source: U.S. Census Bureau, American Community Survey
1-year estimates." This product is not endorsed or certified by the
Census Bureau. Expansion dates: KFF, Status of State Action on the
Medicaid Expansion Decision (State Health Facts), status as of October
2026,
<https://www.kff.org/affordable-care-act/state-indicator/state-activity-around-expanding-medicaid-under-the-affordable-care-act/>.

## Details

The ACS starts asking about health insurance in 2008, but the 2008
Summary File has no table by ratio of income to poverty nor a Medicaid
table, so the panel starts in 2009. It stops in 2019 because the Census
Bureau did not release standard 2020 ACS 1-year estimates (only
experimental ones), after the pandemic disrupted data collection. The
rates are computed from the published counts of tables B27001 and B27016
(state level), so they reproduce the published percentages; see
`data-raw/medicaid.R` in the package's GitHub repository for the build
and its validation checks.

## Teaching use

- Did the expansion reduce uninsurance among low-income adults?
  Staggered DiD with
  [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md)
  (`treat = "adopt"`, never-treated or not-yet-treated controls),
  aggregated by cohort and event time with
  [`summary.panglm_attgt()`](https://elbadisyimad.com/panglm/reference/summary.panglm_attgt.md).
  Weights must be constant within a state: use the base-year population
  `povbase` (`povpop` in 2013).

- How large is the bias of the single-coefficient two-way fixed effects
  regression? Compare
  [`did()`](https://elbadisyimad.com/panglm/reference/did.md) (or
  [`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md) with
  `model = "within"`, `effect = "twoways"` and the `expanded` dummy)
  with [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md)
  and [`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md):
  the effect grows with time since expansion, so early cohorts serve as
  controls for later ones.

- Are pre-trends parallel, and how does the effect build up?
  [`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md)
  with leads and lags, and `summary(type = "dynamic")` of
  [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md).
  Louisiana (July 2016, coded 2017) is partly treated in 2016: compare
  `anticipation = 1`.

- Inference with 51 clusters: the TWFE coefficient from
  [`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md) with
  `model = "within"` and `effect = "twoways"` has a classical standard
  error less than half the cluster-robust one (`vcov = "cluster"`), the
  default of
  [`did()`](https://elbadisyimad.com/panglm/reference/did.md).

- Did Medicaid enrolment (`medicaidpov`, `medicaidadult`) rise, and is
  the fall in uninsurance among all ages (`uninsured`) smaller, as
  expected from a policy aimed at low-income adults? Same estimators
  with other outcomes (a "first stage" and a dilution check).

- Robustness: drop the states with prior coverage (`prior != "none"`),
  restrict to the ACA income cut (`uninslow`, 2011-2019), or use
  population weights versus none. Weighted fixed-effects fits are not
  available in
  [`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md); a
  weighted TWFE is `model = "pooling"` with state and year dummies.

## References

Miller, S., Johnson, N. and Wherry, L. R. (2021). Medicaid and
mortality: new evidence from linked survey and administrative data.
*Quarterly Journal of Economics*, 136(3), 1783-1829.
[doi:10.1093/qje/qjab004](https://doi.org/10.1093/qje/qjab004)

Courtemanche, C., Marton, J., Ukert, B., Yelowitz, A. and Zapata, D.
(2017). Early impacts of the Affordable Care Act on health insurance
coverage in Medicaid expansion and non-expansion states. *Journal of
Policy Analysis and Management*, 36(1), 178-210.
[doi:10.1002/pam.21961](https://doi.org/10.1002/pam.21961)

Sommers, B. D., Arntson, E., Kenney, G. M. and Epstein, A. M. (2013).
Lessons from early Medicaid expansions under health reform: interviews
with Medicaid officials. *Medicare and Medicaid Research Review*, 3(4).
[doi:10.5600/mmrr.003.04.a02](https://doi.org/10.5600/mmrr.003.04.a02)

Smith, J. C. and Medalia, C. (2015). *Health Insurance Coverage in the
United States: 2014*. Current Population Reports P60-253, US Census
Bureau.

Callaway, B. and Sant'Anna, P. H. C. (2021). Difference-in-differences
with multiple time periods. *Journal of Econometrics*, 225(2), 200-230.
[doi:10.1016/j.jeconom.2020.12.001](https://doi.org/10.1016/j.jeconom.2020.12.001)

Goodman-Bacon, A. (2021). Difference-in-differences with variation in
treatment timing. *Journal of Econometrics*, 225(2), 254-277.
[doi:10.1016/j.jeconom.2021.03.014](https://doi.org/10.1016/j.jeconom.2021.03.014)

## Examples

``` r
data(medicaid)
str(medicaid)
#> 'data.frame':    561 obs. of  20 variables:
#>  $ state        : Factor w/ 51 levels "Alabama","Alaska",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ year         : int  2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 ...
#>  $ abb          : Factor w/ 51 levels "AK","AL","AR",..: 2 2 2 2 2 2 2 2 2 2 ...
#>  $ fips         : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ region       : Factor w/ 4 levels "Northeast","Midwest",..: 3 3 3 3 3 3 3 3 3 3 ...
#>  $ pop          : int  4616028 4702769 4720067 4741698 4755381 4767383 4781070 4782529 4793512 4810094 ...
#>  $ uninsured    : num  13.8 14.6 14.3 13.3 13.6 ...
#>  $ adultpop     : int  2832614 2895534 2907812 2902592 2910291 2898691 2892886 2882442 2814970 2807923 ...
#>  $ uninsadult   : num  19.9 21.2 21 20.1 20.3 ...
#>  $ medicaidadult: num  9.79 10.13 9.94 10.94 10.95 ...
#>  $ povpop       : int  459298 508921 517524 510882 511039 526239 511475 467307 445307 447100 ...
#>  $ povbase      : int  511039 511039 511039 511039 511039 511039 511039 511039 511039 511039 ...
#>  $ uninspov     : num  41.6 44.1 42.8 42.8 41.4 ...
#>  $ medicaidpov  : num  30.5 28.6 27.4 29.1 29.4 ...
#>  $ lowpop       : int  NA NA 749349 749124 741692 741218 720302 677523 643151 633544 ...
#>  $ uninslow     : num  NA NA 41.8 40.1 40.1 ...
#>  $ adopt        : int  NA NA NA NA NA NA NA NA NA NA ...
#>  $ expdate      : Date, format: NA NA ...
#>  $ expanded     : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ prior        : Factor w/ 3 levels "none","early",..: 1 1 1 1 1 1 1 1 1 1 ...
# adoption cohorts (first treated year; NA = not adopted by 2026)
table(unique(medicaid[, c("state", "adopt")])$adopt, useNA = "ifany")
#> 
#> 2014 2015 2016 2017 2019 2020 2021 2022 2024 <NA> 
#>   26    3    2    1    2    2    1    2    2   10 
ix <- c("state", "year")

# Callaway-Sant'Anna ATT(g, t) for the uninsured rate of adults below
# the poverty line, never-treated controls, states weighted equally
cs <- attgt(uninspov ~ 1, data = medicaid, index = ix, treat = "adopt",
            control = "never")
summary(cs)
#> 
#> Callaway and Sant'Anna group-time ATT
#> Aggregation: overall (simple weighted average)
#> Control group: never treated  Estimator: unconditional  Individuals: 51 
#> 
#> ATT: -6.54 (SE 1.295, 95% CI -9.078 to -4.002)
#> 
summary(cs, type = "dynamic")   # leads test pre-trends, lags the build-up
#> 
#> Callaway and Sant'Anna group-time ATT
#> Aggregation: event study (by time since treatment)
#> Control group: never treated  Estimator: unconditional  Individuals: 51 
#> 
#> ATT (mean of e >= 0): -6.562 (SE 1.317, 95% CI -9.143 to -3.98)
#> 
#>  event estimate     se     lower   upper
#>     -9  0.97529 0.8140  -0.62013  2.5707
#>     -8 -0.09853 0.6997  -1.46999  1.2729
#>     -7  1.16922 1.0103  -0.81093  3.1494
#>     -6 -1.12859 1.1047  -3.29372  1.0365
#>     -5  1.48015 0.7483   0.01347  2.9468
#>     -4 -0.63116 0.4610  -1.53479  0.2725
#>     -3 -0.29785 0.5167  -1.31050  0.7148
#>     -2  0.40279 0.3635  -0.30965  1.1152
#>     -1 -0.38889 0.5469  -1.46071  0.6829
#>      0 -5.15564 0.9010  -6.92165 -3.3896
#>      1 -6.39822 1.2670  -8.88146 -3.9150
#>      2 -7.08836 1.3521  -9.73834 -4.4384
#>      3 -7.59013 1.4220 -10.37726 -4.8030
#>      4 -6.52274 1.5487  -9.55810 -3.4874
#>      5 -6.61563 1.8257 -10.19386 -3.0374
#> 

# weighted by the pre-expansion (2013) low-income adult population,
# with not-yet-treated controls
cs2 <- attgt(uninspov ~ 1, data = medicaid, index = ix, treat = "adopt",
             control = "notyet", weights = "povbase")
summary(cs2)
#> 
#> Callaway and Sant'Anna group-time ATT
#> Aggregation: overall (simple weighted average)
#> Control group: not yet treated  Estimator: unconditional  Individuals: 51 
#> 
#> ATT: -7.018 (SE 1.428, 95% CI -9.815 to -4.22)
#> 

# extended TWFE and the single-coefficient TWFE for comparison (the
# latter can be biased under staggered adoption with dynamic effects)
ew <- etwfe(uninspov ~ 1, data = medicaid, index = ix, treat = "adopt")
summary(ew)
#> 
#> Extended TWFE (Wooldridge), linear model, control group: not yet treated
#> Observations: 561  Individuals: 51  vcov: cluster 
#> Fixed effects: cohort and period
#> 
#> ATT (response scale): -7.116 (SE 1.4, 95% CI -9.86 to -4.373, p = 3.705e-07)
#> 
tw <- did(uninspov ~ 1, data = medicaid, index = ix, treat = "adopt")
coef(tw)["treated"]
#>   treated 
#> -7.206276 

# event-study regression with a joint pre-trend test
eventstudy(uninspov ~ 1, data = medicaid, index = ix, treat = "adopt",
           leads = 4, lags = 5)
#> 
#> Event study (two-way fixed effects, gaussian)
#> n = 561  groups = 51  vcov: cluster 
#> Reference event time: -1 (end points binned) 
#> 
#>  event estimate     se    lower  upper
#>     -4   0.7173 0.7083  -0.6709  2.106
#>     -3   0.7217 0.5023  -0.2628  1.706
#>     -2   0.2804 0.4149  -0.5327  1.094
#>     -1   0.0000 0.0000   0.0000  0.000
#>      0  -5.5588 0.8460  -7.2170 -3.901
#>      1  -7.4610 1.1993  -9.8115 -5.110
#>      2  -7.6228 1.3074 -10.1853 -5.060
#>      3  -7.7790 1.5098 -10.7382 -4.820
#>      4  -6.5449 1.6404  -9.7601 -3.330
#>      5  -6.3989 2.0020 -10.3228 -2.475
#> 
#> Pre-trend test: Chi-squared = 2.334, df = 3, p-value = 0.506

# first stage: Medicaid coverage of adults below the poverty line
summary(attgt(medicaidpov ~ 1, data = medicaid, index = ix,
              treat = "adopt", control = "never"))
#> 
#> Callaway and Sant'Anna group-time ATT
#> Aggregation: overall (simple weighted average)
#> Control group: never treated  Estimator: unconditional  Individuals: 51 
#> 
#> ATT: 10.24 (SE 1.406, 95% CI 7.48 to 12.99)
#> 

# \donttest{
# the same TWFE coefficient with classical and cluster-robust SEs
fe <- panglm(uninspov ~ expanded, data = medicaid, index = ix,
             model = "within", effect = "twoways")
fecl <- update(fe, vcov = "cluster")
cbind(classical = sqrt(diag(vcov(fe))), cluster = sqrt(diag(vcov(fecl))))
#>          classical  cluster
#> expanded 0.5604576 1.202028

# weighted TWFE: a pooled fit with state and year dummies
tww <- panglm(uninspov ~ expanded + state + factor(year), data = medicaid,
              index = ix, model = "pooling", weights = "povbase",
              vcov = "cluster")
coef(summary(tww))["expanded", ]
#>      Estimate    Std. Error       z value      Pr(>|z|) 
#> -7.022593e+00  1.496188e+00 -4.693658e+00  2.683625e-06 

# robustness: drop the states with pre-2014 coverage of these adults
summary(attgt(uninspov ~ 1, data = subset(medicaid, prior == "none"),
              index = ix, treat = "adopt", control = "never"))
#> 
#> Callaway and Sant'Anna group-time ATT
#> Aggregation: overall (simple weighted average)
#> Control group: never treated  Estimator: unconditional  Individuals: 41 
#> 
#> ATT: -8.878 (SE 1.255, 95% CI -11.34 to -6.418)
#> 
# }
```
