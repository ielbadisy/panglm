# Synthetic COPD follow-up panel

A simulated longitudinal (panel) dataset of chronic obstructive
pulmonary disease (COPD) patients followed over several clinic visits,
used throughout this package's examples and vignette. It is entirely
synthetic (not real patient data), generated to have a realistic panel
structure (time-invariant baseline covariates plus an individual-level
frailty or random intercept), so it exercises pooled, fixed-effects and
random-effects fits for continuous, count and binary outcomes without
depending on any suggested package. The generating code is in
`data-raw/copd.R` in the package's GitHub repository.

## Usage

``` r
copd
```

## Format

A data frame with 300 rows (60 patients x 5 visits) and 10 columns:

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

## Source

Simulated; see `data-raw/copd.R` in the package's GitHub repository.

## Examples

``` r
data(copd)
str(copd)
#> 'data.frame':    300 obs. of  10 variables:
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
table(copd$visit)
#> 
#>  1  2  3  4  5 
#> 60 60 60 60 60 
```
