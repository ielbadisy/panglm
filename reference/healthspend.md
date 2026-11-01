# Health expenditure, income and population health: a country-year panel

A real public panel of 193 countries observed every year from 2000 to
2022 (23 years), assembled from the World Bank World Development
Indicators (WDI). It links current health expenditure (from the WHO
Global Health Expenditure Database, as redistributed in WDI), its
financing mix (out-of-pocket and domestic government shares), real
income, population ageing, life expectancy and under-5 mortality. The
panel is balanced in its rows (every country has all 23 years) but
unbalanced in its values: some country-years are missing (`NA`) for some
indicators, so fitted models use an unbalanced panel. Countries with no
health expenditure data in any year (mostly small territories) and all
regional and income aggregates are excluded. The build script, with the
download calls and the validation checks, is `data-raw/healthspend.R` in
the package's GitHub repository.

## Usage

``` r
healthspend
```

## Format

A data frame with 4,439 rows (193 countries x 23 years) and 15 columns.
The panel index pair is (`iso3c`, `year`).

- iso3c:

  factor, ISO 3166-1 alpha-3 country code (193 levels); the panel
  individual index

- year:

  integer calendar year, 2000-2022; the panel time index

- country:

  character, World Bank country name

- region:

  factor, World Bank region (7 levels: East Asia & Pacific; Europe &
  Central Asia; Latin America & Caribbean; Middle East, North Africa,
  Afghanistan & Pakistan; North America; South Asia; Sub-Saharan
  Africa); time invariant

- income:

  factor, World Bank income group with levels `low`, `lowermiddle`,
  `uppermiddle`, `high`; the classification current at download
  (2026-11-01), so it is time invariant here and does not track past
  reclassifications

- chepc:

  numeric, current health expenditure (CHE) per capita, PPP, current
  international dollars (WDI `SH.XPD.CHEX.PP.CD`), 2 decimals

- chereal:

  numeric, real CHE per capita in constant 2021 international dollars,
  derived as `chegdp / 100 * gdppc`; use this (not `chepc`, which is in
  current prices) together with `gdppc` for income elasticities

- chegdp:

  numeric, CHE as a percentage of GDP (`SH.XPD.CHEX.GD.ZS`)

- oop:

  numeric, out-of-pocket expenditure as a percentage of CHE
  (`SH.XPD.OOPC.CH.ZS`), 0-100; divide by 100 for a proportion

- gghe:

  numeric, domestic general government health expenditure as a
  percentage of CHE (`SH.XPD.GHED.CH.ZS`), 0-100

- gdppc:

  numeric, GDP per capita, PPP, constant 2021 international dollars
  (`NY.GDP.PCAP.PP.KD`), 1 decimal

- lifeexp:

  numeric, life expectancy at birth, total, years (`SP.DYN.LE00.IN`);
  period estimates include crisis years, e.g. 14.7 for the Central
  African Republic in 2009

- u5mr:

  numeric, under-5 mortality rate per 1,000 live births (`SH.DYN.MORT`)

- pop:

  integer, total population (`SP.POP.TOTL`)

- pop65:

  numeric, population aged 65 and above, percentage of total
  (`SP.POP.65UP.TO.ZS`)

Missing values: `chepc` 92, `chereal` 217, `chegdp`, `oop` and `gghe`
91, `gdppc` 174, `u5mr` 23; `lifeexp`, `pop` and `pop65` are complete.

## Source

World Bank, World Development Indicators, API v2
(<https://api.worldbank.org/v2/country/all/indicator/SH.XPD.CHEX.PP.CD?date=2000:2023&format=json&per_page=20000>
and the same call for the other indicator codes; country metadata from
<https://api.worldbank.org/v2/country?format=json&per_page=400>), WDI
release last updated 2026-10-08, downloaded 2026-11-01. Health
expenditure series: World Health Organization, Global Health Expenditure
Database (update of 12 December 2025), via WDI. Life expectancy and
population: UN World Population Prospects; under-5 mortality: UN
Inter-agency Group for Child Mortality Estimation; GDP: World Bank
International Comparison Program; all via WDI.

Licence: Creative Commons Attribution 4.0 (CC BY 4.0). The WDI catalog
entry
(<https://datacatalog.worldbank.org/search/dataset/0037712/World-Development-Indicators>)
states "This dataset is licensed under Creative Commons Attribution
4.0", and the metadata of each of the nine series used carries
`License_Type` "CC BY-4.0". Attribution: "World Bank, World Development
Indicators (CC BY 4.0)". Changes from the source: country subset,
renamed columns, rounding, and the derived column `chereal`.

## Teaching use

The panel is suited to the classic macro health-economics questions.
Take logs first, e.g.
`hs <- transform(healthspend, lche = log(chereal), lgdp = log(gdppc))`.

- Is health care a luxury good? Income elasticity of real health
  spending, `lche ~ lgdp`: pooled vs between vs within vs random effects
  ([`panglm()`](https://elbadisyimad.com/panglm/reference/panglm.md)
  with `model = "pooling"`, `"between"`, `"within"`, `"random"`),
  correlated random effects (`cre = TRUE`) and
  [`panhausman()`](https://elbadisyimad.com/panglm/reference/panhausman.md)
  (classical and Mundlak). Year effects (`effect = "twoways"`) remove
  common shocks and lower the within elasticity from about 1.15 to about
  0.84.

- How persistent is health spending, and what is the long-run income
  elasticity? Dynamic model with a lagged dependent variable,
  [`pangmm()`](https://elbadisyimad.com/panglm/reference/pangmm.md)
  (difference GMM `"d"` and system GMM `"ld"`, with `collapse = TRUE` to
  keep the instrument count below the number of countries); long-run
  elasticity = b / (1 - rho).

- What do time-invariant characteristics add once income is allowed to
  be correlated with the country effect? Hausman-Taylor with income
  endogenous and population ageing plus region exogenous:
  `paniv(..., model = "hausmantaylor")`.

- Does out-of-pocket financing fall as countries grow richer? The
  out-of-pocket share `oop / 100` is a proportion in (0, 1):
  `panglm(family = "beta")`, pooled with clustered standard errors or
  random effects with `cre = TRUE`.

- Does the income elasticity differ across countries or across the
  spending distribution? Random slopes (`random = ~ lgdp`, with `lgdp`
  centred) and fixed-effects quantile regression
  ([`panqr()`](https://elbadisyimad.com/panglm/reference/panqr.md),
  `"mmqr"` or `"canay"`).

- Does spending buy health? Life expectancy or under-5 mortality on
  `lche` and `lgdp` with two-way fixed effects (descriptive: spending is
  endogenous to health needs).

## References

World Bank (2026). World Development Indicators. Washington, DC: The
World Bank. <https://datacatalog.worldbank.org/search/dataset/0037712>

World Health Organization (2025). Global Health Expenditure Database.
<https://apps.who.int/nha/database>

Getzen, T. E. (2000). Health care is an individual necessity and a
national luxury: applying multilevel decision models to the analysis of
health care expenditures. Journal of Health Economics, 19(2), 259-270.

Baltagi, B. H. and Moscone, F. (2010). Health care expenditure and
income in the OECD reconsidered: evidence from panel data. Economic
Modelling, 27(4), 804-811.

Hausman, J. A. and Taylor, W. E. (1981). Panel data and unobservable
individual effects. Econometrica, 49(6), 1377-1398.

## Examples

``` r
data(healthspend)
hs <- transform(healthspend, lche = log(chereal), lgdp = log(gdppc))
ix <- c("iso3c", "year")

## income elasticity of real health spending
fe <- panglm(lche ~ lgdp, data = hs, index = ix, model = "within")
re <- panglm(lche ~ lgdp, data = hs, index = ix, model = "random")
fe2 <- panglm(lche ~ lgdp, data = hs, index = ix, model = "within",
              effect = "twoways", vcov = "cluster")
cbind(within = coef(fe), random = coef(re)[-1], twoways = coef(fe2))
#>       within   random   twoways
#> lgdp 1.14926 1.135463 0.8396004
panhausman(fe, re)
#> 
#>  Hausman specification test (fixed vs. random effects)
#> 
#> chisq = 3.735, df = 1, p-value = 0.05329
#> alternative hypothesis: one model is inconsistent
#> 
panhausman(fe, method = "mundlak", data = hs)
#> 
#>  Mundlak specification test (fixed vs. random effects), vcov: cluster
#> 
#> F = 1.463, df = 1, 185, p-value = 0.228
#> alternative hypothesis: one model is inconsistent
#> 

## persistence and long-run elasticity: system GMM, collapsed instruments
gs <- pangmm(lche ~ lag(lche, 1) + lgdp | lag(lche, 2:99) + lag(lgdp, 2:99),
             data = hs, index = ix, transformation = "ld",
             model = "twostep", collapse = TRUE)
summary(gs)
#> 
#> Call:
#> pangmm(formula = lche ~ lag(lche, 1) + lgdp | lag(lche, 2:99) + 
#>     lag(lgdp, 2:99), data = hs, index = ix, transformation = "ld", 
#>     model = "twostep", collapse = TRUE)
#> 
#> Estimator: two-way two-step system GMM (Blundell-Bond) 
#> n = 7886  groups = 193  T = 23  instruments = 66 
#> vcov: robust (Windmeijer) 
#> 
#> Coefficients:
#>              Estimate Std. Error z value Pr(>|z|)    
#> lag(lche, 1)  0.79999    0.06043  13.238   <2e-16 ***
#> lgdp          0.15476    0.09996   1.548    0.122    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Arellano-Bond AR(1) test: z = -11.51, p = < 2.2e-16
#> Arellano-Bond AR(2) test: z = 0.418, p = 0.676
#> Sargan test (not robust): chi2(42) = 161.4, p = 6.752e-16
#> Hansen J test (robust): chi2(42) = 65.94, p = 0.0106
#> Wald test for coefficients: chi2(2) = 365, p = < 2.2e-16
#> 
b <- coef(gs)
b[["lgdp"]] / (1 - b[["lag(lche, 1)"]])   # long-run income elasticity
#> [1] 0.7737543

## Hausman-Taylor: income endogenous; ageing and region exogenous
ht <- paniv(lche ~ lgdp + pop65 + region | pop65 + region,
            data = hs, index = ix, model = "hausmantaylor")
coef(ht)
#>                                             (Intercept) 
#>                                            -3.380823494 
#>                                                    lgdp 
#>                                             1.019289147 
#>                                                   pop65 
#>                                             0.043322428 
#>                             regionEurope & Central Asia 
#>                                            -0.080872877 
#>                         regionLatin America & Caribbean 
#>                                             0.047087247 
#> regionMiddle East, North Africa, Afghanistan & Pakistan 
#>                                            -0.051324970 
#>                                     regionNorth America 
#>                                             0.477931406 
#>                                        regionSouth Asia 
#>                                            -0.245957148 
#>                                regionSub-Saharan Africa 
#>                                             0.002156185 

## out-of-pocket share, a proportion: beta regression
hs$oopshare <- hs$oop / 100
bp <- panglm(oopshare ~ lgdp + pop65, data = hs, index = ix,
             model = "pooling", family = "beta", vcov = "cluster")
coef(bp)
#>   (Intercept)          lgdp         pop65 
#>  1.041328e+00 -1.844259e-01 -3.140689e-05 

## quantile heterogeneity of the elasticity
coef(panqr(lche ~ lgdp, data = hs, index = ix, tau = c(0.1, 0.5, 0.9),
           model = "mmqr"))
#>       tau=0.1  tau=0.5  tau=0.9
#> lgdp 1.179659 1.149445 1.118256

# \donttest{
## random-effects beta with Mundlak terms, every fifth year (slower)
hs5 <- subset(hs, year %% 5 == 0)
br <- panglm(oopshare ~ lgdp + pop65, data = hs5, index = ix,
             model = "random", family = "beta", cre = TRUE,
             vcov = "cluster", R = 7)
summary(br)
#> 
#> Call:
#> panglm(formula = oopshare ~ lgdp + pop65, data = hs5, index = ix, 
#>     model = "random", family = "beta", vcov = "cluster", cre = TRUE, 
#>     R = 7)
#> 
#> Family: beta  Link: logit  Model: random  Effect: individual 
#> n = 916  groups = 186  vcov: cluster 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.77971    0.68937   2.582 0.009833 ** 
#> lgdp        -0.19501    0.10056  -1.939 0.052476 .  
#> pop65       -0.04780    0.01323  -3.613 0.000302 ***
#> mean_lgdp   -0.08168    0.12642  -0.646 0.518234    
#> mean_pop65   0.04828    0.01724   2.801 0.005101 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: 911.1 
#> Dispersion parameter: 0.02356 
#> Random-intercept SD: 0.9156 
#> Mundlak test (individual means jointly zero): chi2(2) = 8.734, p = 0.01269
#> Residual df: 909  Iterations: 84
#> 

## country-specific income slopes (centre the regressor first)
hs$lgdpc <- hs$lgdp - 9
rs <- panglm(lche ~ lgdpc, data = hs, index = ix, model = "random",
             random = ~ lgdpc, R = 7)
summary(rs)
#> 
#> Call:
#> panglm(formula = lche ~ lgdpc, data = hs, index = ix, model = "random", 
#>     random = ~lgdpc, R = 7)
#> 
#> Family: gaussian  Link: identity  Model: random  Effect: individual 
#> n = 4222  groups = 186  vcov: classical 
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  6.06809    0.06365   95.34   <2e-16 ***
#> lgdpc        1.18760    0.04754   24.98   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: 1144 
#> Dispersion parameter: 0.02321 
#> Random effects (SD): (Intercept) 0.8076, lgdpc 0.5966 
#> Random-effects correlation: lgdpc ~ (Intercept) -0.04519 
#> Residual df: 4220  Iterations: 12
#> 
# }
```
