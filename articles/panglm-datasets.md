# Health-economics datasets for teaching panel models

`panglm` ships one simulated and five real panel datasets so that every
estimator family can be taught on data where it makes substantive sense.
The simulated `copd` panel has known data-generating values, so students
can first check that an estimator recovers the truth; the real panels
then pose the same questions on public data, and their build scripts
(`data-raw/` in the GitHub repository) reproduce published numbers to
show that the data were assembled correctly. Each help page (for example
[`?docvisits`](https://elbadisyimad.com/panglm/reference/docvisits.md))
lists every variable, the source and licence, the validation, and a
“Teaching use” section with research questions.

``` r

library(panglm)
```

## Overview

| Dataset | Panel | Outcome types | Estimators | Research question |
|----|----|----|----|----|
| `copd` (simulated) | 60 patients x 5 visits | continuous, count, binary, cost with zeros, proportion, ordinal | all families, [`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md), [`panqr()`](https://elbadisyimad.com/panglm/reference/panqr.md), [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md), [`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md), [`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md) | Do the estimators recover the known effects, including that of a staggered care-management programme on cost? |
| `docvisits` | 7,293 persons x 1984-1994 (7 waves, unbalanced) | count, binary, ordinal (0-10) | Poisson / negbin within and random, [`panhausman()`](https://elbadisyimad.com/panglm/reference/panhausman.md), probit with `dynamic = TRUE`, FE ordered logit, [`panfmm()`](https://elbadisyimad.com/panglm/reference/panfmm.md), [`panhurdle()`](https://elbadisyimad.com/panglm/reference/panhurdle.md), [`paniv()`](https://elbadisyimad.com/panglm/reference/paniv.md) Hausman-Taylor | Does add-on insurance raise doctor visits; is doctor use state dependent? |
| `meps` | 3,000 US adults x 2018-2019 | expenditure with zeros, counts, binary, ordinal | [`twopart()`](https://elbadisyimad.com/panglm/reference/twopart.md), `family = "gamma"`, `"tobit"`, [`panqr()`](https://elbadisyimad.com/panglm/reference/panqr.md), [`panselect()`](https://elbadisyimad.com/panglm/reference/panselect.md), negbin, `cre = TRUE` | How much less do the uninsured spend, and how much of the gap is selection? |
| `medicaid` | 51 states x 2009-2019 | rates (continuous) | [`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md), [`etwfe()`](https://elbadisyimad.com/panglm/reference/etwfe.md), [`eventstudy()`](https://elbadisyimad.com/panglm/reference/eventstudy.md), [`did()`](https://elbadisyimad.com/panglm/reference/did.md), two-way within | Did the ACA Medicaid expansion reduce uninsurance? |
| `healthspend` | 193 countries x 2000-2022 | continuous (log spending), shares in (0, 1) | within / random / CRE, [`panhausman()`](https://elbadisyimad.com/panglm/reference/panhausman.md), [`pangmm()`](https://elbadisyimad.com/panglm/reference/pangmm.md), [`paniv()`](https://elbadisyimad.com/panglm/reference/paniv.md) Hausman-Taylor, `family = "beta"`, random slopes, [`panqr()`](https://elbadisyimad.com/panglm/reference/panqr.md) | Is health care a luxury good? |
| `cigar` | 46 US states x 1963-1992 | continuous (log demand) | pooled vs within vs [`pangmm()`](https://elbadisyimad.com/panglm/reference/pangmm.md), `dynamic = TRUE` | How addictive is cigarette demand, and what is the long-run price elasticity? |

Each exercise below states a question, fits the model and interprets the
numbers it prints.

## copd: recovering a known programme effect

**Question.** Clinics 1 and 2 adopted a care-management programme at
visit 3, clinics 3 and 4 at visit 4, and clinics 5 and 6 never did. By
construction, the programme cuts the cost of users by `1 - exp(-0.5)`,
about 39 percent, and its sample average effect on the treated visits is
-406 euros (see
[`?copd`](https://elbadisyimad.com/panglm/reference/copd.md)). Do the
estimators find it?

``` r

idx <- c("id", "visit")
att <- attgt(cost ~ 1, data = copd, index = idx, treat = "adopt",
             control = "notyet")
summary(att)
#> 
#> Callaway and Sant'Anna group-time ATT
#> Aggregation: overall (simple weighted average)
#> Control group: not yet treated  Estimator: unconditional  Individuals: 60 
#> 
#> ATT: -184 (SE 199.1, 95% CI -574.2 to 206.2)
#> 
#> Pre-trend test (pre-treatment ATT(g, t) = 0): Chi-squared = 2.116, df = 3, p-value = 0.5487

tp <- twopart(cost ~ exacerbations + hospitalized + smoker + treatment + crp +
                visit + programme + factor(clinic),
              binary_formula = ~ exacerbations + hospitalized + smoker +
                treatment + crp,
              data = copd, index = idx, positive = "lognormal")
coef(tp$positive)[["programme"]]
#> [1] -0.5178754
```

**Interpretation.** With only 60 patients the Callaway and Sant’Anna
estimate is imprecise (about -184 euros, standard error about 200): its
confidence interval covers the true -406 euros, and the pre-trend test
does not reject parallel trends, but a single small panel cannot pin the
effect down. The two-part model is far more precise because it
conditions on exacerbations and admissions, which explain most of the
cost variation: its programme coefficient on the log positive cost is
about -0.52 (standard error 0.13) against a true value of -0.5.
Rerunning the design on a large simulated panel (as `data-raw/copd.R`
does) brings both estimators close to the truth; that is the point of a
known-truth teaching dataset.

## docvisits: state dependence in doctor use

**Question.** People who saw a doctor last year are more likely to see
one this year. Is that habit (true state dependence) or persistent ill
health (unobserved heterogeneity)? We use the consecutive waves 1984 to
1988 and the Wooldridge initial-conditions approach (`dynamic = TRUE`).

``` r

early <- subset(docvisits, year <= 1988)
f <- doctor ~ age + hhninc + hhkids + married + working
pooled <- panglm(f, data = early, index = c("id", "year"), model = "pooling",
                 family = "probit", dynamic = TRUE, vcov = "cluster")
random <- panglm(f, data = early, index = c("id", "year"), model = "random",
                 family = "probit", dynamic = TRUE, R = 8)
cbind(pooled = coef(pooled)[c("lag_doctor", "initial_doctor")],
      random = coef(random)[c("lag_doctor", "initial_doctor")])
#>                   pooled    random
#> lag_doctor     0.6453130 0.2824974
#> initial_doctor 0.5457217 0.9103464
```

**Interpretation.** The pooled dynamic probit gives a lagged-visit
coefficient of about 0.65, which mixes true state dependence with
persistent unobserved health. The random-effects probit with initial
conditions shrinks it to about 0.28 (roughly 8 percentage points on the
probability scale, see
[`partialeffects()`](https://elbadisyimad.com/panglm/reference/partialeffects.md)),
while the initial-condition coefficient and the random-effect variance
absorb the rest: most of the apparent persistence is heterogeneity, not
habit.

## meps: the spending gap of the uninsured

**Question.** How much less do uninsured adults spend on health care,
and how much of the cross-sectional gap survives once we compare the
same person before and after a change of coverage?

``` r

tp <- twopart(totexp ~ insurance + I(age / 10) + sex + nchronic + povcat,
              data = meps, index = c("id", "year"), positive = "gamma")
partialeffects(tp, variables = "insurance")
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

pos <- subset(meps, totexp > 0)
pooled <- panglm(totexp ~ insurance + povcat + I(age / 10) + sex + nchronic,
                 data = pos, index = c("id", "year"), family = "gamma",
                 model = "pooling")
within <- panglm(totexp ~ insurance + povcat, data = pos,
                 index = c("id", "year"), family = "gamma", model = "within")
exp(c(pooled = coef(pooled)[["insuranceuninsured"]],
      within = coef(within)[["insuranceuninsured"]])) - 1
#>     pooled     within 
#> -0.6795660 -0.3096233
```

**Interpretation.** The two-part model puts the average partial effect
of being uninsured rather than privately insured at about -6,200 dollars
per person-year, while public-only coverage differs little from private.
Among people with positive spending, the pooled gamma model says the
uninsured spend about 68 percent less; the within estimate, identified
only from people who gained or lost coverage between 2018 and 2019,
shrinks this to about 31 percent less. Much of the cross-sectional gap
reflects who is uninsured rather than the effect of losing coverage.

## medicaid: staggered adoption of the Medicaid expansion

**Question.** Did the Affordable Care Act Medicaid expansion reduce
uninsurance among adults below the poverty line? Adoption was staggered
(2014 to 2019 within the panel; later adopters and the ten non-adopters
serve as controls), so the single-coefficient two-way fixed effects
regression may compare late adopters with already-treated early
adopters.
[`attgt()`](https://elbadisyimad.com/panglm/reference/attgt.md) builds
clean comparisons with never-treated states.

``` r

ix <- c("state", "year")
cs <- attgt(uninspov ~ 1, data = medicaid, index = ix, treat = "adopt",
            control = "never")
summary(cs, type = "dynamic")
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

tw <- did(uninspov ~ 1, data = medicaid, index = ix, treat = "adopt")
ew <- etwfe(uninspov ~ 1, data = medicaid, index = ix, treat = "adopt")
c(attgt = summary(cs)$overall$estimate, etwfe = summary(ew)$overall$estimate,
  twfe = coef(tw)[["treated"]])
#>     attgt     etwfe      twfe 
#> -6.539795 -7.116279 -7.206276

# first stage: Medicaid coverage of the same adults
summary(attgt(medicaidpov ~ 1, data = medicaid, index = ix,
              treat = "adopt", control = "never"))$overall
#>   estimate      se    lower   upper
#> 1  10.2355 1.40579 7.480201 12.9908
```

**Interpretation.** Leads before expansion are small and scattered
around zero, which supports parallel trends. The uninsured share of
adults below the poverty line falls by about 5 percentage points in the
year of expansion and by 6 to 8 points after two to three years; the
overall group-time ATT is about -6.5 points (SE about 1.3). Medicaid
coverage of the same adults rises by about 10 points, more than
uninsurance falls: part of the gain replaces other coverage (crowd-out).
Here the two-way fixed effects estimate (about -7.2) stays close to the
robust ones (extended TWFE about -7.1), because most treated states
adopted in 2014 and the effect is fairly flat after two years; with more
spread-out adoption and stronger dynamics the gap can be large.

## healthspend: is health care a luxury?

**Question.** Does health spending rise more than proportionally with
national income (an income elasticity above 1)?

``` r

hs <- transform(healthspend, lche = log(chereal), lgdp = log(gdppc))
ix <- c("iso3c", "year")
fe <- panglm(lche ~ lgdp, data = hs, index = ix, model = "within")
re <- panglm(lche ~ lgdp, data = hs, index = ix, model = "random")
fe2 <- panglm(lche ~ lgdp, data = hs, index = ix, model = "within",
              effect = "twoways", vcov = "cluster")
cbind(within = coef(fe), random = coef(re)[-1], twoways = coef(fe2))
#>       within   random   twoways
#> lgdp 1.14926 1.135463 0.8396004
panhausman(fe, method = "mundlak", data = hs)
#> 
#>  Mundlak specification test (fixed vs. random effects), vcov: cluster
#> 
#> F = 1.463, df = 1, 185, p-value = 0.228
#> Wald chisq = 1.463, df = 1, p-value = 0.2264
#> alternative hypothesis: one model is inconsistent
```

**Interpretation.** Within countries, real health spending grows by
about 1.15 percent for each 1 percent of income, which looks like a
luxury good, and the Mundlak test does not reject random effects. Adding
year effects brings the elasticity down to about 0.84: much of the
apparent luxury effect is a common upward trend in spending across all
countries.

## cigar: habits and the long-run price elasticity

**Question.** Cigarette demand is addictive: how strong is the habit,
and what is the long-run price elasticity? The lagged dependent variable
makes pooled OLS biased upwards and the within estimator biased
downwards (Nickell bias); difference GMM should lie in between.

``` r

ix <- c("state", "year")
f <- lnc ~ lag(lnc, 1) + lnp + lnpn + lny
ols <- panglm(f, data = cigar, index = ix, model = "pooling")
fe <- panglm(f, data = cigar, index = ix, model = "within", effect = "twoways")
ab <- pangmm(lnc ~ lag(lnc, 1) + lnp + lnpn + lny | lag(lnc, 2:99),
             data = cigar, index = ix, transformation = "d", model = "onestep")
#> Warning: 437 instruments for 46 individuals: too many instruments bias GMM
#> toward OLS and weaken the Hansen test; use collapse = TRUE or a smaller maxlag
tab <- rbind(OLS = coef(ols)[-1], Within = coef(fe), GMM = coef(ab)[1:4])
round(tab, 3)
#>        lag(lnc, 1)    lnp   lnpn    lny
#> OLS          0.969 -0.090  0.024 -0.031
#> Within       0.833 -0.299  0.034  0.100
#> GMM          0.843 -0.377 -0.016  0.139
round(tab[, "lnp"] / (1 - tab[, 1]), 2)  # long-run price elasticity
#>    OLS Within    GMM 
#>  -2.96  -1.79  -2.40
```

**Interpretation.** The habit coefficient from difference GMM (0.84)
lies between the within (0.83) and OLS (0.97) estimates, as theory
predicts, and the three rows reproduce Table 8.1 of Baltagi’s
*Econometric Analysis of Panel Data*. The short-run price elasticity is
about -0.38 and the long-run one about -2.4.
